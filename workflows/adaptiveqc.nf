/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryMap; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowAdaptiveqc.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FQCRS_BAM                                       } from '../modules/local/fqcrs_bam/main'
include { DORADO_BASECALLER                               } from '../modules/local/dorado/basecaller/main'
include { DORADO_BASECALLER_DEMUX                         } from '../modules/local/dorado/basecaller_demux/main'
include { POD5_CONVERT                                    } from '../modules/local/pod5/convert/main'
include { ONT_SEQUENCING_SUMMARY_TO_PARQUET               } from '../modules/local/ont_sequencing_summary_to_parquet/main'
include { FQCRS_TO_HIVE as FQCRS_TO_HIVE_REBASECALL       } from '../modules/local/fqcrs_to_hive/main'
include { SAMTOOLS_VIEW_FASTQ as SAMTOOLS_VIEW_FASTQ_PASS } from '../modules/local/samtools/view_fastq/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { ALIGN_READS                } from '../subworkflows/local/align_reads'
include { PREPARE_GENOME             } from '../subworkflows/local/prepare_genome'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_SAMTOOLS_PASS } from '../modules/nf-core/samtools/fastq/main'
include { MOSDEPTH as MOSDEPTH_REBASECALL                } from '../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_500                       } from '../modules/nf-core/mosdepth/main'
include { MULTIQC                                        } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                    } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ADAPTIVEQC {

    //TODO: Samplesheet header validation
    // Run demux if specified in samplesheet

    ch_versions = Channel.empty()

    // Chech mandatory input files
    ch_sample         = Channel.fromSamplesheet('input', immutable_meta: false)
    ch_fasta          = Channel.fromPath(params.fasta).map { it -> [it.simpleName, it] }.collect()
    ch_bed            = Channel.fromPath(params.bed).map { it -> [it.simpleName, it] }.collect()

    modified_paths = ch_sample.map{ sample, dir -> "${dir}**" }
    all_files_channel = modified_paths.flatMap { path -> file(path).collect() }

    ch_mod_bases = []
    ch_kit = []

    if(!params.skip_demux) {
        if(params.barcode_kit) { ch_kit = params.barcode_kit } // todo, exit on not set
    }

    if(params.mod_bases) { ch_mod_bases = params.mod_bases} // todo, exit on not set

    // Index genome
    if(!params.skip_mapping) {

        PREPARE_GENOME( ch_fasta )
        ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

        // Gather indices
        fasta = ch_fasta
        fai   = PREPARE_GENOME.out.fai
        mmi   = PREPARE_GENOME.out.mmi
    }

    // Start by putting all files into channels
    all_files_channel
        .map { it -> [ it.getParent(), it.getName(), it] }
        .branch{ dir, name, path ->
                 sample_sheet: name =~ /^sample_sheet_.*\.csv$/
                 pod5: path.getExtension() == "pod5"
                 bam: path.getExtension() == "bam"
                 fastq: name =~  /.*\.fastq\.gz$/ // ok?
                 fast5: path.getExtension() == "fast5"
        }
        .set { all_files_branched }

    // Extract the base dir for experiment - based on samplesheet location, might change in updated MinKNOW version
    all_files_branched.sample_sheet
        .map{ dir, name, path ->
            run_id = dir.toString().tokenize("/").takeRight(3).join("/")
            base_dir = dir.toString().replaceAll(run_id, '')
            experiment_dir = dir.toString().replaceAll(base_dir, '')
            (experiment, sample, run_id) = experiment_dir.tokenize("/")
            new_meta = [
                base_dir:base_dir,
                experiment_dir:experiment_dir,
                experiment:experiment,
                sample:sample,
                run_id:run_id,
            ]
            [new_meta]
        }
        .set{ base_dir }

    // Then we can look for fast5-files that needs to be converted to pod5
    base_dir
        .combine(all_files_branched.fast5)
        // Keep only the correct base_dir match (sample_sheet) for the run
        // by checking path for run_id. Almost like a bad join
        .filter{ sample_sheet_meta, dir, name, path ->
            path =~ /.*${sample_sheet_meta.run_id}.*/
        }
        .map{ new_fast5_sample_sheet_meta, new_fast5_dir, new_fast5_name, new_fast5_path ->
            new_fast5_meta =[
                id:new_fast5_name,
                base_dir:new_fast5_sample_sheet_meta.base_dir,
                experiment_dir:new_fast5_sample_sheet_meta.experiment_dir,
                experiment:new_fast5_sample_sheet_meta.experiment,
                sample:new_fast5_sample_sheet_meta.sample,
                run_id:new_fast5_sample_sheet_meta.run_id,
            ]
            [new_fast5_meta, new_fast5_path]
          }
        .set{ ch_pod5_convert_in }

    // Convert fast5 files to pod5

    POD5_CONVERT( ch_pod5_convert_in )

    // This then needs to be joined with the other pod5's,
    POD5_CONVERT.out.pod5
        .map{ meta, path ->
            converted_meta = [
                //id:meta.id,
                //base_dir:meta.base_dir,
                //experiment_dir:meta.experiment_dir,
                experiment:meta.experiment,
                sample:meta.sample,
                run_id:meta.run_id//,
                //type:'fast5_converted'
            ]
        [converted_meta, path]
        }
        .set{ ch_converted_pod5 }

    // Extract sequencing summaries -> to_parquet,
    // Could split these files to reduce RAM-usage below 128 GB
    all_files_channel
        .map { it -> [ ['dir':it.getParent()], it] }
        .filter { meta, file -> file =~ /sequencing_summary_/ }
        .map { id, summary ->
            (experiment, sample, run_id) = id.dir.toString().tokenize("/").takeRight(3)
            id = [experiment, sample, run_id].join(".")
                meta = [
                    id:id,
                    experiment:experiment,
                    sample:sample,
                    run_id:run_id,
                ]
                [meta, summary]
        }
        .set { ch_make_sequencing_summary_hive_in }

    ONT_SEQUENCING_SUMMARY_TO_PARQUET( ch_make_sequencing_summary_hive_in )

    base_dir.combine(all_files_branched.pod5)
        .filter{ sample_sheet_meta, dir, name, path ->
            path =~ /.*${sample_sheet_meta.run_id}.*/
        }
            // This should at least work for the pod5's
            // I see no better way, unless we specify each and every file in the samplesheet
        .map{ sample_sheet_meta, dir, name, path ->
            group = [sample_sheet_meta.run_id, sample_sheet_meta.experiment, sample_sheet_meta.sample].join("f2c471d8-5cdc-4194-a723-ce2a39ce8683")
            [group, path]
        }
        .concat(
               ch_converted_pod5.map{ pod5_converted_meta, pod5_converted_path ->
                    pod5_group = [pod5_converted_meta.run_id, pod5_converted_meta.experiment, pod5_converted_meta.sample].join("f2c471d8-5cdc-4194-a723-ce2a39ce8683")
                    [pod5_group, pod5_converted_path]
                }
               )
        .groupTuple()
        .map{ group, paths ->
            (run_id, experiment, sample) = group.split("f2c471d8-5cdc-4194-a723-ce2a39ce8683")
            meta = [
               id:run_id,
               experiment:experiment,
               sample:sample,
               run_id:run_id,
            ]
            [meta, paths]
        }
        .set{ ch_dorado_basecall_in }

    if(!params.skip_demux) {

        // Directly demux basecalled BAM's to reduce temporary files
         DORADO_BASECALLER_DEMUX( ch_dorado_basecall_in, params.dorado_model, ch_mod_bases, params.barcode_kit)


        DORADO_BASECALLER_DEMUX.out.barcode_bam // This process will always output multiple files..
            .transpose(remainder:true) // Transpose works but have to wait for all to finish - nedd to get it to group by barcode
            .map{ meta_demux, bam_demux ->
                (run_id, barcode_demux, extension_demux) = bam_demux.toString().tokenize(".")
                meta_demux_new = [
                    id:[meta_demux.experiment, meta_demux.sample, barcode_demux].join("."),
                    experiment:meta_demux.experiment,
                    sample:meta_demux.sample,
                    run_id:meta_demux.run_id,
                    type:'rebasecalled_all',
                    barcode:barcode_demux,
                    extension:extension_demux,
                    file:bam_demux.getName()
                ]
                [meta_demux_new, bam_demux]
            }
        .set{ ch_samtools_fastq_rebasecall_in}

    } else {

        DORADO_BASECALLER( ch_dorado_basecall_in, params.dorado_model, ch_mod_bases)

        DORADO_BASECALLER.out.ubam
            .map{ meta, ubam ->
                meta_non_demux = [
                    id:meta.id,
                    experiment:meta.experiment,
                    sample:meta.sample,
                    run_id:meta.run_id,
                    type:'rebasecalled_all',
                    barcode:'non_barcoded',
                    extension:ubam.getExtension(),
                    file:ubam.getName()
                ]
                [meta_non_demux, ubam]
            }
        .set{ ch_samtools_fastq_rebasecall_in}
    }

    // Extract pass fastq based on qs-tag for input to other pipelines
    SAMTOOLS_VIEW_FASTQ_PASS( ch_samtools_fastq_rebasecall_in.map{meta, bam -> [meta, bam, []]}, ch_fasta, [] )

    //Run fqcrs
    FQCRS_BAM( ch_samtools_fastq_rebasecall_in )
    // Make fqcrs-hive
    FQCRS_TO_HIVE_REBASECALL( FQCRS_BAM.out.res )

    // Run mapping workflow
    if(!params.skip_mapping) {
        ALIGN_READS( SAMTOOLS_VIEW_FASTQ_PASS.out.fastq, mmi )

        bam_bai = ALIGN_READS.out.bam_bai
        bam     = ALIGN_READS.out.bam
        bai     = ALIGN_READS.out.bai

        // Run mosdepth per 500 bp and per bed-file
        MOSDEPTH_REBASECALL( bam_bai.combine(ch_bed.map{ meta, bed -> bed }), ch_fasta )
        MOSDEPTH_500( bam_bai.map{ meta, bam, bai -> [meta, bam, bai, []] }, ch_fasta )

    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowAdaptiveqc.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowAdaptiveqc.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
