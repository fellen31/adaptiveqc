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

include { FQCRS                             } from '../modules/local/fqcrs/main'
include { FQCRS as FQCRS_REBASECALL         } from '../modules/local/fqcrs/main'
include { DORADO_BASECALLER                 } from '../modules/local/dorado/basecaller/main'
include { DORADO_DEMUX                      } from '../modules/local/dorado/demux/main'
include { POD5_CONVERT                      } from '../modules/local/pod5/convert/main'
include { ONT_SEQUENCING_SUMMARY_TO_PARQUET } from '../modules/local/ont_sequencing_summary_to_parquet/main'
include { FQCRS_TO_HIVE as FQCRS_TO_HIVE_ORIGINAL   } from '../modules/local/fqcrs_to_hive/main'
include { FQCRS_TO_HIVE as FQCRS_TO_HIVE_REBASECALL } from '../modules/local/fqcrs_to_hive/main'

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

include { SAMTOOLS_CAT                } from '../modules/nf-core/samtools/cat/main'
include { SAMTOOLS_SORT               } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_VIEW as SAMTOOLS_PASS } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX              } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FASTQ              } from '../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_SAMTOOLS_PASS } from '../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_REBASECALL } from '../modules/nf-core/samtools/fastq/main'
include { MOSDEPTH as MOSDEPTH_REBASECALL } from '../modules/nf-core/mosdepth/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

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

    /* Main workflow

       1. Do rebasecalling first if params.basecalling
       2. Mix channels if not rebasecalling

    */

    ch_versions = Channel.empty()

    // Chech mandatory input files
    ch_sample         = Channel.fromSamplesheet('input', immutable_meta: false)
    ch_fasta          = Channel.fromPath(params.fasta).map { it -> [it.simpleName, it] }.collect()
    ch_bed            = Channel.fromPath(params.bed).map { it -> [it.simpleName, it] }.collect()

    modified_paths = ch_sample.map{ sample, dir -> "${dir}**" }
    all_files_channel = modified_paths.flatMap { path -> file(path).collect() }

    ch_mod_bases = []
    ch_kit = []
    if(params.barcode_kit) { ch_kit = params.barcode_kit }


    // Index genome
    if(!params.skip_mapping) {

        PREPARE_GENOME( ch_fasta )
        ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

        // Gather indices
        fasta = ch_fasta
        fai   = PREPARE_GENOME.out.fai
        mmi   = PREPARE_GENOME.out.mmi
    }

    // From here starts "new" (from ONT samplesheet)
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

    // Extract the base dir for experiment - based on samplesheet location
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

    // Then we can
    base_dir
        .combine(all_files_branched.fast5)
        // Keep only the correct base_dir match (sample_sheet) for the run
        // by checking path for run_id. Almost like a bad join
        .filter{ sample_sheet_meta, dir, name, path ->
            path =~ /.*${sample_sheet_meta.run_id}.*/
        }
        .set{ ch_pod5_convert_in }

    ch_pod5_convert_in
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
        // Think this might work?
        .set{ ch_pod5_convert_in }

    // Convert fast5 files to pod5
    POD5_CONVERT( ch_pod5_convert_in )
    // This then needs to be joined with the other pod5's
    POD5_CONVERT.out.pod5
        .map{ meta, path ->
            converted_meta = [
                id:meta.id,
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

    // Ok, redo fastq here base_dir
    base_dir
        .combine(all_files_branched.fastq)
        // Keep only the correct base_dir match (sample_sheet) for the run
        // by checking path for run_id. Almost like a bad join
        .filter{ sample_sheet_meta, dir, name, path ->
            path =~ /.*${sample_sheet_meta.run_id}.*/
        }
        // Problem if we have the 'barcode' in name maybe... TODO: safest to specify in samplesheet?
        // Yes, especially if we can't use name and have to use path
        .branch{ sample_sheet_meta, dir, name, path ->
            barcoded: path =~ /barcode|unclassified|mixed/         // We can't make this better either because what if we don't
            non_barcoded: !(path =~ /barcode|unclassified|mixed/)  // Divide the BAMs into pass and fail?
                                                                   // TODO: Add check for pass or fail division or not
        }
        .set{ new_fastq_channel }

    new_fastq_channel.barcoded
          // Should also test that sample_sheet_meta.experiment == experiment from path
          .map{ new_fastq_sample_sheet_meta, new_fastq_dir, new_fastq_name, new_fastq_path ->
            // This will be [fastq_fail, barcode19] if two dirs exist, can this be contiditonal somehow, that barcode = '' if
            // the barcode dir does not exist?
            // Not sure where or when unique naming is needed but somewhere it messes up if I don't have it
            (new_fastq_type, new_fastq_barcode) = new_fastq_dir.toString()
                .replaceAll(new_fastq_sample_sheet_meta.base_dir, '')
                .replaceAll(new_fastq_sample_sheet_meta.experiment_dir, '')
                .tokenize("/")
            new_fastq_meta = [
                id:new_fastq_name,
                base_dir:new_fastq_sample_sheet_meta.base_dir,
                experiment_dir:new_fastq_sample_sheet_meta.experiment_dir,
                experiment:new_fastq_sample_sheet_meta.experiment,
                sample:new_fastq_sample_sheet_meta.sample,
                run_id:new_fastq_sample_sheet_meta.run_id,
                type:new_fastq_type,
                barcode:new_fastq_barcode
            ]
            [new_fastq_meta, new_fastq_path]
          }
        // Think this might work?
        .set{ new_fastq_channel_barcoded_meta }

// new_fastq_channel_barcoded_meta.view() // why is this just sometimes working...

    // Can we just do this?
    new_fastq_channel.non_barcoded
          // Should also test that sample_sheet_meta.experiment == experiment from path
          .map{ new_fastq_non_barcoded_sample_sheet_meta, new_fastq_non_barcoded_dir, new_fastq_non_barcoded_name, new_fastq_non_barcoded_path ->
            // This will be [fastq_fail, barcode19] if two dirs exist, can this be contiditonal somehow, that barcode = '' if
            // the barcode dir does not exist?
            // Not sure where or when unique naming is needed but somewhere it messes up if I don't have it
            (new_fastq_non_barcoded_type, new_fastq_non_barcoded_barcode) = new_fastq_non_barcoded_dir.toString()
                .replaceAll(new_fastq_non_barcoded_sample_sheet_meta.base_dir, '')
                .replaceAll(new_fastq_non_barcoded_sample_sheet_meta.experiment_dir, '')
                .tokenize("/")
            new_fastq_non_barcoded_meta = [
                id:new_fastq_non_barcoded_name,
                base_dir:new_fastq_non_barcoded_sample_sheet_meta.base_dir,
                experiment_dir:new_fastq_non_barcoded_sample_sheet_meta.experiment_dir,
                experiment:new_fastq_non_barcoded_sample_sheet_meta.experiment,
                sample:new_fastq_non_barcoded_sample_sheet_meta.sample,
                run_id:new_fastq_non_barcoded_sample_sheet_meta.run_id,
                type:new_fastq_non_barcoded_type,
                barcode:'non_barcoded'
            ]
            [new_fastq_non_barcoded_meta, new_fastq_non_barcoded_path]
          }
        // Think this might work?
        .set{ new_fastq_non_barcoded_channel_non_barcoded_meta }

        //new_fastq_non_barcoded_channel_non_barcoded_meta.view()
final_fastq_channel = new_fastq_channel_barcoded_meta.concat(new_fastq_non_barcoded_channel_non_barcoded_meta)

    // Try this for BAM files
    base_dir
        .combine(all_files_branched.bam)
        // Keep only the correct base_dir match (sample_sheet) for the run
        // by checking path for run_id. Almost like a bad join
        .filter{ sample_sheet_meta, dir, name, path ->
            path =~ /.*${sample_sheet_meta.run_id}.*/
        }
        // Problem if we have the 'barcode' in name maybe... TODO: safest to specify in samplesheet?
        .branch{ sample_sheet_meta, dir, name, path ->
            barcoded: name =~ /barcode|unclassified|mixed/         // We can't make this better either because what if we don't
            non_barcoded: !(name =~ /barcode|unclassified|mixed/)  // Divide the BAMs into pass and fail?
                                                                   // TODO: Add check for pass or fail division or not
        }
        .set{ new_bam_channel }

    new_bam_channel.barcoded
          // Should also test that sample_sheet_meta.experiment == experiment from path
          .map{ new_bam_barcoded_sample_sheet_meta, new_bam_barcoded_dir, new_bam_barcoded_name, new_bam_barcoded_path ->

            (new_bam_barcoded_type, new_bam_barcoded_barcode) = new_bam_barcoded_dir.toString()
                .replaceAll(new_bam_barcoded_sample_sheet_meta.base_dir, '')
                .replaceAll(new_bam_barcoded_sample_sheet_meta.experiment_dir, '')
                .tokenize("/")
            new_bam_barcoded_meta = [
                id:new_bam_barcoded_name,
                base_dir:new_bam_barcoded_sample_sheet_meta.base_dir,
                experiment_dir:new_bam_barcoded_sample_sheet_meta.experiment_dir,
                experiment:new_bam_barcoded_sample_sheet_meta.experiment,
                sample:new_bam_barcoded_sample_sheet_meta.sample,
                run_id:new_bam_barcoded_sample_sheet_meta.run_id,
                type:new_bam_barcoded_type,
                barcode:new_bam_barcoded_barcode
            ]
            [new_bam_barcoded_meta, new_bam_barcoded_path]
          }
        // Think this might work?
        .set{ new_bam_channel_barcoded_meta }

    // Can we just do this?
    new_bam_channel.non_barcoded
          // Should also test that sample_sheet_meta.experiment == experiment from path
          .map{ new_bam_non_barcoded_sample_sheet_meta, new_bam_non_barcoded_dir, new_bam_non_barcoded_name, new_bam_non_barcoded_path ->
            // This will be [fastq_fail, barcode19] if two dirs exist, can this be contiditonal somehow, that barcode = '' if
            // the barcode dir does not exist?
            // Not sure where or when unique naming is needed but somewhere it messes up if I don't have it
            (new_bam_non_barcoded_type, new_bam_non_barcoded_barcode) = new_bam_non_barcoded_dir.toString()
                .replaceAll(new_bam_non_barcoded_sample_sheet_meta.base_dir, '')
                .replaceAll(new_bam_non_barcoded_sample_sheet_meta.experiment_dir, '')
                .tokenize("/")
            new_bam_non_barcoded_meta = [
                id:new_bam_non_barcoded_name,
                base_dir:new_bam_non_barcoded_sample_sheet_meta.base_dir,
                experiment_dir:new_bam_non_barcoded_sample_sheet_meta.experiment_dir,
                experiment:new_bam_non_barcoded_sample_sheet_meta.experiment,
                sample:new_bam_non_barcoded_sample_sheet_meta.sample,
                run_id:new_bam_non_barcoded_sample_sheet_meta.run_id,
                type:new_bam_non_barcoded_type,
                barcode:'non_barcoded'
            ]
            [new_bam_non_barcoded_meta, new_bam_non_barcoded_path]
          }
        // Think this might work?
        .set{ new_bam_channel_non_barcoded_meta }

    final_bam_channel = new_bam_channel_barcoded_meta.concat(new_bam_channel_non_barcoded_meta)
        // We can't groupTuple on items within meta, so we have to temporarily
        // take them out and refer to them by key
    final_bam_channel
        .map{ hej_meta, hej_path ->
            hej_group = [hej_meta.base_dir,
                         hej_meta.experiment_dir,
                         hej_meta.experiment,
                         hej_meta.sample,
                         hej_meta.run_id,
                         hej_meta.type,
                         hej_meta.barcode
                         ].join("f2c471d8-5cdc-4194-a723-ce2a39ce8683")

            [hej_group, hej_path]
        }
        .groupTuple(by: 0) // Now this is always 98, but relies on "." to join and split
        // After grouping, we can then put the meta back together, with a new id
        .map{ group, paths ->
            //id = [experiment, sample, type, barcode].join(".")
            (base_dir_i, experiment_dir_i, experiment_i, sample_i, run_id_i, type_i, barcode_i) = group.toString().split("f2c471d8-5cdc-4194-a723-ce2a39ce8683")

            new_meta_all = [
                id:[experiment_i, sample_i, type_i, barcode_i].join("."),
                base_dir:base_dir_i,
                experiment_dir:experiment_dir_i,
                experiment:experiment_i,
                sample:sample_i,
                run_id:run_id_i,
                type:type_i,
                barcode:barcode_i,
            ]
            [new_meta_all, paths]
        }
        .set{ all_bam_grouped }

    final_fastq_channel
        .map{ final_fastq_meta, final_fastq_path ->
            final_fastq_group = [final_fastq_meta.base_dir,
                         final_fastq_meta.experiment_dir,
                         final_fastq_meta.experiment,
                         final_fastq_meta.sample,
                         final_fastq_meta.run_id,
                         final_fastq_meta.type,
                         final_fastq_meta.barcode
                         ].join("f2c471d8-5cdc-4194-a723-ce2a39ce8683")

            [final_fastq_group, final_fastq_path]
        }
        .groupTuple(by: 0) // Now this is always 98, but relies on "." to join and split
        // After grouping, we can then put the meta back together, with a new id
        .map{ group_final_fastq, paths_final_fastq ->
            //id = [experiment, sample, type, barcode].join(".")
            (base_dir_g, experiment_dir_g, experiment_g, sample_g, run_id_g, type_g, barcode_g) = group_final_fastq.toString().split("f2c471d8-5cdc-4194-a723-ce2a39ce8683")

            new_meta_final_fastq = [
                id:[experiment_g, sample_g, type_g, barcode_g].join("."),
                base_dir:base_dir_g,
                experiment_dir:experiment_dir_g,
                experiment:experiment_g,
                sample:sample_g,
                run_id:run_id_g,
                type:type_g,
                barcode:barcode_g,
            ]
            [new_meta_final_fastq, paths_final_fastq]
        }
        .set{ all_fastq_grouped }


    if(params.skip_basecall) {

        if(params.bam) {
            SAMTOOLS_CAT( all_bam_grouped )
            SAMTOOLS_SORT( SAMTOOLS_CAT.out.bam )
            SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

            SAMTOOLS_SORT.out.bam
                .join(SAMTOOLS_INDEX.out.bai)
                .groupTuple()
                .set{ bam_bai_non_rebasecalled }

            SAMTOOLS_FASTQ( bam_bai_non_rebasecalled.map{ meta, bam, bai -> [meta, bam] } )

            // Want to do FQCRS on all reads, including failed reads
            FQCRS ( SAMTOOLS_FASTQ.out.fastq )
            // Resuming might be problematic when/if not overwriting?
            FQCRS_TO_HIVE_ORIGINAL ( FQCRS.out.res )

        } else if (params.fastq) {
            FQCRS ( all_fastq_grouped )
            FQCRS_TO_HIVE_ORIGINAL ( FQCRS.out.res )
        }

    }

    // Extract sequencing summaries -> to_parquet, this should be able to be done in parallel right?
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
        .set { ch_make_fqcrs_hive_in }


    ONT_SEQUENCING_SUMMARY_TO_PARQUET( ch_make_fqcrs_hive_in )

    // Do we even need to know about barcodes at all to just run rebasecall? No?
    base_dir.combine(all_files_branched.pod5)
        .filter{ sample_sheet_meta, dir, name, path ->
            path =~ /.*${sample_sheet_meta.run_id}.*/
        }
            // This should at least work for the pod5's
            // I see no better way, unless we specify each and every file in the samplesheet
        .map{ sample_sheet_meta, dir, name, path ->
            group = [sample_sheet_meta.run_id, sample_sheet_meta.experiment, sample_sheet_meta.sample, sample_sheet_meta.run_id].join("f2c471d8-5cdc-4194-a723-ce2a39ce8683")
            [group, path]
        }
        .groupTuple()
        .map{ group, paths ->
            (id, experiment, sample, run_id ) = group.split("f2c471d8-5cdc-4194-a723-ce2a39ce8683")
            meta = [
               id:id,
               experiment:experiment,
               sample:sample,
               run_id:run_id,
            ]
            [meta, paths]
        }
        // For fastq and bam we might need other strategies when it comes to identifying barcodes
        // Can we look for 'barcode' or 'unclassified' in path?
        // Branch on this
        // And if we don't find any barcode set 'barcode':'no_barcode'?
        // Add pod5 coverted files
        .concat(ch_converted_pod5)
        .groupTuple()
        .set{ ch_dorado_basecall_in }


    if(params.mod_bases) { ch_mod_bases = params.mod_bases}

    if(!params.skip_basecall) {

        DORADO_BASECALLER( ch_dorado_basecall_in, params.dorado_model, ch_mod_bases)

        if(!params.skip_demux) {

            DORADO_DEMUX ( DORADO_BASECALLER.out.ubam, params.barcode_kit)

            DORADO_DEMUX.out.barcode_bam // This process will always output multiple files..
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
        }
        if(params.skip_demux) {

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

    SAMTOOLS_PASS( ch_samtools_fastq_rebasecall_in.map{meta, bam -> [meta, bam, []]}, ch_fasta, [] )

    // This should be combined into fcqrs-process
    SAMTOOLS_FASTQ_REBASECALL( ch_samtools_fastq_rebasecall_in )

    SAMTOOLS_FASTQ_SAMTOOLS_PASS( SAMTOOLS_PASS.out.bam )

    SAMTOOLS_FASTQ_REBASECALL.out.fastq
        .set{ ch_fqcrs_rebasecall_in }

    if(!params.skip_mapping) {
        ALIGN_READS( SAMTOOLS_FASTQ_SAMTOOLS_PASS.out.fastq, mmi )

        bam_bai = ALIGN_READS.out.bam_bai
        bam     = ALIGN_READS.out.bam
        bai     = ALIGN_READS.out.bai
    }

    FQCRS_REBASECALL( ch_fqcrs_rebasecall_in )
    FQCRS_TO_HIVE_REBASECALL( FQCRS_REBASECALL.out.res )

    }

    if(params.skip_basecall && !params.skip_mapping) {
        if(params.bam)
        {
        MOSDEPTH_REBASECALL( bam_bai_non_rebasecalled.combine(ch_bed.map{ meta, bed -> bed }), ch_fasta )
        }

    } else if (!params.skip_mapping){

        MOSDEPTH_REBASECALL( bam_bai.combine(ch_bed.map{ meta, bed -> bed }), ch_fasta )

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
