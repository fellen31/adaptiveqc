process FQCRS_BAM {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/fellen31/fqcrs_samtools:0.0.2"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "FQCRS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("res/"), emit: res
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p res/${meta.experiment}/${meta.sample}/${meta.run_id}/${meta.type}/${meta.barcode}

    samtools fastq -@ ${task.cpus} ${bam} | fqcrs | pigz -p ${task.cpus} > res/${meta.experiment}/${meta.sample}/${meta.run_id}/${meta.type}/${meta.barcode}/${meta.id}.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fqcrs: \$(echo "no-version" )
    END_VERSIONS
    """
}
