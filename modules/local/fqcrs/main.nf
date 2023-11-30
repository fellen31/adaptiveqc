process FQCRS {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/fellen31/fqcrs:latest"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "FQCRS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(fastq)

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

    cat ${fastq} | fqcrs > res/${meta.experiment}/${meta.sample}/${meta.run_id}/${meta.type}/${meta.barcode}/${meta.id}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fqcrs: \$(echo "no-version" )
    END_VERSIONS
    """
}
