process SAMTOOLS_VIEW_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/fellen31/fqcrs_samtools:0.0.2"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "FQCRS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    path qname

    output:
    tuple val(meta), path("*.fastq.gz"),  emit: fastq
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def readnames = qname ? "--qname-file ${qname}": ""
    """
    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        ${reference} \\
        ${readnames} \\
        $args \\
        -O BAM \\
        $input |\\
    samtools \\
        fastq \\
        $args2 \\
        --threads ${task.cpus-1} \\
        | pigz -p ${task.cpus} > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
