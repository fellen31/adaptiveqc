process SAMTOOLS_VIEW_TARGET_READS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(bed)

    output:
    tuple val(meta), path("*.txt.gz"),  emit: read_ids
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def target_file = bed ? "--target-file ${bed}" : ""

    // Could just add a small r-script that makes the txt to parquet conversion here, and then set appropriate hive-structure in modules.config...
    """
    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        ${target_file} \\
        ${reference} \\
        $args \\
        $input \\
        $args2 | awk '{print \$1}' | gzip -c > ${meta.id}.read_id.txt.gz

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
