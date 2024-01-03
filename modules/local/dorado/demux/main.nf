process DORADO_DEMUX {
    tag "$meta.id"
    label 'process_high' // Something like if modified bases - process high?

    container "docker.io/fellen31/dorado:0.5.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "DORADO_DEMUX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(reads)
    val(kit)

    output:
    tuple val(meta), path("*.bam"), emit: barcode_bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    // This is not necessary so should be in config, not here...?
    def kit_name = kit ? "--kit-name ${kit}" : ''

    """
    dorado demux \\
        ${reads} \\
        $kit_name \\
        --output-dir . \\
        --threads ${task.cpus}

    for file in *_barcode*.*;
    do
        echo \$file
        new_name=\$(echo \$file | sed 's/${kit}_/${meta.id}./g')
        mv \$file \$new_name
    done

    mv unclassified.bam ${meta.id}.unclassified.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado -vv 2>&1 | sed -n 's/^dorado: *\\(.*\\)/\\1/p')
    END_VERSIONS
    """
}
