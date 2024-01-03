process DORADO_BASECALLER_DEMUX {
    tag "$meta.id"
    label 'process_gpu'
    label 'process_high' // Something like if modified bases - process high?

    container "docker.io/fellen31/dorado:0.5.1"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "DORADO_BASECALLER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(pod5)
    val(model)
    val(modified_bases)
    val(kit)

    output:
    tuple val(meta), path("*.bam"), emit: barcode_bam
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"

    def mod_bases = modified_bases ? "--modified-bases ${modified_bases}" : ''
    def kit_name = kit ? "--kit-name ${kit}" : ''

    """
    dorado basecaller \\
        /models/${model} \\
        ./ \\
        $mod_bases |\\
        dorado demux \\
        $kit_name \\
        --output-dir . \\
        --threads ${task.cpus}

    for file in *_barcode*.bam;
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
