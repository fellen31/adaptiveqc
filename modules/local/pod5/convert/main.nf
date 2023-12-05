process POD5_CONVERT {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/fellen31/pod5:0.3.2"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "POD5_CONVERT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.pod5"), emit: pod5
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // TODO: {fast5,from_fast5,to_fast5} in args
    """
    pod5 convert fast5 \\
        $args \\
        -o ${prefix}.pod5 \\
        -t ${task.cpus} \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pod5: \$(pod5 -v 2>&1 | sed -n 's/^Pod5 version: //g')
    END_VERSIONS
    """
}
