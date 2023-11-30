process ONT_SEQUENCING_SUMMARY_TO_PARQUET {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/fellen31/ont_sequencing_summary_to_parquet:231117"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "sequencing_summary_to_parquet module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(sequencing_summary)

    output:
    tuple val(meta), path("*"), emit: hive
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    """
    ont_sequencing_summary_to_parquet.R -s ${sequencing_summary} --experiment ${meta.experiment} --run_id ${meta.run_id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ont_sequecing_summary_to_parquet: \$(echo "no-version" )
    END_VERSIONS
    """
}
