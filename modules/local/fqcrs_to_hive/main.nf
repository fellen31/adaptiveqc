process FQCRS_TO_HIVE {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/fellen31/ont_sequencing_summary_to_parquet:231223"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "fqcrs_to_hive module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(fqcrs)

    output:
    tuple val(meta), path("*.parquet"), emit: hive
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"

    """
    fqcrs_to_hive.R \\
        --input ${fqcrs} \\
        --experiment ${meta.experiment} \\
        --sample ${meta.sample} \\
        --run_id ${meta.run_id} \\
        --dir ${meta.type} \\
        --barcode ${meta.barcode} \\
        --out .

    echo ${meta.experiment} ${meta.sample} ${meta.run_id} ${meta.type} ${meta.barcode} >> TEST

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fqcrs_to_hive.R: \$(echo "no-version" )
    END_VERSIONS
    """
}
