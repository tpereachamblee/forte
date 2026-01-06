process FUSION_FILTER {
    tag "$meta.id"
    label "process_single"

/// must be using singularity 3.7+
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/rocker-org/tidyverse:4.4.2' :
        'ghcr.io/rocker-org/tidyverse:4.4.2' }"

    input:
    tuple val(meta), path(cff), path(starfusion), path(fusioncatcher), path(arriba)
    path clinical_genes
    path fusioncatcher_ref
    tuple val(meta2), path(gtf)
    path cis_sage_allow

    output:
    tuple val(meta), path("*_filtered_fusions.tsv")                  , emit: filtered_fusions
    tuple val(meta), path("*_cis_sage_fusions.tsv")                  , emit: cis_sage_fusions
    tuple val(meta), path("*_cvr.tsv")                               , emit: cvr_fusions
    tuple val(meta), path("*_iannotatesv_input.tsv")                 , emit: iannotatesv_input
    tuple val(meta), path("*_iannotatesv_canoncicalTranscripts.tsv") , emit: iannotatesv_candidates
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def sample = "${meta.sample}"
    """
    fusion_filtering.R \\
        --cff ${cff} \\
        --starfusion ${starfusion} \\
        --fusioncatcher ${fusioncatcher} \\
        --arriba ${arriba} \\
        --clinical_genes ${clinical_genes} \\
        --out_prefix ${sample} \\
        --fc_reference_dir  ${fusioncatcher_ref} \\
        --gtf ${gtf} \\
        --cis_sage_allow ${cis_sage_allow} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        fusion_filtering.R: 0.1.0
    END_VERSIONS
    """
}
