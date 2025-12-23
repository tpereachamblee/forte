process FUSVIZ {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' ?
        'docker://blancojmskcc/target_fusviz:7.4.0':
        'blancojmskcc/target_fusviz:7.4.0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(tsv)
    path(cytobands)
    path(annotation)
    path(chromosomes)
    path(protein_domains)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    FusViz \\
        --fusions=${tsv} \\
        --alignments=${bam} \\
        --cytobands=${cytobands} \\
        --annotation=${annotation} \\
        --chromosomes=${chromosomes} \\
        --output=${prefix}_FusViz.pdf \\
        --proteinDomains=${protein_domains} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: \$(echo \$(FusViz --version 2>&1) | sed 's/^.*FusViz //;  s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_FusViz.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusviz: \$(echo \$(FusViz --version 2>&1) | sed 's/^.*FusViz //;  s/ .*\$//')
    END_VERSIONS
    """
}
