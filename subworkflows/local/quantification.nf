include { HTSEQ_COUNT                                 } from '../../modules/local/htseq/count/main'
include { KALLISTO_QUANT                              } from '../../modules/nf-core/kallisto/quant/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_GENE } from '../../modules/nf-core/subread/featurecounts/main'
include { COUNT_FEATURES                              } from '../../modules/local/count_features/main'
include { CUSTOM_TX2GENE                              } from '../../modules/nf-core/custom/tx2gene/main'
include { TXIMETA_TXIMPORT                            } from '../../modules/nf-core/tximeta/tximport/main'


workflow QUANTIFICATION {
    take:
    bam
    bai
    gtf
    reads
    kallisto_idx

    main:

    ch_versions = Channel.empty()


    HTSEQ_COUNT(
        bam.join(bai,by:[0]),
        gtf.map{it[1]}
    )
    ch_versions   = ch_versions.mix(HTSEQ_COUNT.out.versions)


    FEATURECOUNTS_GENE(
        bam,
        gtf.map{it[1]}
    )
    ch_versions = ch_versions.mix(FEATURECOUNTS_GENE.out.versions)

    KALLISTO_QUANT(
        reads,
        kallisto_idx,
        [],
        []
    )
    ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions)

    ch_kallisto_grouped = KALLISTO_QUANT.out.abundance
        .map { meta, abundance ->
            def meta_sample = [id: meta.sample]
            [meta_sample, abundance]
        }
        .groupTuple()

    CUSTOM_TX2GENE(
        gtf,
        ch_kallisto_grouped,
        "kallisto",
        "gene_id",
        "gene_name"
    )
    ch_versions = ch_versions.mix(CUSTOM_TX2GENE.out.versions)

    TXIMETA_TXIMPORT(
        ch_kallisto_grouped,
        CUSTOM_TX2GENE.out.tx2gene,
        "kallisto"
    )
    ch_versions = ch_versions.mix(TXIMETA_TXIMPORT.out.versions)

    COUNT_FEATURES(
        KALLISTO_QUANT.out.abundance,
        gtf.map{it[1]}
    )


    emit:
    htseq_counts               = HTSEQ_COUNT.out.counts
    htseq_summary              = HTSEQ_COUNT.out.summary
    featurecounts_gene_counts  = FEATURECOUNTS_GENE.out.counts
    featurecounts_gene_summary = FEATURECOUNTS_GENE.out.summary
    kallisto_log               = KALLISTO_QUANT.out.log
    kallisto_count_feature     = COUNT_FEATURES.out.kallisto_count_feature
    tximport_gene_tpm          = TXIMETA_TXIMPORT.out.tpm_gene
    tximport_gene_counts       = TXIMETA_TXIMPORT.out.counts_gene
    tximport_gene_counts_ls    = TXIMETA_TXIMPORT.out.counts_gene_length_scaled
    tximport_gene_counts_s     = TXIMETA_TXIMPORT.out.counts_gene_scaled
    tximport_gene_lengths      = TXIMETA_TXIMPORT.out.lengths_gene
    tximport_transcript_tpm    = TXIMETA_TXIMPORT.out.tpm_transcript
    tximport_transcript_counts = TXIMETA_TXIMPORT.out.counts_transcript
    tximport_transcript_lengths= TXIMETA_TXIMPORT.out.lengths_transcript
    ch_versions                = ch_versions
}
