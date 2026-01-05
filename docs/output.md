# mskcc/forte: Output

## Introduction

This document describes the output produced by the FORTE pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Read Preprocessing](#read-preprocessing)
- [Alignment](#alignment)
- [Quantification](#quantification)
- [Fusion Calling](#fusion-calling)
- [Fusion Merging and Annotation](#fusion-merging-and-annotation)
- [Fusion Visualization](#fusion-visualization)
- [Splicing](#splicing)
- [QC](#qc)
- [Fillouts](#fillouts)
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Read Preprocessing

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/fastp/`
  - `*.fastp.html`
  - `*.fastp.json`
  - `*.fastp.log`
  - `*.fastp.fastq.gz`
- `analysis/<sample>/umitools/extract/`
  - `logs/.umi_extract.log`

</details>

[FastP](https://github.com/OpenGene/fastp) gives general quality metrics about your sequenced reads and also trims the reads according to base quality and presence of adapter sequences.

[UMI-tools extract](https://umi-tools.readthedocs.io/en/latest/reference/extract.html) removes UMI sequences from reads and adds it to the read header. As a result, aligners do not attempt to align the UMI sequence and the aligned reads will be ready for deduplication.

### Alignment

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/STAR/`
  - `*.Aligned.sortedByCoord.out.bam`
  - `*.Aligned.sortedByCoord.out.bam.bai`
  - `log/`
    - `*.Log.out`
    - `*.Log.final.out`
    - `*.Log.progress.out`
    - `*.SJ.out.tab`
- `analysis/<sample>/umitools/dedup/`
  - `*.dedup.bam`
  - `*.dedup.bam.bai`
  - `logs/`
    - `*.dedup_edit_distance.tsv`
    - `*.dedup_per_umi_per_position.tsv`
    - `*.dedup_per_umi.tsv`

</details>

[STAR](https://github.com/alexdobin/STAR) is an ultrafast universal RNA-seq aligner.

[UMI-tools dedup](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html) deduplicates reads based on the mapping co-ordinate and the UMI attached to the read.

### Quantification

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/featurecounts/`
  - `*.gene.featureCounts.txt`
- `analysis/<sample>/kallisto/`
  - `abundance.h5`
  - `abundance.tsv`
  - `run_info.json`
  - `*.log.txt`

</details>

[featureCounts](https://subread.sourceforge.net/featureCounts.html) takes a file with aligned sequencing reads, plus a list of genomic features and counts how many reads map to each feature.

[Kallisto](http://pachterlab.github.io/kallisto/) quantifies abundances of transcripts from RNA-Seq data using high-throughput sequencing reads.

### Fusion Calling

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/arriba/`
  - `*.fusions.discarded.tsv`
  - `*.fusions.tsv`
  - `*_arriba.cff`
- `analysis/<sample>/fusioncatcher/`
  - `*.fusioncatcher.fusion-genes.hg19.txt`
  - `*.fusioncatcher.fusion-genes.txt`
  - `*.fusioncatcher.log`
  - `*.fusioncatcher.summary.txt`
  - `*_fusioncatcher.cff`
  - `*.supporting-reads_gene-fusions*.zip`
- `analysis/<sample>/starfusion/`
  - `*.starfusion.abridged.coding_effect.tsv`
  - `*.starfusion.abridged.tsv`
  - `*.starfusion.fusion_predictions.tsv`
  - `*_starfusion.cff`
  - `STAR/`
    - `*.Chimeric.out.junction`
    - `log/`
      - `*.Log.final.out`
      - `*.Log.out`
      - `*.Log.progress.out`
      - `*.SJ.out.tab`

</details>

[Arriba](https://arriba.readthedocs.io/en/latest/) uses the STAR aligner to detect of gene fusions from RNA-Seq data.

[FusionCatcher](https://github.com/ndaniel/fusioncatcher) searches for novel/known somatic fusion genes, translocations, and chimeras in RNA-seq data.

[STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) uses the STAR aligner to identify candidate fusion transcripts supported by Illumina reads.

[CommonFusionFormat (CFF)](https://github.com/ccmbioinfo/MetaFusion/wiki/metafusion-file-formats) file is originally created by the developers of the [MetaFusion](https://github.com/mskcc/MetaFusion) tool. Forte generates `*.cff` files that can be used with MetaFusion.

### Fusion Merging and Annotation

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/metafusion`
  - `*.final.cff`
  - `*.unfiltered.cff`
  - `*_filtered_fusions.tsv`
  - `*_filtered_fusions_cvr.tsv`
  - `*_cis_sage_fusions.tsv`
  - `*_iannotatesv_input.tsv`
  - `*_iannotatesv_canoncicalTranscripts.tsv`
  - `intermediates/`
    - `cis-sage.cluster`
    - `*.cff.cleaned_chr.renamed.reann.WITH_SEQ.exons`
    - `*_metafusion_cluster.unfiltered.cff`
    - `final.n1.cluster`
    - `problematic_chromosomes.cff`
- `analysis/<sample>/agfusion_clinical`
  - `<sample>/`
  - `*.expanded_agfusion_transcripts.tsv`
- `analysis/<sample>/agfusion`
  - `<sample>/`

</details>

FORTE uses a custom fork of [Metafusion](https://github.com/mskcc/MetaFusion) to filter, cluster and annotate the fusion calls. Several `intermediate` files are included in the output, [see wiki for detailed information](https://github.com/mskcc/forte/wiki/Metafusion-Output).

`Fusion_effect` information is added using a custom fork of [AGFusion](https://github.com/anoronh4/AGFusion).

`*.expanded_agfusion_transcripts.tsv` file is available for manual review of all possible frame status's across all transcript combinations for a select list of clinical genes within `assets/clinical_genes.txt`.

`*_filtered_fusions.tsv`, `*_filtered_fusions_cvr.tsv`,`*_cis_sage_fusions.tsv`, `*_iannotatesv_canoncicalTranscripts.tsv` and `*_iannotatesv_input.tsv` are all generated by `bin/fusion_filtering.R`. This script takes in the final formated CFF with `Fusion_effect` column and generated a actionability column, as well as selects a single fusion per cluster. `*_filtered_fusions.tsv` and `*_cis_sage_fusions.tsv` files have identical columns, with the latter only containing cis_sage clusters which are not in our `assets/cis_sage_allowlist.txt`. The columns are as follows:

| Column          | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| --------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| sample          | Sample ID                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| cluster         | assigned cluster according to metafusion, this will map back to CFF output                                                                                                                                                                                                                                                                                                                                                                                                   |
| tool            | Which of the fusion callers called this fusion within the cluster, A = Arriba, F = FusionCatcher, S = StarFusion                                                                                                                                                                                                                                                                                                                                                             |
| fusion          | Fusion name in format of 5'Gene::3'Gene                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| total_support   | For the selected breakpoint, what was the max read support found across the tools                                                                                                                                                                                                                                                                                                                                                                                            |
| breakpoint      | Selected breakpoint                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| action          | REPORT - report this fusion clinically, NOVEL - this fusion is reportable but does not include a gene from `assets/clinical_genes.txt`, READ_THROUGH - this fusion is a cis_sage fusion, drop - this fusion did not pass filtering                                                                                                                                                                                                                                           |
| reason          | LOW_SUPPORT - all breakpoint in cluster have less any 5 reads supporting, ONE_CALLER - a single caller identified this fusion, NO_SIG_GENE - this fusion does not include a gene from `assets/clinical_genes.txt`, FP - this fusion has a known false positive flag set in `fc_flags` or `sf_flags`, CIS_SAGE - this fusion was labelled as cis_sage by metafusion                                                                                                           |
| fp_tools        | fc - fusioncatcher false positive, sf - starfusion false positive                                                                                                                                                                                                                                                                                                                                                                                                            |
| fp_flags        | output of false positive flags                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| frame_status_br | Tools which label the selected breakpoint as in frame. A = Arriba, F = FusionCatcher, S = StarFusion, G = AGFusion                                                                                                                                                                                                                                                                                                                                                           |
| frame_status_cl | Tools which label any breakpoint in the cluster as in frame. A = Arriba, F = FusionCatcher, S = StarFusion, G = AGFusion                                                                                                                                                                                                                                                                                                                                                     |
| tx5             | 5' Transcript ID                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| tx3             | 3' Transcript ID                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Fusion_Effect   | AGFusion frame status result for the selected breakpoint. in-frame: in frame fusion. out-of-frame: out of frame fusion. in-frame (with mutation): fusion is in frame but the fusion results in a mutation between the two genes. Typically happens when the 5' gene last exon in the fusion and the 3' gene first exon in the fusion straddles a codon leading to a functional condon when fused. exon-exon: a non-coding transcript is utilized for 5'gene, 3' gene or both |
| somatic_flags   | Known true positive somatic flags from fusioncatcher. Descriptions available [here](https://illumina.github.io/NirvanaDocumentation/data-sources/fusioncatcher/#somatic)                                                                                                                                                                                                                                                                                                     |

`*_filtered_fusions_cvr.tsv` is filtered to fusions with the REPORT action from `*_filtered_fusions.tsv`. This file follows designations from the annotation script which will be described in detail later. This file is a placeholder for upload to CVR while we work on the annotation script.

`*_iannotatesv_canoncicalTranscripts.tsv` and `*_iannotatesv_input.tsv` are inputs into [iAnnotateSV:msk-target branch](https://github.com/rhshah/iAnnotateSV/tree/msk-target). See iAnntateSV-msktarget readme for more information.

[FusionAnnotator.py from the oncokb-annotator](https://github.com/oncokb/oncokb-annotator/blob/master/FusionAnnotator.py) is also run and added to the final cff file.

### Fusion Visualization

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/fusviz/`
  - `<sample>_FusViz.pdf`

</details>

[FusViz](https://hub.docker.com/repository/docker/blancojmskcc/target_fusviz) is a Pythonic FUSion VIsualiZation tool, built and tailored for MSK-TARGET panel data. It generates PDF visualizations of detected gene fusions, showing the fusion breakpoints, gene structures, and supporting read alignments. The visualization includes cytoband information, chromosome context, and protein domain annotations to aid in the interpretation of fusion events.

### Splicing

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/portcullis/`
  - `*.portcullis.log`
  - `portcullis_filtered.pass.junctions.bed`
  - `portcullis_filtered.pass.junctions.tab`
  - `portcullis_filtered.pass.junctions.exon.gff3`
  - `portcullis_filtered.pass.junctions.intron.gff3`
  - `*_oncogenic_isoforms.txt`
  - `*_oncogenic_isoforms_dropped.txt`

</details>

[Portcullis](https://portcullis.readthedocs.io/en/latest/) analyzes and quantifies splice junctions from a BAM file.

FORTE uses a custom script to calculate the percentage of oncogenic isoforms found from Potcullis for EGFR, ARv7 and MET exon14 deletion. These oncogenic isoforms include deletion of EGFR exons 2 through 7, deletion of EGFR exons 14 and 15, AR variant 7 novel transcript which involves splicing of exon 3 to a cryptic exon 3 and MET exon 14 skipping.

### QC

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/picard/`
  - `*.rna_metrics`
  - `*.CollectHsMetrics.coverage_metrics`
- `analysis/<sample>/rseqc/`
  - `*.bam_stat.txt`
  - `*.DupRate_plot.pdf`
  - `*.DupRate_plot.r`
  - `*.infer_experiment.txt`
  - `*.inner_distance_freq.txt`
  - `*.inner_distance_mean.txt`
  - `*.inner_distance_plot.pdf`
  - `*.inner_distance_plot.r`
  - `*.inner_distance.txt`
  - `*.junction_annotation.log`
  - `*.junction.bed`
  - `*.junction.Interact.bed`
  - `*.junction_plot.r`
  - `*.junctionSaturation_plot.pdf`
  - `*.junctionSaturation_plot.r`
  - `*.junction.xls`
  - `*.pos.DupRate.xls`
  - `*.read_distribution.txt`
  - `*.seq.DupRate.xls`
  - `*.splice_events.pdf`
  - `*.splice_junction.pdf`
- `analysis/<sample>/multiqc/`
  - `dedupbam_multiqc_report_data/`
    - `*.json`
    - `*.log`
    - `*.txt`
  - `dedupbam_multiqc_report.html`
  - `dedupbam_multiqc_report_plots/`
    - `pdf/*.pdf`
    - `png/*.png`
    - `svg/*.svg`
  - `dupbam_multiqc_report_data/`
    - `*.json`
    - `*.log`
    - `*.txt`
  - `dupbam_multiqc_report.html`
  - `dupbam_multiqc_report_plots/`
    - `pdf/*.pdf`
    - `png/*.png`
    - `svg/*.svg`

</details>

[Picard's CollectHsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard-) collects hybrid-selection (HS) metrics for a SAM or BAM file. This is only produced if baitset is indicated in the samplesheet.

[Picard's CollectRnaSeqMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037057492-CollectRnaSeqMetrics-Picard-) produces RNA alignment metrics for a SAM or BAM file.

[RSeQC](https://rseqc.sourceforge.net/) provides a number of useful modules that can comprehensively evaluate high throughput RNAseq data.

[MultiQC](https://multiqc.info/) is a visualization tool that searches a given directory for analysis/qc logs and compiles a HTML report. Most of the pipeline QC results are visualized in the report and further statistics are available in the report data directory. FORTE produces a second MultiQC report for each sample that has UMI. FORTE also produces 1-2 reports under the `multiqc/` folder where all samples are aggregated together, one for non-deduplicated results and the other for deduplicated results.

### Fillouts

<details markdown="1">
<summary>Output files</summary>

- `analysis/<sample>/fillouts`
  - `*.fillout.maf`

</details>

[GetBaseCountsMultiSample (GBCMS)](https://github.com/zengzheng123/GetBaseCountsMultiSample) calculates the base counts in a given BAM file for all the sites in a given MAF file

FORTE uses a custom script to output a MAF file that combines all original columns and new columns from GBCMS.

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
