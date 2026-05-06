# Expected Outputs

Upon successful completion, CARLISLE generates a comprehensive directory structure under `WORKDIR/results`. Each subdirectory contains outputs corresponding to specific stages of the workflow — from raw alignment statistics to annotated peak results and quality control summaries.

---

## Directory Overview

- **`alignment_stats/`** – Contains detailed alignment reports for each sample, including mapping efficiency, read depth, and spike-in alignment metrics.

- **`bam/`** – Stores sorted and indexed **BAM** files for all samples. This directory also includes per-sample and spike-in alignment statistics useful for downstream normalization and QC.

- **`bedgraph/`** – Includes **BEDGRAPH** coverage tracks summarizing read density across the genome. These files serve as intermediates for visualization and peak-calling validation.

- **`bigwig/`** – Contains **BigWig** files generated from normalized coverage data, suitable for visualization in genome browsers such as **[UCSC Genome Browser](https://genome.ucsc.edu/)** or **[IGV](https://igv.org/)**.

- **`fragments/`** – Stores fragment length distributions and deduplicated fragment data (particularly important for **CUT&RUN** and **CUT&Tag** experiments). Useful for assessing fragment size enrichment and MNase digestion efficiency.

- **`peaks/`** – The core results directory containing called peaks, differential comparisons, and annotations.

  - Subdirectories are organized by **quality thresholds** (e.g., `0.05`, `0.01`), representing the significance cutoffs applied during peak calling.
  - Each quality threshold directory includes:
    - **`contrasts/`** – Contains results of differential binding analyses defined in the contrast manifest, including:
      - Differential enrichment results from DESeq2 (AUC-based and fragment-based)
      - 3-column BED files for up-regulated peaks in each group (`up_group1.bed`, `up_group2.bed`)
      - **`homer_deg/`** – HOMER motif enrichment analysis for differentially enriched peaks
      - **`ame_deg/`** – AME (Analysis of Motif Enrichment) results for differentially enriched peaks using HOCOMOCO motifs
    - **`<peak_caller>/`** – Subdirectories for each peak caller (e.g., `macs2`, `seacr`, `gopeaks`). Each includes raw peak calls and annotated results.
      - **`peak_output/`** – Raw peak calls organized by control mode:
        - **`individual/`** – Peaks called using individual replicate controls (present for all analyses)
        - **`pooled/`** – Peaks called using merged high-depth controls (present when `pool_controls: true`)
      - **`annotation/`** – Contains enriched feature and pathway analyses, organized by control mode:
        - **`go_enrichment/`** – Results from **[ChIP-Enrich](https://chipenrich.med.umich.edu/)** gene set enrichment, generated when `run_go_enrichment: true` is enabled.
        - **`homer/`** – Output from **[HOMER](http://homer.ucsd.edu/homer/)** motif discovery and annotation.
        - **`rose/`** – Output from **[ROSE](https://bitbucket.org/young_computation/rose/src/master/)** super-enhancer analysis, generated when `run_rose: true` is specified.

- **`qc/`** – Centralized quality control directory containing comprehensive **[MultiQC](https://multiqc.info/)** summaries, **FastQC** metrics, and spike-in normalization reports (when applicable).

---

## Example Directory Layout

Below is an example of the CARLISLE output structure for a typical CUT&RUN experiment:

```
results/
├── alignment_stats/
├── bam/
│   └── pooled_controls/           # Merged control BAMs (when pool_controls: true)
├── bedgraph/
│   └── pooled_controls/           # Pooled control bedgraphs (when pool_controls: true)
├── bigwig/
├── fragments/
│   └── pooled_controls/           # Pooled control fragments (when pool_controls: true)
├── peaks/
│   ├── 0.05/
│   │   ├── contrasts/
│   │   │   ├── homer_deg/         # DEG motif enrichment (HOMER)
│   │   │   ├── ame_deg/           # DEG motif enrichment (AME/HOCOMOCO)
│   │   │   ├── up_group1.bed      # 3-column BED for group1 up-regulated peaks
│   │   │   └── up_group2.bed      # 3-column BED for group2 up-regulated peaks
│   │   ├── gopeaks/
│   │   │   ├── peak_output/
│   │   │   │   ├── individual/    # Peaks with individual controls
│   │   │   │   └── pooled/        # Peaks with pooled controls
│   │   │   └── annotation/
│   │   │       ├── individual/
│   │   │       │   ├── go_enrichment/
│   │   │       │   ├── homer/
│   │   │       │   └── rose/
│   │   │       └── pooled/
│   │   │           ├── go_enrichment/
│   │   │           ├── homer/
│   │   │           └── rose/
│   │   ├── macs2/
│   │   │   ├── peak_output/
│   │   │   │   ├── individual/
│   │   │   │   └── pooled/
│   │   │   └── annotation/
│   │   │       ├── individual/
│   │   │       │   ├── go_enrichment/
│   │   │       │   ├── homer/
│   │   │       │   └── rose/
│   │   │       └── pooled/
│   │   │           ├── go_enrichment/
│   │   │           ├── homer/
│   │   │           └── rose/
│   │   └── seacr/
│   │       ├── peak_output/
│   │       │   ├── individual/
│   │       │   └── pooled/
│   │       └── annotation/
│   │           ├── individual/
│   │           │   ├── go_enrichment/
│   │           │   ├── homer/
│   │           │   └── rose/
│   │           └── pooled/
│   │               ├── go_enrichment/
│   │               ├── homer/
│   │               └── rose/
│   └── 0.01/
│       └── ...
└── qc/
    ├── fastqc_raw/
    └── fqscreen_raw/
```

---

> 🧭 **Tip:** The structure is intentionally hierarchical, enabling automated report generation and simplifying downstream integration with visualization tools and statistical frameworks.
