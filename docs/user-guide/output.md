# Expected Outputs

Upon successful completion, CARLISLE generates a comprehensive directory structure under `WORKDIR/results`. Each subdirectory contains outputs corresponding to specific stages of the workflow — from raw alignment statistics to annotated peak results and quality control summaries.

---

## Directory Overview

- **`alignment_stats/`** – Contains detailed alignment reports for each sample, including mapping efficiency, read depth, and spike-in alignment metrics.

- **`bam/`** – Stores sorted and indexed **BAM** files for all samples. This directory also includes per-sample and spike-in alignment statistics useful for downstream normalization and QC.

- **`bedgraph/`** – Includes **BEDGRAPH** coverage tracks summarizing read density across the genome. These files serve as intermediates for visualization and peak-calling validation.

- **`bigwig/`** – Contains **BigWig** files generated via `bamCoverage --scaleFactor`, where the scale factor reflects the spike-in read count or library size depending on `norm_method`. These are the same normalization values applied to the bedgraph tracks, ensuring genome browser tracks and deepTools heatmaps/profiles are directly comparable. Suitable for visualization in **[UCSC Genome Browser](https://genome.ucsc.edu/)** or **[IGV](https://igv.org/)**.

  > ℹ️ **Note:** There are no separately RPGC-normalized bigwigs. All bigwigs — including those consumed by deepTools heatmap and profile rules — come from this single `bamCoverage` step.

- **`fragments/`** – Stores fragment length distributions and deduplicated fragment data (particularly important for **CUT&RUN** and **CUT&Tag** experiments). Useful for assessing fragment size enrichment and MNase digestion efficiency.

- **`peaks/`** – The core results directory containing called peaks, differential comparisons, and annotations.

  - Subdirectories are organized by **quality thresholds** (e.g., `0.05`, `0.01`), representing the significance cutoffs applied during peak calling.
  - Each quality threshold directory includes:

    - **`contrasts/`** – Contains results of differential binding analyses defined in the contrast manifest, including:
      - Differential enrichment results from DESeq2 (AUC-based and fragment-based)
      - 3-column BED files for up-regulated peaks in each group (`up_group1.bed`, `up_group2.bed`)
      - **`homer_deg/`** – HOMER motif enrichment analysis for differentially enriched peaks (when `run_motif_enrichment_deg_peaks: true`)
      - **`ame_deg/`** – AME (Analysis of Motif Enrichment) results for differentially enriched peaks using [HOCOMOCO v14 CORE](https://hocomoco11.autosome.org/) motifs (when `run_motif_enrichment_deg_peaks: true`)
      - **`all.peaks.with_control.txt`** / **`all.peaks.without_control.txt`** – Aggregate peak count tables (filename encodes `run_without_controls` setting; older versions used `all.peaks.txt`)
      - **`Peak_counts.with_control.xlsx`** / **`Peak_counts.without_control.xlsx`** – Excel peak count summaries (older versions used `Peak counts.xlsx` with a space — update any downstream scripts accordingly)

    - **`<peak_caller>/`** – Subdirectories for each peak caller (e.g., `macs2`, `seacr`, `gopeaks`). Each includes raw peak calls and annotated results.

      - **`peak_output/`** – Raw peak calls organized by control mode:
        - **`individual/`** – Peaks called using individual replicate controls (present for all analyses)
        - **`pooled/`** – Peaks called using merged high-depth controls (present when `pool_controls: true`)
      - **`annotation/`** – Contains enriched feature and pathway analyses, organized by control mode:

        - **`go_enrichment/`** – Results from **[ChIP-Enrich](https://chipenrich.med.umich.edu/)** gene set enrichment, generated when `run_go_enrichment: true` is enabled.
        - **`homer/`** – Output from **[HOMER](http://homer.ucsd.edu/homer/)** motif discovery and annotation.
        - **`rose/`** – Output from **[ROSE](https://github.com/younglab/ROSE)** super-enhancer analysis, generated when `run_rose: true` is specified.

- **`qc/`** – Centralized quality control directory containing:
  - **`fastqc_raw/`** – Per-sample [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports on raw reads
  - **`fqscreen_raw/`** – [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) contamination reports
  - **`multiqc_report.html`** – Aggregated [MultiQC](https://multiqc.info/) report covering FastQC, alignment, and duplication statistics across all samples
  - **`spikein/`** – Spike-in normalization QC plots and tables (present when `norm_method: spikein`), including per-sample spike-in alignment rates and scaling factors

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
│   │   │   ├── homer_deg/                          # DEG motif enrichment (HOMER)
│   │   │   ├── ame_deg/                            # DEG motif enrichment (AME/HOCOMOCO)
│   │   │   ├── up_group1.bed                       # 3-column BED for group1 up-regulated peaks
│   │   │   ├── up_group2.bed                       # 3-column BED for group2 up-regulated peaks
│   │   │   ├── all.peaks.with_control.txt          # aggregate peak counts (with controls)
│   │   │   └── Peak_counts.with_control.xlsx       # Excel peak counts (with controls)
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
    ├── fqscreen_raw/
    ├── multiqc_report.html
    └── spikein/                   # Spike-in QC plots/tables (when norm_method: spikein)
```

---

> 🧭 **Tip:** The structure is intentionally hierarchical, enabling automated report generation and simplifying downstream integration with visualization tools and statistical frameworks.

---

## Interpreting Results

### QC: What to Check First

**Alignment stats / MultiQC** (`alignment_stats/`, `qc/multiqc_report.html`):

- Genome alignment rate should be >70% for typical CUT&RUN. Rates <50% suggest adapter contamination, wrong genome, or library prep failure.
- Spike-in alignment rate should be consistent across replicates within a condition. High CV in spike-in rates makes spike-in normalization unreliable; fall back to `norm_method: library`.
- Duplication rates for CUT&RUN/CUT&Tag are typically **lower** than ChIP-seq (~10–40%) due to MNase-based cleavage. Very low duplication (<5%) may indicate insufficient library complexity. Very high duplication (>60%) suggests over-amplification or low-complexity input.

**Fragment length distribution** (`fragments/`):

- For histone marks (H3K4me3, H3K27me3, etc.), expect a nucleosomal ladder: a dominant sub-nucleosomal peak (~150–200 bp) for mono-nucleosome protection and smaller peaks at ~300–400 bp and ~500–600 bp. CUT&RUN enriches for this pattern more sharply than ChIP-seq.
- For transcription factors, expect a dominant sub-nucleosomal peak at ~100–150 bp (accessible DNA footprint). A flat or bimodal distribution without a sub-nucleosomal peak indicates poor antibody specificity or failed experiment.

**Spike-in QC** (`qc/spikein/`):

- Examine the scaling factors table. Scaling factors should be inversely proportional to genome-aligned read depth — high-yield samples get smaller scale factors. If scaling factors vary >10-fold within a condition, review whether spike-in input was consistent across samples.

---

### Peak Calls: Choosing and Filtering

**Peak caller selection:**

| Caller | Best for | Output type | Notes |
|---|---|---|---|
| GoPeaks (narrow) | TFs, sharp histone marks (H3K4me3) | narrowPeak | Recommended for most CUT&RUN; explicitly optimized for low-background assays |
| GoPeaks (broad) | Broad histone marks (H3K27me3, H3K9me3) | broadPeak | — |
| MACS2 (narrow) | TFs, sharp marks | narrowPeak | Use when comparing to ChIP-seq datasets; widely adopted |
| MACS2 (broad) | Broad marks | broadPeak | DESeq2 differential analysis often fails on broadPeak due to excessive peak counts |
| SEACR (stringent/relaxed) | CUT&RUN-specific | bedGraph-derived | Uses AUC signal thresholding; no model fitting; less sensitive but very specific |

**Quality thresholds:** Peaks are called at each threshold in `quality_thresholds` independently. For discovery analyses, start with `q=0.05`. Tighten to `q=0.01` for high-confidence sets used in motif enrichment or validation. Do not use a single threshold for all downstream analyses — cross-threshold consistency is itself a signal of peak robustness.

**Individual vs. pooled controls:**

- Individual-mode peaks are more sensitive to replicate-specific variation but can include replicate-specific artifacts.
- Pooled-mode peaks benefit from increased control depth and lower noise floors, but may miss low-occupancy binding events present only in one replicate.
- Use **pooled-mode peaks as the primary call set** when replicates are well-correlated (Pearson r > 0.9 in deepTools); fall back to individual if replicates diverge.

---

### Differential Enrichment: DESeq2 Results

DESeq2 tests for differential AUC (area under the coverage curve) or fragment counts at consensus peak loci between two conditions. Results in `contrasts/` include:

- **`*_AUC_DESeq2_results.tsv`** — AUC-based: integrates signal across the peak; better for broad marks or when peak boundaries shift between conditions.
- **`*_fragments_DESeq2_results.tsv`** — Fragment count-based: more analogous to RNA-seq count models; better for narrow marks with stable peak boundaries.

Key columns to interpret:

| Column | Interpretation |
|---|---|
| `log2FoldChange` | Direction and magnitude; filtered at `contrasts_lfc_cutoff` (default 0.59 ≈ 1.5-fold) |
| `padj` | BH-adjusted p-value; filtered at `contrasts_fdr_cutoff` (default 0.05) |
| `baseMean` | Mean normalized signal across all samples; very low baseMean peaks are often noise |

> Peaks with `baseMean < 5` that pass FDR are suspect — flag them before reporting. DESeq2 can return significant p-values for very low-count features under low dispersion estimates.

**Up-regulated BED files** (`up_group1.bed`, `up_group2.bed`): 3-column BED of peaks with significantly higher enrichment in each group. Use these directly as input for downstream tools (GREAT, HOMER, AME, bedtools intersect with other genomic features).

---

### Motif Enrichment: HOMER and AME

**HOMER** (`homer_deg/`, `homer/`): Performs de novo motif discovery and known motif enrichment against the JASPAR/HOMER database. The top de novo motifs for DEG peaks point to TFs whose binding sites are preferentially accessible or occupied under that condition.

- The `knownResults.html` file ranks enrichment of HOMER's curated motif library. Focus on motifs with `p < 1e-5` and `% of target sequences > 10%`.
- De novo motifs that match known TF families validate antibody specificity and experimental biology; unexpected motifs suggest co-occupancy or indirect recruitment.

**AME** (`ame_deg/`): Tests enrichment of [HOCOMOCO v14 CORE](https://hocomoco11.autosome.org/) TF motifs (human/mouse) in up-regulated peak sequences vs. a shuffled background. HOCOMOCO v14 covers ~900 TF binding models with high positional weight matrix precision.

- `ame.tsv` columns: `motif_ID`, `adj_p-value`, `E-value`, `TP` (true positives), `FP` (false positives). Use `adj_p-value < 0.05` and `TP / (TP + FP) > 0.1` as a conservative filter.
- Cross-reference HOMER and AME hits: motifs appearing in both analyses under the same condition are high-confidence candidates for further validation (e.g., ChIP-seq, luciferase assay).

---

### GO / Pathway Enrichment: ChIP-Enrich

`go_enrichment/` contains [ChIP-Enrich](https://chipenrich.med.umich.edu/) results, which correct for gene locus length when testing peak-to-gene assignment. This avoids the systematic bias in standard Fisher-based enrichment where long genes accumulate more peaks by chance.

- **`chipenrich`** method: assigns peaks to nearest TSS within a configurable window and tests genesets. Best for narrow marks (TFs, H3K4me3) where peaks map to specific promoter/enhancer-gene links.
- **`polyenrich`** method: counts peaks per gene locus rather than presence/absence; better for broad marks.
- Key output file: `*_enrichment_results.tsv`. Columns `FDR` and `Effect` are most informative. Positive `Effect` = geneset enriched relative to background gene set size.

---

### Super-Enhancers: [ROSE](https://github.com/younglab/ROSE)

`rose/` identifies super-enhancers by stitching H3K27ac (or other active mark) peaks within `stitch_distance` (default 12,500 bp), then applying an inflection-point cutoff on the ranked signal curve. The inflection point separates typical enhancers from super-enhancers.

- **`*_SuperEnhancers.bed`**: the final super-enhancer calls. Overlap with known oncogenes or lineage TFs to nominate driver regulatory elements.
- **`*_Enhancers_withSuper.bed`**: all stitched enhancers ranked by signal; use this to inspect the inflection point and evaluate whether the cutoff is biologically reasonable.
- ROSE is only meaningful for **active chromatin marks** (H3K27ac, H3K4me1). Running it on H3K27me3 or TF datasets is not informative.

