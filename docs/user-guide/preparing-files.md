# Preparing Files

The CARLISLE pipeline is configured and controlled through a set of editable configuration and manifest files. After running `carlisle --runmode=init --workdir=/path/to/workdir` (see [Running the Pipeline](run.md)), default templates for these files are automatically generated under `WORKDIR/config/`.

---

## Configuration Files

CARLISLE’s configuration system is modular and designed for both flexibility and transparency. The main configuration files include:

- `config/config.yaml` – global pipeline settings and user parameters.
- `resources/cluster.yaml` – cluster resource specifications for **[Biowulf](https://hpc.nih.gov/)** or other SLURM-based systems.
- `resources/tools.yaml` – software versions, tool paths, and binary locations.

### Cluster Configuration (`cluster.yaml`)

The cluster configuration file defines computational resources such as memory, CPU cores, and runtime limits for each Snakemake rule. Parameters can be adjusted globally or per rule. Edits should be made with caution, as inappropriate resource settings may cause job failures or queuing delays.

### Tools Configuration (`tools.yaml`)

This file specifies which versions of each tool are used during execution. When running on Biowulf, tools are automatically loaded from environment modules, ensuring consistency across users. Once CARLISLE transitions to containers, these version pins will map to container image tags instead of module versions, guaranteeing strict reproducibility.

### Primary Configuration (`config.yaml`)

The main configuration file (`config.yaml`) contains parameters grouped into logical sections:

- **Folders and Paths:** defines input/output directories and manifest file locations.
- **User Parameters:** controls feature-level behavior (e.g., thresholds, normalization methods, peak calling options).
- **References:** specifies genome assemblies, index paths, spike-in references, and species annotations.

> ⚠️ **Important:** Always verify that reference genome paths and spike-in references correspond to accessible Biowulf or shared filesystem locations.

---

## User Parameters

### Spike-in Controls

CARLISLE supports spike-in normalization using reference genomes such as _E. coli_ or _Drosophila melanogaster_. The parameter `spikein_genome` defines the spike-in species, and `spikein_reference` provides the corresponding FASTA path.

Example for _E. coli_ spike-in:

```yaml
run_contrasts: true
norm_method: "spikein"
spikein_genome: "ecoli"
spikein_reference:
  ecoli:
    fa: "PIPELINE_HOME/resources/spikein/Ecoli_GCF_000005845.2_ASM584v2_genomic.fna"
```

Example for _Drosophila_ spike-in:

```yaml
run_contrasts: true
norm_method: "spikein"
spikein_genome: "drosophila"
spikein_reference:
  drosophila:
    fa: "/fdb/igenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa"
```

Example for _Saccharomyces cerevisiae_ spike-in:

```yaml
norm_method: "spikein"
spikein_genome: "saccharomyces"
spikein_reference:
  saccharomyces:
    fa: "$PIPELINE_HOME/resources/spikein/S_cer_S288C_R64.fna"
```

If spike-ins are unavailable or insufficient, normalization can alternatively be performed based on library size. Recommended workflow:

1. Run CARLISLE with `norm_method: spikein` for an initial QC assessment.
2. Evaluate spike-in alignment statistics (found in `alignment_stats/alignment_stats.tsv` in your results directory).
3. Change `norm_method` to `library` in your `config.yaml`.
4. Re-run CARLISLE — pooled control outputs are automatically regenerated when `norm_method` changes.

> ℹ️ **Normalization change behavior:** Pooled control fragment and bedgraph filenames encode the normalization method (e.g., `*.spikein.bedgraph`). Changing `norm_method` in an existing results directory causes Snakemake to detect stale targets and regenerate them automatically — **no manual deletion of intermediate files is required**.

> ℹ️ **Note:** `alignment_stats.tsv` is generated automatically by the pipeline and does not need to be specified in your configuration.

### Duplication Status

Control deduplication behavior using the `dupstatus` parameter:

```yaml
dupstatus: "dedup, no_dedup"
```

> ✅ **Recommendation:** Keep this setting unchanged, let CARLISLIE run with dedup and no_dedup options and then choose which peakSets to use later.

> 🧬 **Note:** Linear deduplication is essential for CUT&RUN and CUT&Tag datasets to avoid PCR bias and ensure accurate read quantification.

### Peak Callers

CARLISLE supports three major peak callers, configurable via the `peaktype` parameter:

1. **[MACS2](https://github.com/macs3-project/MACS)** – supports `narrowPeak` and `broadPeak` modes.
2. **[SEACR](https://seacr.fredhutch.org/)** – supports stringent and relaxed thresholds, for both normalized and non-normalized datasets.
3. **[GoPeaks](https://github.com/maxsonBraunLab/gopeaks)** – optimized for CUT&RUN and CUT&Tag data; recommended for most applications.

> ✅ **Recommendation:** Use GoPeaks for its superior signal detection in sparse chromatin accessibility datasets.

Example configuration:

```yaml
peaktype: "macs2_narrow, gopeaks_narrow"
```

### MACS2 Control Option

Enable control sample usage for MACS2 to improve specificity:

```yaml
macs2_control: "Y"
```

### Optional Analysis Steps

Control execution of computationally intensive annotation steps:

```yaml
run_rose: false                          # ROSE super-enhancer analysis (set to true to enable)
run_go_enrichment: false                 # ChIP-Enrich GO enrichment (set to true to enable)
run_motif_enrichment_called_peaks: false # HOMER motif discovery on all called peaks
run_motif_enrichment_deg_peaks: false    # HOMER + AME motif enrichment on DEG peaks only
```

> ⏱️ **Performance Note:** ROSE, GO enrichment, and motif enrichment are disabled by default due to their computational requirements. Enable them when you need super-enhancer identification, pathway enrichment, or motif discovery.

- **`run_motif_enrichment_called_peaks`**: When `true`, runs HOMER `findMotif` on the full set of called peaks for each sample/condition.
- **`run_motif_enrichment_deg_peaks`**: When `true`, runs both HOMER motif discovery **and** AME (Analysis of Motif Enrichment) against the [HOCOMOCO v14 CORE](https://hocomoco11.autosome.org/) motif database on up-regulated peak BED files (`up_group1.bed`, `up_group2.bed`) from each contrast. Both tools must be enabled together for full DEG motif enrichment output.

When `run_go_enrichment: true`, additional parameters control the enrichment methods and gene sets used:

```yaml
go_enrichment_methods: "chipenrich"  # options: chipenrich, polyenrich, hybridenrich
geneset_id: "GOBP,GOCC,GOMF,kegg_pathway,reactome"
```

Available `geneset_id` values include: `biocarta_pathway`, `ctd`, `cytoband`, `drug_bank`, `GOBP`, `GOCC`, `GOMF`, `hallmark`, `immunologic`, `kegg_pathway`, `mesh`, `metabolite`, `microrna`, `oncogenic`, `panther_pathway`, `pfam`, `protein_interaction_biogrid`, `reactome`, `transcription_factors`.

> ⚠️ **Performance:** `hybridenrich` is significantly slower than `chipenrich` or `polyenrich`. Only add it when its model is specifically required. GO enrichment is only supported for `hg19` and `hg38` samples.

### Pooled Controls

Control whether the pipeline pools control replicates for peak calling:

```yaml
pool_controls: true
```

When enabled (`true`), CARLISLE runs peak calling in **both modes**:

- **Individual mode** – Each treatment replicate is paired with its individual control replicate
- **Pooled mode** – Each treatment replicate is compared against merged high-depth controls from all control replicates

This dual-mode analysis enables comparison of replicate-specific vs merged control strategies. Results are organized in separate `individual/` and `pooled/` subdirectories within peak calling outputs.

> 💡 **Use Case:** Pooled controls provide increased depth and reduced noise but may miss replicate-specific artifacts. Running both modes allows downstream selection of the most appropriate strategy.

> ⚠️ **Note:** If controls have no replicates to pool (each control has only 1 replicate), pooling will have no effect. Consider setting `pool_controls: false` in such cases.

### Singularity Cache Directory

CARLISLE uses Singularity/Apptainer containers for R-based steps (DESeq2, GO enrichment, ROSE). **Most users do not need to configure this.** Loading the `ccbrpipeliner` module automatically sets `SIFCACHE` to the shared CCBR container cache, so pre-pulled images are used immediately:

```bash
module load ccbrpipeliner
# SIFCACHE is now set to the shared cache — no further action required
```

If you need to override the cache location (e.g., for a custom or updated image), you can do so in two ways:

```bash
# Option 1: Pass at runtime (takes precedence over everything)
carlisle --runmode=run --workdir=/path/to/workdir --singcache=/path/to/sif/cache

# Option 2: Set as an environment variable before running
export SIFCACHE=/path/to/sif/cache
carlisle --runmode=run --workdir=/path/to/workdir
```

If neither `--singcache` nor `$SIFCACHE` is set, CARLISLE resolves the cache directory in the following order:

1. `/data/${USER}/.singularity` — if `/data/${USER}/` exists on the filesystem (standard on Biowulf)
2. `${WORKDIR}/.singularity` — fallback when `/data/${USER}/` is not available

> ⚠️ **Warning:** Pointing to a cache directory that does not already contain the required `.sif` files will cause Singularity to pull all container images from Docker Hub. This can take **significant time** (depending on network conditions) and **consume several gigabytes of disk space**. Use the shared CCBR cache whenever possible.

### Control Sample Requirements

By default, **CARLISLE requires a control sample** (e.g., IgG, input DNA) paired with every treatment sample. Each non-control row in the sample manifest must have `controlName` and `controlReplicateNumber` filled in.

If you do not have control samples, **control-free mode is supported** — see [Control-Free Mode](#control-free-mode) below.

> 💡 **No controls?** Options include:
> - Using publicly available IgG controls from similar cell types ([GEO](https://www.ncbi.nlm.nih.gov/geo/), [ENCODE](https://www.encodeproject.org/))
> - Using input DNA as a proxy control
> - Enabling `run_without_controls: true` (see next section)

### Control-Free Mode

When no IgG or antibody control samples are available, set `run_without_controls: true` to run all peak callers without a control:

```yaml
run_without_controls: true
seacr_threshold: 0.01  # numeric FDR threshold used by SEACR in control-free mode
```

When enabled:

- **All treatment samples** are called as peaks against no background control.
- `macs2_control` is automatically forced to `"N"` (no-control MACS2 mode).
- `pool_controls` is automatically forced to `false`.
- **SEACR** uses the numeric `seacr_threshold` value instead of a control bedgraph. The value represents the fraction of the signal distribution used as the peak-calling threshold (e.g., `0.01` = top 1%).
- **GoPeaks** runs without the `-c` control BAM flag.
- The sample manifest **does not** require `controlName` or `controlReplicateNumber` columns to be filled in.

> ⚠️ **Caution:** Control-free peak calling will yield higher false-positive rates. Results should be interpreted with care and ideally validated by comparing to matched control experiments.

### deepTools Annotation Tracks

Control which annotation BED files are used to build deepTools coverage heatmaps and profiles:

```yaml
deeptools_bedtypes: "geneinfo,protein_coding,ca_ctcf,ca_h3k4me3,ca_tf,pls,pels"
```

Available options (comma-separated, no spaces):

cCRE bedtypes (`pls`, `pels`, `dels`, `ca_ctcf`, `ca_h3k4me3`, `ca_tf`) are sourced from the [ENCODE SCREEN database](https://screen.encodeproject.org/).

| Bedtype | Description |
|---|---|
| `geneinfo` | All genes: gene bodies, promoters, intergenic regions |
| `protein_coding` | Protein-coding genes only |
| `pls` | ENCODE SCREEN promoter-like signatures |
| `pels` | ENCODE SCREEN proximal enhancer-like signatures |
| `dels` | ENCODE SCREEN distal enhancer-like signatures |
| `ca_ctcf` | CTCF-bound chromatin accessibility regions (ENCODE SCREEN) |
| `ca_h3k4me3` | H3K4me3-marked chromatin accessibility / active promoters (ENCODE SCREEN) |
| `ca_tf` | Transcription factor-bound chromatin accessibility (ENCODE SCREEN) |

> ⚠️ **Memory Warning:** The `dels` BED file (ENCODE dELS) is very large. Including `dels` in `deeptools_bedtypes` requires `>=240g` memory for the `deeptools_mat` and `deeptools_heatmap` rules. Update the corresponding entries in `cluster.yaml` before enabling it.

---

### Quality Thresholds

Set peak-calling quality thresholds using the `quality_thresholds` parameter:

```yaml
quality_thresholds: "0.1, 0.05, 0.01"
```

Refer to tool-specific defaults:

- MACS2 q-value: [0.05](https://manpages.ubuntu.com/manpages/xenial/man1/macs2_callpeak.1.html)
- GoPeaks p-value: [0.05](https://github.com/maxsonBraunLab/gopeaks#usage)
- SEACR FDR threshold: [1.0](https://github.com/FredHutch/SEACR#usage)

### MACS2 Broad Peak Threshold

For MACS2 broad peak calling (`macs2_broad`), an additional p-value threshold is applied independently of the global `quality_thresholds`:

```yaml
macs2_broad_peak_threshold: "0.01"
```

This maps to the `--broad-cutoff` MACS2 argument and controls the significance cutoff for the broad-region merging step. The default of `0.01` is generally appropriate; however, reduce it (e.g., `0.001`) if broad peaks are excessively numerous or fragmented.

> ℹ️ **Note:** DESeq2 differential analysis frequently fails for broadPeak outputs due to excessive peak counts. For differential analysis, `macs2_narrow` or `gopeaks_narrow` are recommended.

### Differential Analysis Thresholds

DESeq2 significance cutoffs for contrast-based differential enrichment:

```yaml
contrasts_fdr_cutoff: 0.05   # Benjamini-Hochberg adjusted p-value (FDR) threshold
contrasts_lfc_cutoff: 0.59   # log2 fold-change threshold (~1.5-fold change)
```

Both thresholds are applied simultaneously: a peak is considered differentially enriched only if it passes both FDR and log2FC filters. Adjust `contrasts_lfc_cutoff` to `1.0` (2-fold) for more conservative enrichment calls, or lower both thresholds if the experiment has high biological variability.

---

## Reference Files

CARLISLE includes comprehensive reference annotations for supported genomes:

### Built-in Annotations

For each genome (hg38, hg19, hs1/T2T, mm10), the pipeline provides:

- **Gene annotations**: TSS, gene bodies, promoters, intergenic regions (protein-coding and all genes)
- **Blacklisted regions**: ENCODE DAC blacklists for artifact exclusion
- **cCREs (candidate cis-Regulatory Elements)**: From the [ENCODE SCREEN database](https://screen.encodeproject.org/)
  - **PLS** – Promoter-like signatures
  - **pELS** – Proximal enhancer-like signatures
  - **dELS** – Distal enhancer-like signatures
  - **CA-CTCF** – CTCF-bound chromatin accessibility regions
  - **CA-H3K4me3** – H3K4me3-marked chromatin accessibility (active promoters)
  - **CA-TF** – Transcription factor-bound chromatin accessibility
  - **CA** – General chromatin accessibility
  - **TF** – Transcription factor binding sites

These annotations are automatically used by HOMER, GO enrichment, and other annotation tools.

### Custom Genomes

Additional reference genomes can be integrated by defining:

```yaml
species_name:
  fa: "/path/to/species.fa"
  blacklist: "/path/to/blacklistbed/species.bed.gz"
  regions: "chr1 chr2 chr3"
  macs2_g: "hs" # genome shorthand for MACS2
  tss_bed: "/path/to/tss.bed.gz"
  # Add cCRE annotations if available
  ca_pls_bed: "/path/to/cCREs.PLS.bed.gz"
  ca_pels_bed: "/path/to/cCREs.pELS.bed.gz"
  ca_dels_bed: "/path/to/cCREs.dELS.bed.gz"
```

> 🧭 **Best Practice:** Store reference paths under a centralized `/fdb` or `/data` location on Biowulf to ensure accessibility and consistency across users.

---

## Preparing Manifests

CARLISLE uses two manifests:

- `samplemanifest` – required for all analyses.
- `contrasts` – optional, required only for differential analysis with DESeq2.

### Sample Manifest (Required)

Defines sample-level metadata, including sample names, controls, and FASTQ paths.

| sampleName                      | replicateNumber | isControl | controlName                     | controlReplicateNumber | path_to_R1                                              | path_to_R2                                              |
| ------------------------------- | --------------- | --------- | ------------------------------- | ---------------------- | ------------------------------------------------------- | ------------------------------------------------------- |
| 53_H3K4me3                      | 1               | N         | HN6_IgG_rabbit_negative_control | 1                      | <path_to>/53_H3K4me3_1.R1.fastq.gz                      | <path_to>/53_H3K4me3_1.R2.fastq.gz                      |
| 54_H3K4me3                      | 2               | N         | HN6_IgG_rabbit_negative_control | 1                      | <path_to>/54_H3K4me3_1.R1.fastq.gz                      | <path_to>/54_H3K4me3_1.R2.fastq.gz                      |
| HN6_IgG_rabbit_negative_control | 1               | Y         |                                 |                        | <path_to>/HN6_IgG_rabbit_negative_control_1.R1.fastq.gz | <path_to>/HN6_IgG_rabbit_negative_control_2.R2.fastq.gz |

> ℹ️ **Note:** `controlName` and `controlReplicateNumber` are **required** for non-control samples in normal mode. In control-free mode (`run_without_controls: true`), leave these columns blank for all samples and omit control rows entirely.

> ⚠️ **Sample name uniqueness:** Sample names must not be substrings of each other (e.g., having both `H3K4me3` and `H3K4me3_rep1` as `sampleName` values will cause incorrect sample matching). Use fully distinct names for all samples.

### Contrast Manifest (Optional)

Specifies conditions for differential analysis:

| condition1              | condition2           |
| ----------------------- | -------------------- |
| MOC1_siSmyd3_2m_25_HCHO | MOC1_siNC_2m_25_HCHO |

> 📊 **Requirement:** Each condition must have at least two biological replicates to perform DESeq2-based differential analysis.
