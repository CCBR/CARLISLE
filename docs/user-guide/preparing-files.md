# Preparing Files

The CARLISLE pipeline is configured and controlled through a set of editable configuration and manifest files. Upon initialization, default templates for these files are automatically generated under the `/WORKDIR/config` and `/WORKDIR/manifest` directories.

> ⚙️ **Technical Note:** CARLISLE follows a Snakemake-driven workflow architecture where all configuration parameters are read dynamically at runtime. Users are encouraged to version-control configuration files (e.g., via Git) to ensure reproducibility across runs.

> 🚀 **Future Development:** While dependencies are currently module-loaded on the **[Biowulf HPC environment](https://hpc.nih.gov/)**, future releases will adopt containerization using **[Singularity/Apptainer](https://apptainer.org/)** and **[Docker](https://www.docker.com/)**. This shift will provide complete environment encapsulation, allowing consistent execution across HPC and cloud environments.

---

## Configuration Files

CARLISLE’s configuration system is modular and designed for both flexibility and transparency. The main configuration files include:

* `config/config.yaml` – global pipeline settings and user parameters.
* `resources/cluster.yaml` – cluster resource specifications for **[Biowulf](https://hpc.nih.gov/)** or other SLURM-based systems.
* `resources/tools.yaml` – software versions, tool paths, and binary locations.

### Cluster Configuration (`cluster.yaml`)

The cluster configuration file defines computational resources such as memory, CPU cores, and runtime limits for each Snakemake rule. Parameters can be adjusted globally or per rule. Edits should be made with caution, as inappropriate resource settings may cause job failures or queuing delays.


### Tools Configuration (`tools.yaml`)

This file specifies which versions of each tool are used during execution. When running on Biowulf, tools are automatically loaded from environment modules, ensuring consistency across users. Once CARLISLE transitions to containers, these version pins will map to container image tags instead of module versions, guaranteeing strict reproducibility.

### Primary Configuration (`config.yaml`)

The main configuration file (`config.yaml`) contains parameters grouped into logical sections:

* **Folders and Paths:** defines input/output directories and manifest file locations.
* **User Parameters:** controls feature-level behavior (e.g., thresholds, normalization methods, peak calling options).
* **References:** specifies genome assemblies, index paths, spike-in references, and species annotations.

> ⚠️ **Important:** Always verify that reference genome paths and spike-in references correspond to accessible Biowulf or shared filesystem locations.

---

## User Parameters

### Spike-in Controls

CARLISLE supports spike-in normalization using reference genomes such as *E. coli* or *Drosophila melanogaster*. The parameter `spikein_genome` defines the spike-in species, and `spikein_reference` provides the corresponding FASTA path.

Example for *E. coli* spike-in:

```yaml
run_contrasts: true
norm_method: "spikein"
spikein_genome: "ecoli"
spikein_reference:
  ecoli:
    fa: "PIPELINE_HOME/resources/spikein/Ecoli_GCF_000005845.2_ASM584v2_genomic.fna"
```

Example for *Drosophila* spike-in:

```yaml
run_contrasts: true
norm_method: "spikein"
spikein_genome: "drosophila"
spikein_reference:
  drosophila:
    fa: "/fdb/igenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa"
```

If spike-ins are unavailable or insufficient, normalization can alternatively be performed based on library size. Recommended workflow:

1. Run CARLISLE with `norm_method: spikein` for an initial QC assessment.
2. Evaluate spike-in alignment statistics.
3. Add `alignment_stats` to your configuration.
4. Re-run CARLISLE using library-size normalization.

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

### Quality Thresholds

Set peak-calling quality thresholds using the `quality_thresholds` parameter:

```yaml
quality_thresholds: "0.1, 0.05, 0.01"
```

Refer to tool-specific defaults:

* MACS2 q-value: [0.05](https://manpages.ubuntu.com/manpages/xenial/man1/macs2_callpeak.1.html)
* GoPeaks p-value: [0.05](https://github.com/maxsonBraunLab/gopeaks#usage)
* SEACR FDR threshold: [1.0](https://github.com/FredHutch/SEACR#usage)

---

## Reference Files

Additional reference genomes can be integrated by defining:

```yaml
species_name:
  fa: "/path/to/species.fa"
  blacklist: "/path/to/blacklistbed/species.bed"
  regions: "chr1 chr2 chr3"
  macs2_g: "hs" # genome shorthand for MACS2
```

> 🧭 **Best Practice:** Store reference paths under a centralized `/fdb` or `/data` location on Biowulf to ensure accessibility and consistency across users.

---

## Preparing Manifests

CARLISLE uses two manifests:

* `samplemanifest` – required for all analyses.
* `contrasts` – optional, required only for differential analysis with DESeq2.

### Sample Manifest (Required)

Defines sample-level metadata, including sample names, controls, and FASTQ paths.

| sampleName | replicateNumber | isControl | controlName                     | controlReplicateNumber | path_to_R1                                   | path_to_R2                                   |
| ---------- | --------------- | --------- | ------------------------------- | ---------------------- | -------------------------------------------- | -------------------------------------------- |
| 53_H3K4me3 | 1               | N         | HN6_IgG_rabbit_negative_control | 1                      | <path_to>/53_H3K4me3_1.R1.fastq.gz | <path_to>/53_H3K4me3_1.R2.fastq.gz |
| 54_H3K4me3 | 2               | N         | HN6_IgG_rabbit_negative_control | 1                      | <path_to>/54_H3K4me3_1.R1.fastq.gz | <path_to>/54_H3K4me3_1.R2.fastq.gz |
| HN6_IgG_rabbit_negative_control | 1               | Y         |  |                      | <path_to>/HN6_IgG_rabbit_negative_control_1.R1.fastq.gz | <path_to>/HN6_IgG_rabbit_negative_control_2.R2.fastq.gz |


### Contrast Manifest (Optional)

Specifies conditions for differential analysis:

| condition1              | condition2           |
| ----------------------- | -------------------- |
| MOC1_siSmyd3_2m_25_HCHO | MOC1_siNC_2m_25_HCHO |

> 📊 **Requirement:** Each condition must have at least two biological replicates to perform DESeq2-based differential analysis.
