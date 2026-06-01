# Running the Pipeline

## Pipeline Overview

The **CARLISLE** pipeline operates as a modular **[Snakemake](https://snakemake.readthedocs.io/en/stable/)** workflow, designed to support flexible execution on both local and cluster environments. It offers several run modes that control initialization, execution, and management of analysis sessions.

### Required Arguments

```bash
Usage: carlisle -m/--runmode=<RUNMODE> -w/--workdir=<WORKDIR>

1. RUNMODE [string]:
   init      – Initialize the working directory
   run       – Submit jobs to the SLURM cluster (Biowulf)
   reset     – Delete and reinitialize the working directory
   dryrun    – Validate and preview the workflow (no jobs executed)
   unlock    – Unlock the working directory after an aborted run
   runlocal  – Execute the pipeline interactively on a local node
   runtest   – Run the included test dataset on the cluster

2. WORKDIR [string]:
   Absolute or relative path to the desired output directory with write permissions.
```

### Optional Arguments

| Flag              | Description                                                                                                                             |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `--help, -h`      | Display the command-line help message.                                                                                                  |
| `--version, -v`   | Print the current version of CARLISLE.                                                                                                  |
| `--force, -f`     | Force re-execution of all Snakemake rules (overrides cache).                                                                            |
| `--singcache, -c` | Override the Singularity cache directory. See [Singularity Cache Directory](preparing-files.md#singularity-cache-directory) in Preparing Files for the full fallback order. |

### Scheduler Defaults

For cluster execution (`run`, `runtest`), CARLISLE uses scheduler-safe Snakemake defaults:

| Setting | Default | Override env var |
| ------- | ------- | ---------------- |
| Max concurrent jobs | `-j 100` | `CARLISLE_MAX_JOBS` |
| Max job submissions/sec | `--max-jobs-per-second 1` | `CARLISLE_MAX_JOBS_PER_SECOND` |
| Max status checks/sec | `--max-status-checks-per-second 0.1` | `CARLISLE_MAX_STATUS_CHECKS_PER_SECOND` |

---

## Command Descriptions

### Preparation Commands

- **`init` (required)** – Initializes the working directory by copying configuration, manifest, and Snakefiles into place. This step must be performed before any other pipeline action.

  - Use the `-f` or `--force` flag to reinitialize an existing directory.

- **`dryrun` (optional)** – Performs a non-executing validation of the Snakemake DAG, checking for syntax issues, missing files, or permission problems before a full run.

### Processing Commands

- **`runlocal`** – Executes the workflow on a local interactive node. This mode is suitable for quick testing or smaller datasets but should only be used within a Biowulf interactive session (`sinteractive`). A Snakemake HTML report (`report.html`) is generated in the working directory after a successful run.

- **`run`** – Submits the workflow to the **[Biowulf HPC cluster](https://hpc.nih.gov/)** via SLURM. CARLISLE manages job scheduling, dependencies, and notifications. Email alerts are automatically sent for job start, errors, and completion. The Singularity module is loaded automatically before submission — no manual `module load singularity` is required. A Snakemake HTML report (`report.html`) is generated in the working directory after a successful run.

### Maintenance Commands

- **`unlock`** – Unlocks the working directory if Snakemake terminates unexpectedly or a previous job is interrupted. **Use this when you see a "Directory is locked" error.** It is safe to run — it does not delete any files.

- **`reset`** – ⚠️ **Destructive.** Deletes the entire working directory and reinitializes it from scratch. All results, logs, and intermediate files are permanently lost. Only use this if you want to start completely over.

- **`runtest`** – Executes a small, bundled test dataset to verify installation and configuration integrity.

---

## Usage Syntax

All commands follow a consistent syntax:

```bash
carlisle --runmode=<COMMAND> --workdir=/path/to/output/dir
```

For example:

```bash
carlisle --runmode=init --workdir=/data/$USER/project
```

---

## Typical Workflow Example

A standard execution sequence on the Biowulf cluster would include the following steps:

```bash
# Step 1: Initialize working directory
carlisle --runmode=init --workdir=/path/to/output/dir
```

**Step 2: Edit your configuration files** — this is required before running:

| File | What to edit |
|---|---|
| `config/config.yaml` | Set `genome`, `norm_method`, `peaktype`, `run_contrasts`, and other parameters |
| `config/samples.tsv` | Fill in your sample names, replicate numbers, and FASTQ paths |
| `config/contrasts.tsv` | Fill in condition pairs (only needed if `run_contrasts: true`) |

```bash
# Step 3: Perform a dry run to validate before submitting
carlisle --runmode=dryrun --workdir=/path/to/output/dir
```

A successful dry run ends with a job summary table like this:

```
Job stats:
job                count
-----------------  -----
align              9
bam2bg             9
cov_correlation    1
...
total              NNN

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

If you see `Nothing to be done` it means all outputs already exist (re-run with `--force` to recompute). Any line starting with `MissingInputException` or `WorkflowError` indicates a configuration problem to fix before submitting.

```bash
# Step 4: Submit the full workflow to the cluster
carlisle --runmode=run --workdir=/path/to/output/dir
```

> ✅ **Recommendation:** Always perform a dry run before full execution to verify file paths, environment modules, and configuration correctness.

---

## Running in Control-Free Mode

If your experiment has no IgG or antibody control samples, enable control-free mode in `config/config.yaml` before running the pipeline:

```yaml
run_without_controls: true
seacr_threshold: 0.01
```

With this set:

- The sample manifest **does not** need `controlName` or `controlReplicateNumber` filled in — leave those columns blank for all treatment samples.
- `macs2_control` and `pool_controls` are automatically overridden by the pipeline; you do not need to set them manually.
- SEACR will use the numeric `seacr_threshold` (fraction of the signal distribution) instead of a control bedgraph.
- GoPeaks and MACS2 will run without a control BAM/fragment file.

The rest of the workflow is identical:

```bash
# Step 1: Initialize working directory
carlisle --runmode=init --workdir=/path/to/output/dir

# Step 2: Edit config/config.yaml — set run_without_controls: true

# Step 3: Dry run to validate
carlisle --runmode=dryrun --workdir=/path/to/output/dir

# Step 4: Submit to cluster
carlisle --runmode=run --workdir=/path/to/output/dir
```

> ⚠️ **Caution:** Control-free peak calling has higher false-positive rates. Validate results carefully, especially for SEACR where the numeric threshold is the sole background model.

---

## Monitoring a Running Job

After submitting with `run`, CARLISLE itself exits immediately — the pipeline runs as a background SLURM job. To monitor progress:

```bash
# Check your active SLURM jobs
squeue -u $USER

# Watch job status in real time (refreshes every 30 seconds)
watch -n 30 squeue -u $USER

# View the Snakemake master job log (replace JOBID with the number from squeue)
cat /path/to/output/dir/logs/snakemake.log
```

Email notifications are automatically sent to your NIH HPC account email (`$USER@nih.gov`) for:

- **Job start** — confirms the pipeline was accepted by SLURM
- **Job error** — sent if any rule fails; check the log file above for details
- **Job completion** — confirms all rules finished successfully

After a successful run, a `report.html` file is generated in your working directory — open it in a browser for an interactive summary of all pipeline steps and outputs.
