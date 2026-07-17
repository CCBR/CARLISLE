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

| Flag              | Description                                                                                                                                                                 |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--help, -h`      | Display the command-line help message.                                                                                                                                      |
| `--version, -v`   | Print the current version of CARLISLE.                                                                                                                                      |
| `--force, -f`     | Force re-execution of all Snakemake rules (overrides cache).                                                                                                                |
| `--singcache, -c` | Override the Singularity cache directory. See [Singularity Cache Directory](preparing-files.md#singularity-cache-directory) in Preparing Files for the full fallback order. |

### Scheduler Defaults

For cluster execution (`run`, `runtest`), CARLISLE uses scheduler-safe Snakemake defaults:

| Setting                 | Default                              | Override env var                        |
| ----------------------- | ------------------------------------ | --------------------------------------- |
| Max concurrent jobs     | `-j 100`                             | `CARLISLE_MAX_JOBS`                     |
| Max job submissions/sec | `--max-jobs-per-second 1`            | `CARLISLE_MAX_JOBS_PER_SECOND`          |
| Max status checks/sec   | `--max-status-checks-per-second 0.1` | `CARLISLE_MAX_STATUS_CHECKS_PER_SECOND` |

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

| File                   | What to edit                                                                   |
| ---------------------- | ------------------------------------------------------------------------------ |
| `config/config.yaml`   | Set `genome`, `norm_method`, `peaktype`, `run_contrasts`, and other parameters |
| `config/samples.tsv`   | Fill in your sample names, replicate numbers, and FASTQ paths                  |
| `config/contrasts.tsv` | Fill in condition pairs (only needed if `run_contrasts: true`)                 |

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
quality_thresholds: "0.01"
```

With this set:

- The sample manifest **does not** need `controlName` or `controlReplicateNumber` filled in — leave those columns blank for all treatment samples.
- `macs2_control` and `pool_controls` are automatically overridden by the pipeline; you do not need to set them manually.
- SEACR will use each value in `quality_thresholds` as the numeric threshold (fraction of the signal distribution) instead of a control bedgraph.
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

## Understanding Log Output

CARLISLE prints structured log lines to the terminal during `init`, `dryrun`, `run`, `runlocal`, and other modes. Each line is prefixed with a fixed-width tag that indicates its severity or purpose:

| Prefix | Meaning | When you see it |
|--------|---------|-----------------|
| `STEP` | A major pipeline phase is starting | Beginning of `init`, `dryrun`, `run`, Snakemake submission, etc. |
| `INFO` | Informational detail about the current step | Paths, settings, versions, tool availability |
| `OK` | A step completed successfully | Files copied, modules loaded, job submitted without error |
| `WARN` | A non-fatal issue that may need attention | Outdated submit script backed up, optional tool unavailable |
| `ERROR` | A fatal problem — execution will stop | Missing files, failed module load, bad arguments |
| `NEXT` | Suggested next action for the user | What to edit or run after the current step completes |

### Example terminal output

A typical `carlisle --runmode=run` invocation looks like this:

```
------------------------------------------------------------------
STEP  [run] Preparing SLURM submission
INFO  CARLISLE Run Summary
INFO  Mode:       SLURM
INFO  Workdir:    /data/$USER/project
INFO  Partition:  norm
INFO  Max jobs:   100
INFO  Scheduler:  --max-jobs-per-second 1, --max-status-checks-per-second 0.1
INFO  Log file:   /data/$USER/project/logs/snakemake.log
------------------------------------------------------------------
INFO  Tool Versions:
7.32.4
OK    Snakemake version checked
apptainer version 1.3.6-1.bionic
OK    Singularity version checked
------------------------------------------------------------------
OK    Dry-run was successful. Submitting jobs to scheduler.
------------------------------------------------------------------
OK    Job submitted successfully (SLURM job ID: 12345678)
INFO  Submission output: 12345678
NEXT  Monitor:    squeue -u $USER
NEXT  Progress:   tail -f /data/$USER/project/logs/snakemake.log
NEXT  Status:     ls -1 /data/$USER/project/pipeline.*
NEXT  Sidecar:    /data/$USER/project/pipeline.status.json
------------------------------------------------------------------
```

`NEXT` lines tell you exactly what to run after each stage — copy them directly into your terminal.

### Live progress in `pipeline.running`

While a `run` job is active, CARLISLE updates the `pipeline.running` marker file every 60 seconds with a human-readable progress summary parsed from `snakemake.log`:

```
Progress : 47 / 312 steps complete (15%)
Remaining: 265 steps
Updated  : 2026-07-17 14:32:10
```

Read it at any time with:

```bash
cat /path/to/output/dir/pipeline.running
```

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

# Check CARLISLE run-state marker (exactly one should exist)
ls -1 /path/to/output/dir/pipeline.*

# Inspect state metadata
cat /path/to/output/dir/pipeline.status.json
```

`pipeline.status.json` contains the following fields:

| Field | Description |
|---|---|
| `pipeline` | Pipeline name (`CARLISLE`) |
| `version` | CARLISLE version at the time of submission |
| `state` | Current state: `running`, `completed`, `failed`, or `canceled` |
| `reason` | Machine-readable reason for the state (e.g. `snakemake_and_report_succeeded`, `snakemake_failed`) |
| `runmode` | Run mode used (e.g. `run`) |
| `workdir` | Absolute path to the working directory |
| `user` | Username that submitted the pipeline |
| `slurm_job_id` | SLURM job ID of the pipeline coordinator job |
| `host` | Hostname where the state was last written |
| `submission_timestamp_utc` | UTC timestamp when `carlisle --runmode=run` was invoked |
| `start_timestamp_utc` | UTC timestamp when the SLURM job began executing |
| `duration_seconds` | Wall-clock seconds from job start to final state (null while running or on headnode-only writes) |
| `exit_code` | Exit code of Snakemake (0 = success; null for headnode-only writes) |
| `tasks_done` | Number of Snakemake steps completed at final state (null if log is unavailable) |
| `tasks_total` | Total number of Snakemake steps at final state (null if log is unavailable) |
| `snakemake_log` | Absolute path to the Snakemake master log file |
| `timestamp_utc` | UTC timestamp of the most recent state write |

For `runmode=run`, CARLISLE writes exactly one state marker file in the workdir:

- `pipeline.running` — run submitted and in progress
- `pipeline.completed` — run finished successfully
- `pipeline.failed` — submission or runtime failure
- `pipeline.canceled` — run interrupted (for example `scancel`/signal)

When you start `runmode=run` again in the same workdir, CARLISLE removes any existing `pipeline.*` marker and replaces it with `pipeline.running` for the new run.

Email notifications are automatically sent to your NIH HPC account email (`$USER@nih.gov`) for:

- **Job start** — confirms the pipeline was accepted by SLURM
- **Job error** — sent if any rule fails; check the log file above for details
- **Job completion** — confirms all rules finished successfully

After a successful run, a `report.html` file is generated in your working directory — open it in a browser for an interactive summary of all pipeline steps and outputs.
