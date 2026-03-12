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
| `--singcache, -c` | Define a custom **Singularity cache directory**. Defaults to `/data/${USER}/.singularity`, or `${WORKDIR}/.singularity` if unavailable. |
Scheduler-safe defaults for cluster execution:
- `-j 100`
- `--max-jobs-per-second 1`
- `--max-status-checks-per-second 0.1`

Advanced users can override these defaults with environment variables:
- `CARLISLE_MAX_JOBS`
- `CARLISLE_MAX_JOBS_PER_SECOND`
- `CARLISLE_MAX_STATUS_CHECKS_PER_SECOND`

---

## Command Descriptions

### 🧩 Preparation Commands

- **`init` (required)** – Initializes the working directory by copying configuration, manifest, and Snakefiles into place. This step must be performed before any other pipeline action.

  - Use the `-f` or `--force` flag to reinitialize an existing directory.

- **`dryrun` (optional)** – Performs a non-executing validation of the Snakemake DAG, checking for syntax issues, missing files, or permission problems before a full run.

### ⚙️ Processing Commands

- **`runlocal`** – Executes the workflow on a local interactive node. This mode is suitable for quick testing or smaller datasets but should only be used within a Biowulf interactive session (`sinteractive`).

- **`run`** – Submits the workflow to the **[Biowulf HPC cluster](https://hpc.nih.gov/)** via SLURM. CARLISLE manages job scheduling, dependencies, and notifications. Email alerts are automatically sent for job start, errors, and completion.

### 🧰 Maintenance Commands

- **`unlock`** – Unlocks the working directory if Snakemake terminates unexpectedly or a previous job is interrupted.

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

# Step 2: Perform a dry run to validate the workflow
carlisle --runmode=dryrun --workdir=/path/to/output/dir

# Step 3: Submit the full workflow to the cluster
carlisle --runmode=run --workdir=/path/to/output/dir
```

> ✅ **Recommendation:** Always perform a dry run before full execution to verify file paths, environment modules, and configuration correctness.
