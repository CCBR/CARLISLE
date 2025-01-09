# Running the Pipeline

## Pipeline Overview

The Snakemake workflow has a multiple options

### Required arguments

```
Usage: carlisle -m/--runmode=<RUNMODE> -w/--workdir=<WORKDIR>

1.  RUNMODE: [Type: String] Valid options:
    *) init : initialize workdir
    *) run : run with slurm
    *) reset : DELETE workdir dir and re-init it
    *) dryrun : dry run snakemake to generate DAG
    *) unlock : unlock workdir if locked by snakemake
    *) runlocal : run without submitting to sbatch
    *) runtest: run on cluster with included test dataset
2.  WORKDIR: [Type: String]: Absolute or relative path to the output folder with write permissions.
```

### Optional arguments

--help|-h : print this help.
--version|-v : print the version of carlisle.
--force|-f : use the force flag for snakemake to force all rules to run.
--singcache|-c : singularity cache directory. Default is `/data/${USER}/.singularity` if available, or falls back to `${WORKDIR}/.singularity`. Use this flag to specify a different singularity cache directory.

## Commands explained

The following explains each of the command options:

- Preparation Commands
  - init (REQUIRED): This must be performed before any Snakemake run (dry, local, cluster) can be performed. This will copy the necessary config, manifest and Snakefiles needed to run the pipeline to the provided output directory.
    - the -f/--force flag can be used in order to re-initialize a workdir that has already been created
  - dryrun (OPTIONAL): This is an optional step, to be performed before any Snakemake run (local, cluster). This will check for errors within the pipeline, and ensure that you have read/write access to the files needed to run the full pipeline.
- Processing Commands
  - local: This will run the pipeline on a local node. NOTE: This should only be performed on an interactive node.
  - run: This will submit a master job to the cluster, and subsequent sub-jobs as needed to complete the workflow. An email will be sent when the pipeline begins, if there are any errors, and when it completes.
- Other Commands (All optional)
  - unlock: This will unlock the pipeline if an error caused it to stop in the middle of a run.
  - runtest: This will run a test of the pipeline with test data

To run any of these commands, follow the the syntax:

```
carlisle --runmode=COMMAND --workdir=/path/to/output/dir
```

## Typical Workflow

A typical command workflow, running on the cluster, is as follows:

```
carlisle --runmode=init --workdir=/path/to/output/dir

carlisle --runmode=dryrun --workdir=/path/to/output/dir

carlisle --runmode=run --workdir=/path/to/output/dir
```
