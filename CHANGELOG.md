## CARLISLE development version

### New Features

- **New helper script `_get_pooled_scale.py`**: Computes pooled control scaling factors using named column access from `alignment_stats.tsv`, replacing fragile positional awk column references in `create_pooled_control_bedgraph`. (#230)
- **Control-free analysis mode**: Support for running peak calling analysis without control samples for MACS2, SEACR, and GoPeaks. Enabled via `run_without_controls: true` config flag with `seacr_threshold` parameter for SEACR numeric threshold specification. (#224, #225, #226)
- **Documentation versioning**: Added mike plugin for version-specific documentation selector with dropdown menu on ReadTheDocs. Automatically injects version selector into all documentation pages.
- **Singularity cache configuration**: Explicit `SIFCACHE` environment variable support throughout wrapper script and Snakemake configuration for flexible container storage location management.

### Bug Fixes

- **Pooled control bedgraph: reversed scaling factor columns**: `create_pooled_control_bedgraph` used `$NF` (col 11, `dedup_nreads_spikein`) for LIBRARY and `$(NF-1)` (col 10, `dedup_nreads_genome`) for SPIKEIN — exactly swapped. Fixed with named column access via `_get_pooled_scale.py`. (#230)
- **Pooled control bedgraph: SPIKEIN used pre-filter read counts**: Corrected to `no_dedup_nreads_spikein` (fragment-length + mapq filtered, not deduplicated) to match the individual-replicate `bam2bg` rule. (#230)
- **Pooled control bedgraph: averaged instead of summed read counts**: The merged pooled BAM contains reads from all replicates; its total depth equals the sum of individual counts, not their average. Fixed aggregation. (#230)
- **Pooled control bedgraph: spurious dupstatus filter**: Old awk filtered `$2 == "{dupstatus}"` on a non-existent column, silently dropping all rows. `alignment_stats.tsv` has one row per replicate. Filter removed. (#230)
- **Pooled control outputs stale on norm_method change**: Fragment and bedgraph output paths for pooled controls contained no reference to the normalization method, causing Snakemake to silently skip regeneration when `norm_method` changed in an existing results directory. SEACR received a control bedgraph from the previous normalization run with no warning. Fixed by embedding `NORM_METHOD` in filenames for `create_pooled_control_fragments`, `create_pooled_control_bedgraph`, and all downstream consumers (`macs2_narrow`, `macs2_broad`, `seacr_stringent`, `seacr_relaxed`). GoPeaks rules are unaffected as they consume merged BAM directly. (#232)
- **Unified bigwig generation via `bamCoverage --scaleFactor`**: Replaced `bedGraphToBigWig` in `bam2bg` and removed the separate `deeptools_bw` rule (which used RPGC normalization). Both browser-track bigwigs and deeptools heatmap/profile inputs now come from a single `bamCoverage --scaleFactor` call, ensuring heatmaps/profiles reflect the same spike-in/library normalization as genome browser tracks. `deeptools_prep` and `deeptools_mat` updated to reference `bigwig/` instead of `deeptools/temp/`. (#231)
- **deeptools_prep control-free handling**: Fixed rule to use `TREATMENT_WITHOUTCONTROL_LIST` in control-free mode instead of empty `TREATMENT_CONTROL_LIST`, preventing MissingOutputException. Properly excludes 'nocontrol' sentinel from bigwig file lists.
- **Differential analysis with nocontrol pairs**: Fixed `_control_has_replicate()` function in diff.smk to properly handle control-free sentinel values, ensuring control-free pairs are included in individual mode analysis.
- **Singularity module availability on login node**: Fixed runslurm() to load Singularity module before dryrun execution, preventing WorkflowError when submitting to cluster scheduler.
- **Singularity prefix configuration consistency**: Added explicit `--singularity-prefix` to all Snakemake invocations to respect SIFCACHE environment variable settings across all execution modes.
- **Enhanced error messaging**: Improved error messages for file validation and control sample checks in init.smk for clearer troubleshooting.
- **Version and module diagnostics**: Added print_versions() function to display Snakemake and Singularity versions and added module loading diagnostics to SBATCH scripts for better troubleshooting of cluster execution issues.
- **ROSE container: truncated prep script in v1**: `Dockerfile.rose` used an incorrect `COPY` source path (`_prep_rose_input.py` at build-context root) instead of `workflow/scripts/_prep_rose_input.py`, resulting in an empty/truncated file being baked into the image. At runtime this caused a `SyntaxError: EOL while scanning string literal` in `/opt/ROSE/_prep_rose_input.py`, failing all ROSE jobs. Fixed in container `ccbr_rose:v2`.
- **Peak count outputs overwritten between control/no-control runs**: `all.peaks.txt` and `Peak counts.xlsx` had fixed names, so running with `run_without_controls: true` after a `false` run (or vice versa) in the same results directory would silently overwrite the previous aggregated counts. Output filenames now encode the run mode: `all.peaks.with_control.txt` / `all.peaks.without_control.txt` and `Peak counts.with_control.xlsx` / `Peak counts.without_control.xlsx`. `_plot_peak_counts.R` updated to accept the xlsx filename as an optional 3rd argument.
- **ROSE control BAM passed for nocontrol sentinel**: In control-free mode (`run_without_controls: true`), the `rose` rule incorrectly constructed `--control-bam nocontrol.dedup.bam` for the `run-prep-rose` call, causing an immediate failure (`[ERROR] control does not exist`). Added an explicit check for the `nocontrol` sentinel so ROSE runs without `-c` in control-free mode.
- **Peak counts xlsx not written in control-free mode**: The `count_peaks` rule passed the xlsx filename as a shell-quoted string argument to `_plot_peak_counts.R`, but the space in the filename (`Peak counts.without_control.xlsx`) caused the argument to be split by bash when invoked through Singularity, resulting in R falling back to the default filename `Peak counts.xlsx`. Fixed by passing `{output.peak_table:q}` (Snakemake shell quoting) and having the R script use the received argument as a full absolute path.

## CARLISLE 2.7.6

### Improvements

- **Scheduler-safe Snakemake defaults**: Replaced hardcoded `-j 500` with safer defaults (`-j 100`, `--max-jobs-per-second 1`, `--max-status-checks-per-second 0.1`) for cluster-friendly submission and status polling behavior. (@kopardev)

## CARLISLE 2.7.5

### New Features

- **ROSE containerization**: Containerized the ROSE workflow, added a dedicated prep script, and simplified dependencies by removing annotation-folder/refseq coupling for supported genomes (`hg19`, `hg38`, `mm10`). (#215, @kopardev)
- **GO enrichment workflow split**: Separated GO enrichment table generation from dotplot generation to improve rerun behavior and failure isolation. (#215, @kopardev)

### Improvements

- **ROSE output streamlining**: Reduced ROSE outputs to only required deliverables and adjusted cluster resource requests accordingly. (#215, @kopardev)
- **GO enrichment execution hardening**: Improved GO enrichment and dotplot logging and scheduling defaults for cluster execution. (#210, #211, #212, #213, #215, @kopardev)

### Bug Fixes

- **GO enrichment robustness**: Handle empty BED/TSV inputs without hard failure and improve fallback checks in dotplot generation. (#212, #215, @kopardev)
- **Dotplot label handling**: Fix duplicate wrapped enrichment labels in GO dotplot output. (#215, @kopardev)
- **ROSE empty-input handling**: Prevent hard failures when ROSE prep receives empty peak inputs. (#215, @kopardev)

## CARLISLE 2.7.4

### Bug Fixes

- **HOMER annotation outputs with control modes**: Fix `rule all` expectations to include the `control_mode` subdirectory for HOMER annotation plots and combined q-value tables, preventing MissingInputException during dryruns. (#209, @kopardev)

## CARLISLE 2.7.3

### New Features

- **Pooled control mode support**: Complete pipeline implementation for pooled control analysis across MACS2, SEACR, GoPeaks, ROSE, HOMER, and differential analysis workflows. Enables comparison of replicate-specific vs merged high-depth controls. (#206, @kopardev)
- **ROSE dual control mode**: ROSE enhancer analysis now runs separately for both individual and pooled control modes, generating separate super-enhancer calls for each approach. (#206, @kopardev)
- **Reference file compression**: All reference BED files (blacklists, TSS, gene annotations, cCREs) now stored as `.bed.gz` with automatic decompression during analysis, significantly reducing storage requirements. (#206, @kopardev)
- **cCRE annotations**: Added comprehensive candidate cis-Regulatory Element (cCRE) annotations from ENCODE SCREEN database for all supported genomes (hg38, hg19, mm10, hs1) (#206, @kopardev), including:
  - Promoter-like signatures (PLS)
  - Proximal enhancer-like signatures (pELS)
  - Distal enhancer-like signatures (dELS)
  - Chromatin accessibility regions (CA-CTCF, CA-H3K4me3, CA-TF)
- **Motif enrichment for DEG peaks**: Added HOMER motif discovery and AME motif enrichment analysis specifically for differentially enriched peaks (both AUC-based and fragments-based, for up-regulated peaks in each group). (#206, @kopardev)
- **Enhanced differential analysis outputs**: Added 3-column BED files for up-regulated peaks in each group (`up_group1.bed`, `up_group2.bed`) for downstream enrichment analyses. (#206, @kopardev)
- **HOCOMOCO v14 CORE motifs**: Added complete HOCOMOCO v14 CORE motif database in both HOMER and MEME formats for comprehensive motif enrichment analysis. (#206, @kopardev)

### Improvements

- **Increased resource allocations**: Enhanced memory and thread allocation for computationally intensive rules (#206, @kopardev):
  - ROSE: 96GB memory, 16 threads
  - DESeq2: 96GB memory
  - GO enrichment: 32GB memory, 8 threads
  - deepTools tasks: Increased memory allocations
- **Enhanced GO enrichment**: Improved GO enrichment workflow with contrast file support and parallel processing options for better performance. (#206, @kopardev)
- **deepTools optimization**: Added temporary directory handling and expanded bedtype support for improved coverage analysis. (#206, @kopardev)
- **Improved logging and reruns**: Better handling of rerun scenarios and enhanced logging throughout the pipeline. (#206, @kopardev)

### Bug Fixes

- **ROSE environment isolation**: Fixed Python library conflicts between Snakemake and ROSE environments by explicitly managing `PYTHONPATH` and unsetting conda variables. (#206, @kopardev)
- **ROSE chromosome filtering**: Properly filter NC\_ chromosomes (alternative scaffolds, unplaced contigs) from both treatment and control BAM files before enhancer stitching to prevent analysis failures. (#206, @kopardev)
- **BED file decompression**: Implemented consistent decompression handling for all compressed BED files across init, alignment, and annotation rules. (#206, @kopardev)

## CARLISLE 2.7.2

- Fix how singularity bind paths are set. (#187, @kelly-sovacool)
- Fix the DESeq rule when using the hs1 T2T genome. (#187, @kelly-sovacool, @wong-nw)
- Fix incorrect jobby version (needs v0.4). (#194, @kelly-sovacool)
- Fix jobby usage to output stdout/stderr to log file. (#196, @kopardev)
- Fix numpy error. (#198, @kelly-sovacool)
- Fix logic error with scaling factor calcutions. (#195, @kopardev)
- Run MACS2 with control by default. (#195, @kopardev)
- Make sure `submit_slurm.sbatch` is not overwritten. (#195, @kopardev)
- Handle case when sample names are subsets of each other. (#195, @kopardev)

## CARLISLE 2.7.1

- Remove deprecated 'ccr' partition. (#176, @kopardev)

## CARLISLE 2.7.0

- Now depends on ccbr_tools v0.4 for updated jobby & spooker utilities. (#168, @kelly-sovacool)

## CARLISLE 2.6.4

- Fix driver script: do not try to load the shared conda environment. (#166, @kelly-sovacool)
- Now depends on ccbr_tools v0.4 for updated jobby & spooker utilities. (#168, @kelly-sovacool)

## CARLISLE 2.6.3

- Minor documentation update to use readthedocs theme. (#162, @kelly-sovacool)
- Fix deprecation warning during initialization. (#163, @kelly-sovacool)

## CARLISLE 2.6.2

- Documentation improvements. (#154, @kelly-sovacool)
- If jobby & spooker are not available, try adding them to the path on workflow completion. (#155, @kelly-sovacool)
- Fix bug that flipped library normalization scaling factor (#157, @epehrsson)

## CARLISLE 2.6.1

- Load the module for snakemake v7, but do not specify the minor and patch versions. (#149, @kelly-sovacool)

## CARLISLE 2.6.0

### Bug fixes

- Bug fixes for DESeq (#127, @epehrsson)
  - Removes single-sample group check for DESeq.
  - Increases memory for DESeq.
  - Ensures control replicate number is an integer.
  - Fixes FDR cutoff misassigned to log2FC cutoff.
  - Fixes `no_dedup` variable names in library normalization scripts.
- Fig bug that added nonexistent directories to the singularity bind paths. (#135, @kelly-sovacool)
- Containerize rules that require R (`deseq`, `go_enrichment`, and `spikein_assessment`) to fix installation issues with common R library path. (#129, @kelly-sovacool)
  - The `Rlib_dir` and `Rpkg_config` config options have been removed as they are no longer needed.

### New features

- New visualizations: (#132, @epehrsson)
  - New rules `cov_correlation`, `homer_enrich`, `combine_homer`, `count_peaks`
  - Add peak caller to MACS2 peak xls filename
- New parameters in the config file to make certain rules optional: (#133, @kelly-sovacool)
  - GO enrichment is controlled by `run_go_enrichment` (default: `false`)
  - ROSE is controlled by `run_rose` (default: `false`)
- New `--singcache` argument to provide a singularity cache dir location. The singularity cache dir is automatically set inside `/data/$USER/` or `$WORKDIR/` if `--singcache` is not provided. (#143, @kelly-sovacool)

### Misc

- The singularity version is no longer specified, per request of the biowulf admins. (#139, @kelly-sovacool)
- Minor documentation updates. (#146, @kelly-sovacool)

## CARLISLE v2.5.0

- Refactors R packages to a common source location (#118, @slsevilla)
- Adds a --force flag to allow for re-initialization of a workdir (#97, @slsevilla)
- Fixes error with testrun in DESEQ2 (#113, @slsevilla)
- Decreases the number of samples being run with testrun, essentially running tinytest as default and removing tinytest as an option (#115, @slsevilla)
- Reads version from VERSION file instead of github repo link (#96, #112, @slsevilla)
- Added a CHANGELOG (#116, @slsevilla)
- Fix: RNA report bug, caused by hard-coding of PC1-3, when only PC1-2 were generated (#104, @slsevilla)
- Minor documentation improvements. (#100, @kelly-sovacool)
- Fix: allow printing the version or help message even if singularity is not in the path. (#110, @kelly-sovacool)

## CARLISLE v2.4.1

- Add GitHub Action to add issues/PRs to personal project boards by @kelly-sovacool in #95
- Create install script by @kelly-sovacool in #93
- feat: use summits bed for homer input; save temporary files; fix deseq2 bug by @slsevilla in #108
- docs: adding citation and DOI to pipeline by @slsevilla in #107
- Test a dryrun with GitHub Actions by @kelly-sovacool in #94

## CARLISLE v2.4.0

- Feature- Merged Yue's fork, adding DEEPTOOLS by @slsevilla in #85
- Feature- Added tracking features from SPOOK by @slsevilla in #88
- Feature - Dev test run completed by @slsevilla in #89
- Bug - Fixed bugs related to Biowulf transition

## CARLISLE v2.1.0

- enhancement
- update gopeaks resources
- change SEACR to run "norm" without spikein controls, "non" with spikein controls
- update docs for changes; provide extra troubleshooting guidance
- fix GoEnrich bug for failed plots

## CARLISLE v2.0.1

- fix error when contrasts set to "N"
- adjust goenrich resources to be more efficient

## CARLISLE 2.0.0

- Add a MAPQ filter to samtools (rule align)
- Add GoPeaks MultiQC module
- Allow for library normalization to occur during first pass
- Add --broad-cutoff to MACS2 broad peak calling for MACS2
- Create a spike in QC report
- Reorganize file structure to help with qthreshold folder
- Update variable names of all peak caller
- Merge rules with input/output/wildcard congruency
- Convert the "spiked" variable to "norm_method
- Add name of control used to MACS2 peaks
- Running extra control:sample comparisons that are not needed
- improved resource allocation
- test data originally included 1475 jobs, this version includes 1087 jobs (reduction of 25%) despite including additional features
- moved ~12% of all jobs to local deployment (within SLURM submission)

## CARLISLE 1.2.0

- merge increases to resources; update workflow img, contributions

## CARLISLE 1.1.1

- patch for gz bigbed bug

## CARLISLE 1.1.0

- add broad-cutoff to macs2 broad peaks param settings
- add non.stringent and non.relaxed to annotation options
- merge DESEQ and DESEQ2 rules together
- identify some files as temp

## CARLISLE 1.0.1

- contains patch for DESEQ error with non hs1 reference samples
