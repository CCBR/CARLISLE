# 2. Preparing Files
The pipeline is controlled through editing configuration and manifest files. Defaults are found in the /WORKDIR/config and /WORKDIR/manifest directories, after initialization.

## 2.1 Configs
The configuration files control parameters and software of the pipeline. These files are listed below:

- config/config.yaml
- resources/cluster.yaml
- resources/tools.yaml

### 2.1.1 Cluster Config
The cluster configuration file dictates the resouces to be used during submission to Biowulf HPC. There are two differnt ways to control these parameters - first, to control the default settings, and second, to create or edit individual rules. These parameters should be edited with caution, after significant testing.

### 2.1.2 Tools Config
The tools configuration file dictates the version of each software or program that is being used in the pipeline.

### 2.1.3 Config YAML
There are several groups of parameters that are editable for the user to control the various aspects of the pipeline. These are :

- Folders and Paths
    - These parameters will include the input and ouput files of the pipeline, as well as list all manifest names.
- User parameters
    - These parameters will control the pipeline features. These include thresholds and whether to perform processes.
- References
    - These parameters will control the location of index files, spike-in references, adaptors and species calling information.

#### 2.1.3.1 User Parameters 
##### 2.1.3.1.1 (Spike in Controls)
The pipeline allows for the use of a species specific spike-in control, or the use of normalization via library size. The parameter `spikein_genome` should be set to the species term used in `spikein_reference`.

For example for ecoli spike-in:
```
run_contrasts: "Y"
norm_method: "spikein"
spikein_genome: "ecoli"
spikein_reference:
  ecoli:
    fa: "PIPELINE_HOME/resources/spikein/Ecoli_GCF_000005845.2_ASM584v2_genomic.fna"

```

For example for drosophila spike-in:
```
run_contrasts: "Y"
norm_method: "spikein"
spikein_genome: "drosophila"
spikein_reference:
  drosophila:
    fa: "/fdb/igenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa"

```

If it's determined that the amount of spike-in is not sufficient for the run, a library normaliaztion can be performed.
1. Complete a CARLISLE run with spike-in set to "Y". This will allow for the complete assessment of the spike-in.
2. Run inital QC analysis on the output data
3. Add the alignment_stats dir to the configuration file.
4. Re-run the CARLISLE pipeline 

##### 2.1.3.1.2 Duplication Status
Users can select duplicated peaks (dedup) or non-deduplicated peaks (no_dedup) through the user parameter.
```
dupstatus: "dedup, no_dedup" 
```

##### 2.1.3.1.3 Peak Caller
Three peak callers are available for deployment within the pipeline, with different settings deployed for each caller.

1. [MACS2](https://github.com/macs3-project/MACS) is available with two peak calling options: narrowPeak or broadPeak. NOTE: DESeq step generally fails for broadPeak; generally has too many calls.
```
peaktype: "macs2_narrow, macs2_broad,"
```
2. [SEACR](https://github.com/FredHutch/SEACR) is available with four peak calling options: stringent or relaxed parameters, to be paired with "norm" for samples without a spike-in control and "non" for samples with a spikein control
```
peaktype: "seacr_stringent, seacr_relaxed"
```
3. [GOPEAKS](https://github.com/maxsonBraunLab/gopeaks) is available with two peak calling options: narrowpeaks or broadpeaks
```
peaktype: "gopeaks_narrow, gopeaks_broad"
```
A complete list of the available peak calling parameters and the recommended list of parameters is provided below:

| Peak Caller | Narrow              | Broad             | Normalized, Stringent | Normalized, Relaxed   | Non-Normalized, Stringent | Non-Normalized, Relaxed |
| --- | --- | --- |--- |--- |--- |--- |
| Macs2 | AVAIL | AVAIL | NA | NA | NA | NA |
| SEACR | NA | NA | AVAIL w/o SPIKEIN | AVAIL w/o SPIKEIN | AVAIL w/ SPIKEIN | AVAIL w/ SPIKEIN |
| GoPeaks | AVAIL | AVAIL | NA | NA | NA | NA |

```
# Recommended list
### peaktype: "macs2_narrow, macs2_broad, gopeaks_narrow, gopeaks_broad"

# Available list
### peaktype: "macs2_narrow, macs2_broad, seacr_norm_stringent, seacr_norm_relaxed, seacr_non_stringent, seacr_non_relaxed, gopeaks_narrow, gopeaks_broad"
```

##### 2.1.3.1.3.1 Macs2 additional option
MACS2 can be run with or without the control. adding a control will increase peak specificity
Selecting "Y" for the `macs2_control` will run the paired control sample provided in the sample manifest

##### 2.1.3.1.4 Quality Tresholds
Thresholds for quality can be controled through the `quality_tresholds` parameter. This must be a list of comma separated values. minimum of numeric value required.

- default MACS2 qvalue is 0.05 https://manpages.ubuntu.com/manpages/xenial/man1/macs2_callpeak.1.html
- default GOPEAKS pvalue is 0.05 https://github.com/maxsonBraunLab/gopeaks/blob/main/README.md
- default SEACR FDR threshold 1 https://github.com/FredHutch/SEACR/blob/master/README.md
```
#default values
quality_thresholds: "0.1, 0.05, 0.01"
```

#### 2.1.3.2 References
Additional reference files may be added to the pipeline, if other species were to be used. 

The absolute file paths which must be included are:

1. fa: "/path/to/species.fa"
2. blacklist: "/path/to/blacklistbed/species.bed"

The following information must be included:

1. regions: "list of regions to be included; IE chr1 chr2 chr3"
2.  macs2_g: "macs2 genome shorthand; IE mm IE hs"

## 2.2 Preparing Manifests
There are two manifests, one which required for all pipeliens and one that is only required if running a differential analysis. These files describe information on the samples and desired contrasts. The paths of these files are defined in the snakemake_config.yaml file. These files are:

- samplemanifest
- contrasts

### 2.2.1 Samples Manifest (REQUIRED)
This manifest will include information to sample level information. It includes the following column headers:

- sampleName: the sample name WITHOUT replicate number (IE "SAMPLE")
- replicateNumber: the sample replicate number (IE "1")
- isControl: whether the sample should be identified as a control (IE "Y")
- controlName: the name of the control to use for this sample (IE "CONTROL")
- controlReplicateNumber: the replicate number of the control to use for this sample (IE "1")
- path_to_R1: the full path to R1 fastq file (IE "/path/to/sample1.R1.fastq")
- path_to_R2: the full path to R1 fastq file (IE "/path/to/sample2.R2.fastq")

An example sampleManifest file is shown below:


| sampleName| replicateNumber| isControl| controlName| controlReplicateNumber| path_to_R1| path_to_R2
| --- |--- |--- |--- |--- |--- |--- |
| 53_H3K4me3| 1| N| HN6_IgG_rabbit_negative_control| 1| PIPELINE_HOME/.test/53_H3K4me3_1.R1.fastq.gz| PIPELINE_HOME/.test/53_H3K4me3_1.R2.fastq.gz
| 53_H3K4me3| 2| N| HN6_IgG_rabbit_negative_control| 1| PIPELINE_HOME/.test/53_H3K4me3_2.R1.fastq.gz| PIPELINE_HOME/.test/53_H3K4me3_2.R2.fastq.gz
| HN6_H3K4me3| 1| N| HN6_IgG_rabbit_negative_control| 1| PIPELINE_HOME/.test/HN6_H3K4me3_1.R1.fastq.gz| PIPELINE_HOME/.test/HN6_H3K4me3_1.R2.fastq.gz
| HN6_H3K4me3| 2| N| HN6_IgG_rabbit_negative_control| 1| PIPELINE_HOME/.test/HN6_H3K4me3_2.R1.fastq.gz| PIPELINE_HOME/.test/HN6_H3K4me3_2.R2.fastq.gz
| HN6_IgG_rabbit_negative_control| 1| Y| -| -| PIPELINE_HOME/.test/HN6_IgG_rabbit_negative_control_1.R1.fastq.gz| PIPELINE_HOME/.test/HN6_IgG_rabbit_negative_control_1.R2.fastq.gz


### 2.2.2 Contrast Manifest (OPTIONAL)
This manifest will include sample information to performed differential comparisons.

An example contrast file:

| condition1 | condition2 |
| --- | --- |
| MOC1_siSmyd3_2m_25_HCHO | MOC1_siNC_2m_25_HCHO |
