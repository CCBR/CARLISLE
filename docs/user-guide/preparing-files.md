# 2. Preparing Files
The pipeline is controlled through editing configuration and manifest files. Defaults are found in the /output/dir/config and /output/dir/manifest directories, after initialization.

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
- Advanced parameters
  - These parameters generally do not need to be edited, as they control the outputs to peak calling and duplication type.
- References
  - These parameters will control the location of index files, spike-in references, adaptors and species calling information.

The pipeline allows for the use of a species specific spike-in control, or the use of normalization via library size. The parameter `spikein_genome` should be set to the species term used in `spikein_reference`.

For example for ecoli spike-in:
```
run_contrasts: "Y"

spikein_genome: "ecoli"

spikein_reference:
  ecoli:
    fa: "PIPELINE_HOME/resources/spikein/Ecoli_GCF_000005845.2_ASM584v2_genomic.fna"

```

For example for drosophila spike-in:
```
run_contrasts: "Y"

spikein_genome: "drosophila"

spikein_reference:
  drosophila:
    fa: "/fdb/igenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa"

```

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

```
sampleName	replicateNumber	isControl	controlName	controlReplicateNumber	path_to_R1	path_to_R2
MOC1_siNC_2m_25_HCHO	1	N	MOC1_siNC_2m_25_HCHO_IgG_Ctrl	1	PIPELINE_HOME/.test/MOC1_siNC_2m_25_HCHO_1.1m.R1.fastq.gz	PIPELINE_HOME/.test/MOC1_siNC_2m_25_HCHO_1.1m.R2.fastq.gz
MOC1_siSmyd3_2m_25_HCHO	1	N	MOC1_siNC_2m_25_HCHO_IgG_Ctrl	1	PIPELINE_HOME/.test/MOC1_siSmyd3_2m_25_HCHO_1.1m.R1.fastq.gz	PIPELINE_HOME/.test/MOC1_siSmyd3_2m_25_HCHO_1.1m.R2.fastq.gz
MOC1_siSmyd3_2m_25_HCHO	2	N	MOC1_siNC_2m_25_HCHO_IgG_Ctrl	1	PIPELINE_HOME/.test/MOC1_siSmyd3_2m_25_HCHO_2.1m.R1.fastq.gz	PIPELINE_HOME/.test/MOC1_siSmyd3_2m_25_HCHO_2.1m.R2.fastq.gz
```

### 2.2.2 Contrast Manifest (OPTIONAL)
This manifest will include sample information to performed differential comparisons.

An example contrast file:
```
condition1	condition2
MOC1_siSmyd3_2m_25_HCHO	MOC1_siNC_2m_25_HCHO
```