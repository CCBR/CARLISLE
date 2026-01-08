#########################################################
# IMPORT PYTHON LIBRARIES HERE
#########################################################
import sys
import os
import pandas as pd
import yaml
import pprint
import shutil
import uuid
import re
pp = pprint.PrettyPrinter(indent=4)
#########################################################

#########################################################
# FILE-ACTION FUNCTIONS
#########################################################
def check_existence(filename):
  if not os.path.exists(filename):
    exit("# File: %s does not exists!"%(filename))
  return True

def check_readaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.R_OK):
    exit("# File: %s exists, but cannot be read!"%(filename))
  return True

def check_writeaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.W_OK):
    exit("# File: %s exists, but cannot be read!"%(filename))
  return True

def get_file_size(filename):
    filename=filename.strip()
    if check_readaccess(filename):
        return os.stat(filename).st_size
#########################################################

#########################################################
# DEFINE CONFIG FILE AND READ IT
#########################################################
CONFIGFILE = str(workflow.overwrite_configfiles[0])

# set memory limit
# used for sambamba sort, etc
MEMORYG="100G"

# read in various dirs from config file
WORKDIR=config['workdir']
RESULTSDIR=join(WORKDIR,"results")

PIPELINE_HOME=config['pipeline_home']
check_existence(PIPELINE_HOME)
RESOURCESDIR=join(PIPELINE_HOME,"resources")
check_existence(RESOURCESDIR)

RANDOMSTR=uuid.uuid4()

# get scripts folder
try:
    SCRIPTSDIR = config["scriptsdir"]
except KeyError:
    SCRIPTSDIR = join(WORKDIR,"scripts")
check_existence(SCRIPTSDIR)

if not os.path.exists(join(WORKDIR,"fastqs")):
    os.mkdir(join(WORKDIR,"fastqs"))
if not os.path.exists(join(RESULTSDIR)):
    os.mkdir(join(RESULTSDIR))
for f in ["samplemanifest"]:
    check_readaccess(config[f])

#########################################################
# CREATE SAMPLE DATAFRAME
#########################################################
# each line in the samplemanifest is a replicate
# multiple replicates belong to a sample
# currently 1-4 replicates per sample are supported
print("#"*100)
print("# Checking Sample Manifest...")
SAMPLE2REPLICATES = dict()
REPLICATE2SAMPLE = dict()
df=pd.read_csv(config["samplemanifest"],sep="\t",header=0)
df['replicateName']=df.apply(lambda row:row['sampleName']+"_"+str(row['replicateNumber']),axis=1)
REPLICATES = list(df.replicateName.unique())
replicateName2R1 = dict()
replicateName2R2 = dict()
SKIP_FASTQ_CHECKS = bool(config.get("skip_fastq_checks", False)) or bool(getattr(workflow, "dryrun", False))
for r in REPLICATES:
    sn=df[df['replicateName']==r].iloc[0].sampleName
    if not sn in SAMPLE2REPLICATES:
        SAMPLE2REPLICATES[sn]=[]
    SAMPLE2REPLICATES[sn].append(r)
    REPLICATE2SAMPLE[r]=sn
    r1=df[df['replicateName']==r].iloc[0].path_to_R1
    r1new=join(WORKDIR,"fastqs",r+".R1.fastq.gz")
    if SKIP_FASTQ_CHECKS:
        os.makedirs(os.path.dirname(r1new), exist_ok=True)
        if not os.path.exists(r1new):
            open(r1new, 'a').close()
    else:
        check_readaccess(r1)
        if not os.path.exists(r1new):
            os.symlink(r1,r1new)
    replicateName2R1[r]=r1new
    r2=df[df['replicateName']==r].iloc[0].path_to_R2
    r2new=join(WORKDIR,"fastqs",r+".R2.fastq.gz")
    if SKIP_FASTQ_CHECKS:
        os.makedirs(os.path.dirname(r2new), exist_ok=True)
        if not os.path.exists(r2new):
            open(r2new, 'a').close()
    else:
        check_readaccess(r2)
        if not os.path.exists(r2new):
            os.symlink(r2,r2new)
    replicateName2R2[r]=r2new

print("# Samples and their replicates")
print("# SampleName        Replicates")
for k,v in SAMPLE2REPLICATES.items():
    print("# "+k+"         "+",".join(v))
print("# Read access to all fastq files in confirmed!")
print("# Sample manifest is confirmed!")

print("#"*100)
print("# Checking Contrast Manifest...:")
process_replicates = []
TREATMENTS = []
CONTROLS = []
TREATMENT_CONTROL_LIST=[]
TREATMENT_WITHOUTCONTROL_LIST=[]
TREAT_to_CONTRL_DICT=dict()
for i,t in enumerate(list(df[df['isControl']=="N"]['replicateName'].unique())):
    crow=df[df['replicateName']==t].iloc[0]
    c=crow.controlName+"_"+str(int(crow.controlReplicateNumber))
    if not c in REPLICATES:
        print("# Control NOT found for sampleName_replicateNumber:"+t)
        print("# "+config["samplemanifest"]+" has no entry for sample:"+crow.controlName+"  replicateNumber:"+str(crow.controlReplicateNumber))
        exit()
    print("## "+str(i+1)+") "+t+"        "+c)
    process_replicates.extend([t,c])
    TREATMENTS.append(t)
    CONTROLS.append(c)
    TREATMENT_CONTROL_LIST.append(t+"_vs_"+c)
    TREATMENT_WITHOUTCONTROL_LIST.append(t+"_vs_nocontrol")
    TREAT_to_CONTRL_DICT[t]=c
process_replicates=list(set(process_replicates))
if len(process_replicates)!=len(REPLICATES):
    not_to_process = set(REPLICATES) - set(process_replicates)
    print("# Following replicates will not be processed as they are not part of any Treatment Control combination!")
    for i in not_to_process:
        print("# "+i)
    REPLICATES = process_replicates
print("# Contrast manifest is confirmed!")

# write out the treatment:cntrl
fpath=join(RESULTSDIR,"treatment_control_list.txt")
originalDF = pd.DataFrame(TREAT_to_CONTRL_DICT.items()).rename(columns={0: 'key', 1: 'val'})
split_keysDF = pd.DataFrame(originalDF['key'].str.split(':').tolist())
finalDF = split_keysDF.join(originalDF['val'])
finalDF.to_csv(fpath, sep='\t', header=False, index=False)

# validate pool_controls setting
if config.get("pool_controls", False):
    if not CONTROLS or len(CONTROLS) == 0:
        print("#"*100)
        print("# ERROR: pool_controls is set to 'true' but no control samples are present!")
        print("# When pool_controls is enabled, peaks are called in both 'individual' and 'pooled' modes.")
        print("# Without controls, 'individual' mode would produce peaks without any control,")
        print("# which defeats the purpose of having both modes.")
        print("#")
        print("# Please either:")
        print("#   1. Set pool_controls: false in config.yaml, OR")
        print("#   2. Add control samples to your samples.tsv")
        print("#"*100)
        exit(1)
    # Check if controls have replicates for pooling
    control_samples = list(set([c.rsplit('_',1)[0] for c in CONTROLS]))
    if len(control_samples) == len(CONTROLS):
        print("#"*100)
        print("# WARNING: pool_controls is set to 'true' but controls have no replicates to pool!")
        print("# Each control sample has only 1 replicate:")
        for ctrl in CONTROLS:
            print(f"#   - {ctrl}")
        print("# Pooling will have no effect. Consider setting pool_controls: false")
        print("#"*100)

# Get unique control sample names from samples.tsv where isControl == Y
# This is the authoritative list of base control names (e.g., PBS_IgG, mSTAR_IgG)
# Works with any naming convention including multiple underscores
CONTROL_SAMPLES = list(df[df['isControl']=="Y"]['sampleName'].unique())

# Set control modes based on pool_controls setting
# When pool_controls is true, we run analysis with both individual and pooled controls
# When pool_controls is false, we only run with individual controls
CONTROL_MODES = ["individual", "pooled"] if config.get("pool_controls", False) else ["individual"]

# Create pooled treatment control lists where control replicate numbers are removed
# For example: mSTAR_H3K27me3_1_vs_mSTAR_IgG_1 becomes mSTAR_H3K27me3_1_vs_mSTAR_IgG
TREATMENT_CONTROL_LIST_POOLED = []
for tc in TREATMENT_CONTROL_LIST:
    parts = tc.split("_vs_")
    if len(parts) == 2:
        treatment = parts[0]
        control = parts[1]
        # Remove replicate number from control (last _N pattern)
        import re
        control_no_rep = re.sub(r'_\d+$', '', control)
        TREATMENT_CONTROL_LIST_POOLED.append(f"{treatment}_vs_{control_no_rep}")
    else:
        TREATMENT_CONTROL_LIST_POOLED.append(tc)

# Create regex patterns for wildcard constraints
# Individual mode: both treatment and control end with _\d+
# Pooled mode: control does NOT end with _\d+ (replicate number removed)
INDIVIDUAL_TREATMENT_CONTROL_PATTERN = "|".join([f"^{re.escape(tc)}$" for tc in TREATMENT_CONTROL_LIST])
POOLED_TREATMENT_CONTROL_PATTERN = "|".join([f"^{re.escape(tc)}$" for tc in TREATMENT_CONTROL_LIST_POOLED])

# Function to validate control_mode and treatment_control_list combinations
def is_valid_control_treatment_combo(control_mode, treatment_control_list):
    """
    Validates that treatment_control_list matches the appropriate pattern for control_mode.
    Individual mode: treatment_control_list must be from TREATMENT_CONTROL_LIST (with replicate numbers on both)
    Pooled mode: treatment_control_list must be from TREATMENT_CONTROL_LIST_POOLED (no replicate on control)
    """
    if control_mode == "individual":
        return treatment_control_list in TREATMENT_CONTROL_LIST
    elif control_mode == "pooled":
        return treatment_control_list in TREATMENT_CONTROL_LIST_POOLED
    return False

# set treatment lists depending on whether controls were used for all peak callers
# macs2 allows for with or without controls; all other callers require with controls
if (config["macs2_control"] == "Y"):
    TREATMENT_LIST_M=TREATMENT_CONTROL_LIST
else:
    TREATMENT_LIST_M=TREATMENT_WITHOUTCONTROL_LIST
TREATMENT_LIST_SG=TREATMENT_CONTROL_LIST

# create duplication and peaktype list
DUPSTATUS=config["dupstatus"]
PEAKTYPE=config["peaktype"]
DUPSTATUS=list(map(lambda x:x.strip(),DUPSTATUS.split(",")))
PEAKTYPE=list(map(lambda x:x.strip(),PEAKTYPE.split(",")))

#separate peak types into caller lists
macs_set=[]
g_set=[]
s_set=[]
for pt in PEAKTYPE:
    if pt == "macs2_narrow" or pt == "macs2_broad":
        macs_set.append(pt)
    elif pt == "gopeaks_narrow" or pt == "gopeaks_broad":
        g_set.append(pt)
    elif pt == "seacr_stringent" or pt == "seacr_relaxed":
        s_set.append(pt)
    else:
        print("A peak type combination was used that is non-compatiable")
        print(pt)
        exit()
PEAKTYPE_M=list(set(macs_set))
PEAKTYPE_G=list(set(g_set))
PEAKTYPE_S=list(set(s_set))

# set threshold settings
FDRCUTOFF = config["contrasts_fdr_cutoff"]
LFCCUTOFF = config["contrasts_lfc_cutoff"]
QTRESHOLDS=config["quality_thresholds"]
QTRESHOLDS=list(map(lambda x:x.strip(),QTRESHOLDS.split(",")))

# set contrast settings
if config["run_contrasts"]:
    print("#"*100)
    print("# Checking contrasts to run...")
    contrasts_table = config["contrasts"]
    check_readaccess(contrasts_table)
    contrasts_df=pd.read_csv(contrasts_table,sep="\t",header=0)
    # contrasts_df = contrasts_df.reset_index()  # make sure indexes pair with number of rows
    # print(contrasts_df)
    CONTRASTS=dict()
    CONTRAST_EXCLUDE=dict()
    C1s=[]
    C2s=[]
    DS=[]
    QS=[]
    PT=[]
    CONTRAST_LIST=[]
    SAMPLESINCONTRAST=list()
    for index, row in contrasts_df.iterrows():
        c1 = row['condition1']
        c2 = row['condition2']
        if not c1 in SAMPLE2REPLICATES:
            print(" # %s condition1 in %s has no samples/replicates!"%(c1,contrasts_table))
            exit()
        if not c2 in SAMPLE2REPLICATES:
            print(" # %s condition2 in %s has no samples/replicates!"%(c2,contrasts_table))
            exit()
        # Ensure both conditions have at least one replicate listed
        if len(SAMPLE2REPLICATES[c1]) == 0:
            print(" # %s has no replicates in samples manifest!" % (c1))
            exit()
        if len(SAMPLE2REPLICATES[c2]) == 0:
            print(" # %s has no replicates in samples manifest!" % (c2))
            exit()
        for ds in DUPSTATUS:
            for pt in PEAKTYPE:
                for qt in QTRESHOLDS:
                    QS.append(qt)
                    C1s.append(c1)
                    C2s.append(c2)
                    DS.append(ds)
                    PT.append(pt)
                    contrast_name=c1+"_vs_"+c2+"__"+ds+"__"+pt
                    CONTRASTS[contrast_name]=[c1,c2,ds,pt]
                    CONTRAST_LIST.append(c1+"_vs_"+c2)
        # Optional: outlier/exclude replicates column (comma-separated replicate names)
        if 'exclude_replicates' in contrasts_df.columns:
            excl_val = str(row['exclude_replicates']).strip()
            if (excl_val.lower() != 'nan') and (excl_val != '') and (excl_val.lower() != 'none'):
                excludes = list(map(lambda x: x.strip(), excl_val.split(",")))
                # Validate each excluded replicate
                for ex in excludes:
                    if ex == "":
                        continue
                    m = re.match(r'^(.+?)_(\d+)$', ex)
                    if not m:
                        print(" # Invalid replicate format in exclude_replicates: %s. Expected sampleName_N" % (ex))
                        exit()
                    base_sample = m.group(1)
                    if base_sample not in [c1, c2]:
                        print(" # Excluded replicate %s does not belong to %s or %s" % (ex, c1, c2))
                        exit()
                    if ex not in REPLICATES:
                        print(" # Excluded replicate %s not found in samples manifest replicates" % (ex))
                        exit()
                CONTRAST_EXCLUDE[c1+"_vs_"+c2] = excludes
            else:
                CONTRAST_EXCLUDE[c1+"_vs_"+c2] = []
        else:
            CONTRAST_EXCLUDE[c1+"_vs_"+c2] = []
        SAMPLESINCONTRAST.append(c1)
        SAMPLESINCONTRAST.append(c2)
    SAMPLESINCONTRAST=list(set(SAMPLESINCONTRAST))
    print("# Will run the following contrasts:")
    print("# CONDITION1\tCONDITION2\tDUPSTATUS\tPEAKTYPE")
    i=0
    for k,v in CONTRASTS.items():
        i+=1
        print("# %d) %s\t%s\t%s\t%s"%(i,v[0],v[1],v[2],v[3]))
    # Log any replicate exclusions detected per contrast
    for cl, excl in CONTRAST_EXCLUDE.items():
        if excl:
            print("# Excluding replicates for %s: %s" % (cl, ", ".join(excl)))
    rg_file = join(RESULTSDIR,"replicate_sample.tsv")
    replicates_in_file = []

    if os.path.exists(rg_file) and get_file_size(rg_file) != 0:
        rg_df = pd.read_csv(rg_file,sep="\t",header=None)
        rg_df.columns = ["replicateName","sampleName"]
        replicates_in_file = list(rg_df["replicateName"].unique())

    replicates = []
    for sample in SAMPLESINCONTRAST:
        replicates.extend(SAMPLE2REPLICATES[sample])

    if set(replicates) & set(replicates_in_file) != set(replicates):
        # rg_file needs to be recreated to match the new setting from contrasts.tab .. may you added a new contrast
        if os.path.exists(rg_file): os.remove(rg_file)
        # once file is removed it will be recreated by rule create_replicate_sample_table

# create CONTROL : TREATMENT dictionary
CONTROL_to_TREAT_DICT={}
for t in TREAT_to_CONTRL_DICT:
    c=TREAT_to_CONTRL_DICT[t]
    if c not in CONTROL_to_TREAT_DICT:
        CONTROL_to_TREAT_DICT[c]=[]
    CONTROL_to_TREAT_DICT[c].append(t)

#########################################################
# READ IN TOOLS REQUIRED BY PIPELINE
# THESE INCLUDE LIST OF BIOWULF MODULES (AND THEIR VERSIONS)
# MAY BE EMPTY IF ALL TOOLS ARE DOCKERIZED
#########################################################
## Load tools from YAML file
try:
    TOOLSYAML = config["tools"]
except KeyError:
    TOOLSYAML = join(WORKDIR,"config","tools.yaml")
check_readaccess(TOOLSYAML)
with open(TOOLSYAML) as f:
    TOOLS = yaml.safe_load(f)
#########################################################

#########################################################
# READ CLUSTER PER-RULE REQUIREMENTS
#########################################################
## Load cluster.json
try:
    CLUSTERYAML = config["CLUSTERYAML"]
except KeyError:
    CLUSTERYAML = join(WORKDIR,"config","cluster.yaml")
check_readaccess(CLUSTERYAML)
with open(CLUSTERYAML) as json_file:
    CLUSTER = yaml.safe_load(json_file)

## Create lambda functions to allow a way to insert read-in values
## as rule directives
getthreads=lambda rname:int(CLUSTER[rname]["threads"]) if rname in CLUSTER and "threads" in CLUSTER[rname] else int(CLUSTER["__default__"]["threads"])
getmemg=lambda rname:CLUSTER[rname]["mem"] if rname in CLUSTER and "mem" in CLUSTER[rname] else CLUSTER["__default__"]["mem"]
getmemG=lambda rname:getmemg(rname).replace("g","G")

# Function to get input fastq files for a replicate
def get_input_fastqs(wildcards):
    d = dict()
    d["R1"] = replicateName2R1[wildcards.replicate]
    d["R2"] = replicateName2R2[wildcards.replicate]
    return d
#########################################################

#########################################################
# SET OTHER PIPELINE GLOBAL VARIABLES
#########################################################
print("#"*100)
print("# Pipeline Parameters:")
print("# Working dir :",WORKDIR)
print("# Results dir :",RESULTSDIR)
print("# Scripts dir :",SCRIPTSDIR)
print("# Config YAML :",CONFIGFILE)
print("# Sample Manifest :",config["samplemanifest"])
print("# Cluster YAML :",CLUSTERYAML)

GENOME = config["genome"]
print("# Reference genome : ",GENOME)
GENOMEFA = config["reference"][GENOME]["fa"]
check_readaccess(GENOMEFA)
REGIONS = config["reference"][GENOME]["regions"]

GENOMEBLACKLIST = config["reference"][GENOME]["blacklist"]
check_readaccess(GENOMEBLACKLIST)

NORM_METHOD = config["norm_method"].upper()
print("# Norm method : ",NORM_METHOD)
if NORM_METHOD == "SPIKEIN":
    spikein_genome = config["spikein_genome"]
    SPIKED_GENOMEFA = config["spikein_reference"][spikein_genome]["fa"]
    check_readaccess(SPIKED_GENOMEFA)
    print("# Spike-in genome : ",spikein_genome)
elif NORM_METHOD == "LIBRARY":
    SPIKED_GENOMEFA="LIBRARY"
elif NORM_METHOD == "NONE":
    SPIKED_GENOMEFA=""
else:
    print("User must select from one of the three available norm methods: spikein,library, none")
    exit()

BOWTIE2_INDEX = join(WORKDIR,"bowtie2_index")
refdata = dict()
refdata["genome"] = GENOME
refdata["genomefa"] = GENOMEFA
refdata["blacklistbed"] = GENOMEBLACKLIST
refdata["norm_method"] = NORM_METHOD
refdata["spikein_genome"] = SPIKED_GENOMEFA

# set annotation params
S_DISTANCE=config["stitch_distance"]
GENESET_ID=config["geneset_id"]

#########################################################
# CHECK ACCESS TO OTHER RESOURCES
#########################################################
if GENOME == "hg38" or GENOME == "hg19" or GENOME == "hs1":
    check_readaccess(config["reference"][GENOME]["tss_bed"])
    check_readaccess(config["reference"][GENOME]["rose"])

#########################################################
# DEFINE LOCAL RULES
#########################################################
localrules: create_replicate_sample_table

rule create_replicate_sample_table:
    input:
    output:
        rg_file = join(RESULTSDIR,"replicate_sample.tsv")
    run:
        rg_file = output.rg_file

        replicates = []
        for sample in SAMPLESINCONTRAST:
            replicates.extend(SAMPLE2REPLICATES[sample])
        rg = open(rg_file,'w')
        for r in replicates:
            rg.write("%s\t%s\n"%(r,REPLICATE2SAMPLE[r]))
        rg.close()
        print("# %s file created!"%(rg_file))

rule create_reference:
    input:
        genomefa_source=GENOMEFA,
        blacklist_source=GENOMEBLACKLIST,
    output:
        genomefa = join(BOWTIE2_INDEX,"genome.fa"),
        blacklist = join(BOWTIE2_INDEX,"genome.blacklist.bed"),
        bt2 = join(BOWTIE2_INDEX,"ref.1.bt2"),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        ref_len = join(BOWTIE2_INDEX,"ref.len"),
        spikein_len = join(BOWTIE2_INDEX,"spikein.len") ,
        refjson = join(BOWTIE2_INDEX,"ref.yaml")
    params:
        bt2_base=join(BOWTIE2_INDEX,"ref"),
        bowtie2_dir=BOWTIE2_INDEX,
        use_spikein=NORM_METHOD,
        spiked_source=SPIKED_GENOMEFA,
        spiked_output=join(BOWTIE2_INDEX,"spikein.fa"),
        refdata=refdata,
    envmodules:
        TOOLS["bowtie2"],
        TOOLS["samtools"],
        TOOLS["bedtools"],
    threads: getthreads("create_reference")
    shell:
        """
        set -exo pipefail
        if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then
            TMPDIR="/lscratch/$SLURM_JOB_ID"
        else
            dirname=$(basename $(mktemp))
            TMPDIR="/dev/shm/$dirname"
            mkdir -p $TMPDIR
        fi

        # create dir and links
        if [[ -d {params.bowtie2_dir} ]]; then rm -r {params.bowtie2_dir}; fi
        mkdir -p {params.bowtie2_dir}/ref
        ln -s {input.genomefa_source} {output.genomefa}
        ln -s {input.blacklist_source} {output.blacklist}
        if [[ {params.use_spikein} == "SPIKEIN" ]]; then ln -s {params.spiked_source} {params.spiked_output}; fi

        # create json file and store in "tmp" until the reference is built
        mkdir -p {params.bowtie2_dir}/tmp
        echo {params.refdata} > {params.bowtie2_dir}/tmp/ref.yaml

        # create reference dependent on whether or not a spikein control should be used
        if [[ "{params.use_spikein}" != "SPIKEIN" ]];then
            echo "creating ref without a spike-in control"

            # create faidx for genome and spike fasta
            samtools faidx {output.genomefa}

            # mask genome fa with genome blacklist
            bedtools maskfasta -fi {output.genomefa} -bed {output.blacklist} -fo ${{TMPDIR}}/masked_genome.fa

            # build bowtie index
            bowtie2-build --threads {threads} ${{TMPDIR}}/masked_genome.fa {params.bt2_base}

            # create len files
            cut -f1,2 {output.genomefa}.fai > {output.ref_len}
            cp {output.ref_len} {output.genome_len}
            echo -e "NA\t0" > {output.spikein_len}
        else
            echo "creating ref with a spike-in control"

            # create faidx for genome and spike fasta
            samtools faidx {output.genomefa}
            samtools faidx {params.spiked_output}

            # mask genome fa with genome blacklist
            bedtools maskfasta -fi {output.genomefa} -bed {output.blacklist} -fo ${{TMPDIR}}/masked_genome.fa

            # build bowtie index
            bowtie2-build --threads {threads} ${{TMPDIR}}/masked_genome.fa,{params.spiked_output} {params.bt2_base}

            # create len files
            cut -f1,2 {output.genomefa}.fai > {output.ref_len}
            cp {output.ref_len} {output.genome_len}
            cut -f1,2 {params.spiked_output}.fai >> {output.ref_len}
            cut -f1,2 {params.spiked_output}.fai > {output.spikein_len}
        fi

        # copy ref.yaml only after successfully finishing ref index building
        if [[ -f {output.bt2} ]];then
            mv $(dirname {output.bt2})/tmp/ref.yaml {output.refjson}
        fi

        """
