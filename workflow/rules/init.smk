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
for r in REPLICATES:
    sn=df[df['replicateName']==r].iloc[0].sampleName
    if not sn in SAMPLE2REPLICATES:
        SAMPLE2REPLICATES[sn]=[]
    SAMPLE2REPLICATES[sn].append(r)
    REPLICATE2SAMPLE[r]=sn
    r1=df[df['replicateName']==r].iloc[0].path_to_R1
    check_readaccess(r1)
    r1new=join(WORKDIR,"fastqs",r+".R1.fastq.gz")
    if not os.path.exists(r1new):
        os.symlink(r1,r1new)
    replicateName2R1[r]=r1new
    r2=df[df['replicateName']==r].iloc[0].path_to_R2
    check_readaccess(r2)
    r2new=join(WORKDIR,"fastqs",r+".R2.fastq.gz")
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
    c=crow.controlName+"_"+str(crow.controlReplicateNumber)
    if not c in REPLICATES:
        print("# Control NOT found for sampleName_replicateNumber:"+t)
        print("# "+config["samplemanifest"]+" has no entry for sample:"+crow.controlName+"  replicateNumber:"+crow.controlReplicateNumber)
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
finalDF.to_csv(fpath, '\t', header=False, index=False)

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
    elif pt == "seacr_norm_stringent" or pt == "seacr_norm_relaxed" or pt == "seacr_non_stringent" or pt == "seacr_non_relaxed":
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
if config["run_contrasts"] == "Y":
    print("#"*100)
    print("# Checking constrasts to run...")
    contrasts_table = config["contrasts"]
    check_readaccess(contrasts_table)
    contrasts_df=pd.read_csv(contrasts_table,sep="\t",header=0)
    # contrasts_df = contrasts_df.reset_index()  # make sure indexes pair with number of rows
    # print(contrasts_df)
    CONTRASTS=dict()
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
        SAMPLESINCONTRAST.append(c1)
        SAMPLESINCONTRAST.append(c2)
    SAMPLESINCONTRAST=list(set(SAMPLESINCONTRAST))
    print("# Will run the following contrasts:")
    print("# CONDITION1\tCONDITION2\tDUPSTATUS\tPEAKTYPE")
    i=0
    for k,v in CONTRASTS.items():
        i+=1
        print("# %d) %s\t%s\t%s\t%s"%(i,v[0],v[1],v[2],v[3]))
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