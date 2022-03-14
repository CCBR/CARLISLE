def get_input_bedgraphs(wildcards):
    d=dict()
    for dupstatus in DUPSTATUS:
        t="treatment_bedgraph_"+dupstatus
        c="control_bedgraph_"+dupstatus
        d[t] = join(RESULTSDIR,"bedgraph",wildcards.treatment+"."+dupstatus+".bedgraph")
        d[c] = join(RESULTSDIR,"bedgraph",wildcards.control+"."+dupstatus+".bedgraph")
    return d

rule macs2:
    input:
        fragments_bed = join(RESULTSDIR,"tmp","fragments","{replicate}.{dupstatus}.fragments.bed")
    output:
        narrowPeak = join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.narrowPeak"),
        broadPeak = join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.broadPeak"),
    params:
        replicate = "{replicate}",
        dupstatus = "{dupstatus}",
        macs2_genome = config["reference"][GENOME]["macs2_g"],
        outdir = join(RESULTSDIR,"peaks","macs2","{replicate}")
    threads: getthreads("macs2")
    envmodules:
        TOOLS["macs2"]
    shell:"""
set -exo pipefail
if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then 
    TMPDIR="/lscratch/$SLURM_JOB_ID"
else
    dirname=$(basename $(mktemp))
    TMPDIR="/dev/shm/$dirname"
    mkdir -p $TMPDIR
fi
if [[ ! -d {params.outdir} ]];then mkdir -p {params.outdir};fi
cd {params.outdir}
macs2 callpeak -t {input.fragments_bed} -f BED -g {params.macs2_genome} --keep-dup all -p 1e-5 -n {params.replicate}.{params.dupstatus} --SPMR --shift 0 --call-summits --nomodel
macs2 callpeak -t {input.fragments_bed} -f BED -g {params.macs2_genome} --keep-dup all -p 1e-5 -n {params.replicate}.{params.dupstatus} --SPMR --shift 0 --broad --nomodel
"""

rule seacr:
    input:
        unpack(get_input_bedgraphs)
    output:
        normStringentBed_dedup = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.dedup.norm.stringent.bed"),
        normRelaxedBed_dedup = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.dedup.norm.relaxed.bed"),
        nonStringentBed_dedup = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.dedup.non.stringent.bed"),
        nonRelaxedBed_dedup = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.dedup.non.relaxed.bed"),
        normStringentBed_no_dedup = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.no_dedup.norm.stringent.bed"),
        normRelaxedBed_no_dedup = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.no_dedup.norm.relaxed.bed"),
        nonStringentBed_no_dedup = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.no_dedup.non.stringent.bed"),
        nonRelaxedBed_no_dedup = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.no_dedup.non.relaxed.bed"),
    params:
        treatment = "{treatment}",
        control = "{control}",
        outdir = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}")
    threads: getthreads("seacr")
    envmodules:
        TOOLS["seacr"],
    shell:"""
set -exo pipefail
if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then 
    TMPDIR="/lscratch/$SLURM_JOB_ID"
else
    dirname=$(basename $(mktemp))
    TMPDIR="/dev/shm/$dirname"
    mkdir -p $TMPDIR
fi
cd {params.outdir}

SEACR.sh --bedgraph {input.treatment_bedgraph_dedup} \\
    --control {input.control_bedgraph_dedup} \\
    --normalize norm \\
    --mode stringent \\
    --output {params.treatment}_vs_{params.control}.dedup.norm

SEACR.sh --bedgraph {input.treatment_bedgraph_dedup} \\
    --control {input.control_bedgraph_dedup} \\
    --normalize norm \\
    --mode relaxed \\
    --output {params.treatment}_vs_{params.control}.dedup.norm

SEACR.sh --bedgraph {input.treatment_bedgraph_dedup} \\
    --control {input.control_bedgraph_dedup} \\
    --normalize non \\
    --mode stringent \\
    --output {params.treatment}_vs_{params.control}.dedup.non

SEACR.sh --bedgraph {input.treatment_bedgraph_dedup} \\
    --control {input.control_bedgraph_dedup} \\
    --normalize non \\
    --mode relaxed \\
    --output {params.treatment}_vs_{params.control}.dedup.non

SEACR.sh --bedgraph {input.treatment_bedgraph_no_dedup} \\
    --control {input.control_bedgraph_no_dedup} \\
    --normalize norm \\
    --mode stringent \\
    --output {params.treatment}_vs_{params.control}.no_dedup.norm

SEACR.sh --bedgraph {input.treatment_bedgraph_no_dedup} \\
    --control {input.control_bedgraph_no_dedup} \\
    --normalize norm \\
    --mode relaxed \\
    --output {params.treatment}_vs_{params.control}.no_dedup.norm

SEACR.sh --bedgraph {input.treatment_bedgraph_no_dedup} \\
    --control {input.control_bedgraph_no_dedup} \\
    --normalize non \\
    --mode stringent \\
    --output {params.treatment}_vs_{params.control}.no_dedup.non

SEACR.sh --bedgraph {input.treatment_bedgraph_no_dedup} \\
    --control {input.control_bedgraph_no_dedup} \\
    --normalize non \\
    --mode relaxed \\
    --output {params.treatment}_vs_{params.control}.no_dedup.non

"""

