def get_input_bedgraphs(wildcards):
    d=dict()
    d['treatment_bedgraph'] = join(RESULTSDIR,"bedgraph",wildcards.treatment+".bedgraph")
    d['control_bedgraph'] = join(RESULTSDIR,"bedgraph",wildcards.control+".bedgraph")
    return d

rule seacr:
    input:
        unpack(get_input_bedgraphs)
    output:
        normStringentBed = join(RESULTSDIR,"peaks","{treatment}_vs_{control}","{treatment}_vs_{control}.norm.stringent.bed"),
        normRelaxedBed = join(RESULTSDIR,"peaks","{treatment}_vs_{control}","{treatment}_vs_{control}.norm.relaxed.bed"),
        nonStringentBed = join(RESULTSDIR,"peaks","{treatment}_vs_{control}","{treatment}_vs_{control}.non.stringent.bed"),
        nonRelaxedBed = join(RESULTSDIR,"peaks","{treatment}_vs_{control}","{treatment}_vs_{control}.non.relaxed.bed"),
    params:
        treatment = "{treatment}",
        control = "{control}",
        outdir = join(RESULTSDIR,"peaks","{treatment}_vs_{control}")
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

SEACR.sh --bedgraph {input.treatment_bedgraph} \\
    --control {input.control_bedgraph} \\
    --normalize norm \\
    --mode stringent \\
    --output {params.treatment}_vs_{params.control}.norm

SEACR.sh --bedgraph {input.treatment_bedgraph} \\
    --control {input.control_bedgraph} \\
    --normalize norm \\
    --mode relaxed \\
    --output {params.treatment}_vs_{params.control}.norm

SEACR.sh --bedgraph {input.treatment_bedgraph} \\
    --control {input.control_bedgraph} \\
    --normalize non \\
    --mode stringent \\
    --output {params.treatment}_vs_{params.control}.non

SEACR.sh --bedgraph {input.treatment_bedgraph} \\
    --control {input.control_bedgraph} \\
    --normalize non \\
    --mode relaxed \\
    --output {params.treatment}_vs_{params.control}.non

"""

