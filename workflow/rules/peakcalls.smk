def get_input_bedgraphs(wildcards):
    d=dict()
    for dupstatus in DUPSTATUS:
        t="treatment_bedgraph_"+dupstatus
        c="control_bedgraph_"+dupstatus
        d[t] = join(RESULTSDIR,"bedgraph",wildcards.treatment+"."+dupstatus+".bedgraph")
        d[c] = join(RESULTSDIR,"bedgraph",wildcards.control+"."+dupstatus+".bedgraph")
    return d

def get_input_bams(wildcards):
    d=dict()
    for dupstatus in DUPSTATUS:
        t="treatment_bam_"+dupstatus
        c="control_bam_"+dupstatus
        d[t] = join(RESULTSDIR,"bam",wildcards.treatment+"."+dupstatus+".bam")
        d[c] = join(RESULTSDIR,"bam",wildcards.control+"."+dupstatus+".bam")
    return d

rule gopeaks:
    '''
    ./gopeaks -b /data/CCBR/projects/ccbr1155/CS031014/carlisle_220920/results/bam/53_H3K4me3_1.dedup.bam 
        -c /data/CCBR/projects/ccbr1155/CS031014/carlisle_220920/results/bam/igG_1.dedup.bam -o /data/CCBR/projects/ccbr1155/CS031014/gopeaks/53_H3K4me3_1
    '''
    input:
        unpack(get_input_bams)
    params:
        gopeaks=TOOLS["gopeaks"],
        treatment = "{treatment}",
        control = "{control}",
        dupstatus = "{dupstatus}",
        prefix = join(RESULTSDIR,"peaks","gopeaks","{treatment}_vs_{control}.dedup")
    threads:
        getthreads("gopeaks")
    output:
        narrowPeaks=join(RESULTSDIR,"peaks","gopeaks","{treatment}_vs_{control}.{dupstatus}.narrowGo_peaks.bed"),
        broadPeaks=join(RESULTSDIR,"peaks","gopeaks","{treatment}_vs_{control}.{dupstatus}.broadGo_peaks.bed"),
    shell:
        """
            if [[ {params.dupstatus} == "dedup" ]]; then
                {params.gopeaks} -b {input.treatment_bam_dedup} -c {input.control_bam_dedup} -o {params.prefix}.narrowGo
                {params.gopeaks} -b {input.treatment_bam_dedup} -c {input.control_bam_dedup} -o {params.prefix}.broadGo --broad
            else
                {params.gopeaks} -b {input.treatment_bam_dedup} -c {input.control_bam_no_dedup} -o {params.prefix}.narrowGo
                {params.gopeaks} -b {input.treatment_bam_dedup} -c {input.control_bam_no_dedup} -o {params.prefix}.broadGo --broad
            fi
        """

rule macs2:
    input:
        fragments_bed = rules.bam2bg.output.fragments_bed,
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
        if [[ ! -d {params.outdir} ]];then mkdir -p {params.outdir};fi
        cd {params.outdir}
        macs2 callpeak -t {input.fragments_bed} -f BED -g {params.macs2_genome} --keep-dup all -p 1e-5 -n {params.replicate}.{params.dupstatus} --SPMR --shift 0 --call-summits --nomodel
        macs2 callpeak -t {input.fragments_bed} -f BED -g {params.macs2_genome} --keep-dup all -p 1e-5 -n {params.replicate}.{params.dupstatus} --SPMR --shift 0 --broad --nomodel
        """

rule seacr:
    input:
        unpack(get_input_bedgraphs)
    output:
        normStringentBed = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.norm.stringent.bed"),
        normRelaxedBed = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.norm.relaxed.bed"),
        nonStringentBed = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.non.stringent.bed"),
        nonRelaxedBed = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.non.relaxed.bed"),
    params:
        treatment = "{treatment}",
        control = "{control}",
        outdir = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}")
    threads: getthreads("seacr")
    envmodules:
        TOOLS["seacr"],
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

localrules: peak2bb

rule peak2bb:
    input:
        narrowPeak = join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.narrowPeak"),
        broadPeak = join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.broadPeak"),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
    output:
        narrowbb = join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.narrow.bigbed"),
        broadbb = join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.broad.bigbed"),
    params:
        replicate="{replicate}",
        dupstatus="{dupstatus}",
        memG = getmemG("peak2bb"),
    threads: getthreads("peak2bb")
    envmodules:
        TOOLS["ucsc"]
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
        cut -f1-3 {input.narrowPeak} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.narrow.bed
        bedToBigBed -type=bed3 ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.narrow.bed {input.genome_len} {output.narrowbb}
        cut -f1-3 {input.broadPeak} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.broad.bed
        bedToBigBed -type=bed3 ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.broad.bed {input.genome_len} {output.broadbb}
        """ 

localrules: bed2bb_seacr

rule bed2bb_seacr:
    input:
        normStringentBed = rules.seacr.output.normStringentBed,
        normRelaxedBd = rules.seacr.output.normRelaxedBed,
        nonStringentBed = rules.seacr.output.nonStringentBed,
        nonRelaxedBed = rules.seacr.output.nonRelaxedBed,
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
    output:
        normStringentBed = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.norm.stringent.bigbed"),
        normRelaxedBed= join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.norm.relaxed.bigbed"),
        nonStringentBed = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.non.stringent.bigbed"),
        nonRelaxedBed = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.non.relaxed.bigbed"),
    params:
        outdir = join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}"),
        memG = getmemG("bed2bb_seacr")
    threads: getthreads("bed2bb_seacr")
    envmodules:
        TOOLS["ucsc"]
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
        for fullpath in {input};do
            filename=$(basename -- "$fullpath") 
            extension="${{filename##*.}}" 
            bn="${{filename%.*}}"
            if [[ "$extension" == "bed" ]];then
                cut -f1-3 $fullpath | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/${{bn}}.bed
                bedToBigBed -type=bed3 ${{TMPDIR}}/${{bn}}.bed {input.genome_len} {params.outdir}/${{bn}}.bigbed    
            fi
        done
        """

localrules: bed2bb_gopeaks

rule bed2bb_gopeaks:
    input:
        narrowPeaks=rules.gopeaks.output.narrowPeaks,
        broadPeaks=rules.gopeaks.output.broadPeaks,
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
    output:
        narrowPeaks=join(RESULTSDIR,"peaks","gopeaks","{treatment}_vs_{control}.{dupstatus}.narrowGo_peaks.bigbed"),
        broadPeaks=join(RESULTSDIR,"peaks","gopeaks","{treatment}_vs_{control}.{dupstatus}.broadGo_peaks.bigbed"),
    params:
        outdir = join(RESULTSDIR,"peaks","gopeaks"),
        dupstatus = "{dupstatus}",
        memG = getmemG("bed2bb_gopeaks")
    threads: 
        getthreads("bed2bb_gopeaks")
    envmodules:
        TOOLS["ucsc"]
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
        for fullpath in {input};do
            filename=$(basename -- "$fullpath") 
            extension="${{filename##*.}}" 
            bn="${{filename%.*}}"
            if [[ "$extension" == "bed" ]];then
                cut -f1-3 $fullpath | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/${{bn}}.bed
                bedToBigBed -type=bed3 ${{TMPDIR}}/${{bn}}.bed {input.genome_len} {params.outdir}/${{bn}}.bigbed    
            fi
        done
        """
