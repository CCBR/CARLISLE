def get_input_bams(wildcards):
    d=dict()
    for dupstatus in DUPSTATUS:
        t="treatment_bam"
        c="control_bam"
        d[t] = join(RESULTSDIR,"bam",wildcards.treatment+"."+dupstatus+".bam")
        d[c] = join(RESULTSDIR,"bam",wildcards.control+"."+dupstatus+".bam")
    return d

def get_input_bedgraphs(wildcards):
    d=dict()
    t="treatment_bedgraph"
    c="control_bedgraph"
    d[t] = join(RESULTSDIR,"bedgraph",wildcards.treatment + "." + wildcards.dupstatus + ".bedgraph")
    d[c] = join(RESULTSDIR,"bedgraph",wildcards.control + "." + wildcards.dupstatus + ".bedgraph")
    return d

rule macs2:
    input:
        fragments_bed = rules.bam2bg.output.fragments_bed,
    output:
        narrowPeak = join(RESULTSDIR,"peaks","{qthresholds}","macs2","{replicate}","{replicate}.{dupstatus}_peaks.narrowPeak"),
        broadPeak = join(RESULTSDIR,"peaks","{qthresholds}","macs2","{replicate}","{replicate}.{dupstatus}_peaks.broadPeak"),
    params:
        replicate = "{replicate}",
        dupstatus = "{dupstatus}",
        qthresholds = "{qthresholds}",
        macs2_genome = config["reference"][GENOME]["macs2_g"],
        outdir = join(RESULTSDIR,"peaks","{qthresholds}","macs2","{replicate}")
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
        macs2 callpeak -t {input.fragments_bed} -f BED -g {params.macs2_genome} --keep-dup all -q {params.qthresholds} -n {params.replicate}.{params.dupstatus} --SPMR --shift 0 --call-summits --nomodel
        macs2 callpeak -t {input.fragments_bed} -f BED -g {params.macs2_genome} --keep-dup all -q {params.qthresholds} -n {params.replicate}.{params.dupstatus} --SPMR --shift 0 --broad --nomodel
        """

localrules: peak2bb_macs2

rule peak2bb_macs2:
    input:
        narrowPeak = join(RESULTSDIR,"peaks","{qthresholds}","macs2","{replicate}","{replicate}.{dupstatus}_peaks.narrowPeak"),
        broadPeak = join(RESULTSDIR,"peaks","{qthresholds}","macs2","{replicate}","{replicate}.{dupstatus}_peaks.broadPeak"),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
    output:
        narrowbb = join(RESULTSDIR,"peaks","{qthresholds}","macs2","{replicate}","{replicate}.{dupstatus}_peaks.narrow.bigbed"),
        broadbb = join(RESULTSDIR,"peaks","{qthresholds}","macs2","{replicate}","{replicate}.{dupstatus}_peaks.broad.bigbed"),
    params:
        replicate="{replicate}",
        dupstatus="{dupstatus}",
        memG = getmemG("peapeak2bb_macs2k2bb"),
    threads: getthreads("peak2bb_macs2")
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

rule seacr:
    input:
        unpack(get_input_bedgraphs)
    output:
        normStringentBed = join(RESULTSDIR,"peaks","{qthresholds}","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.norm.stringent.bed"),
        normRelaxedBed = join(RESULTSDIR,"peaks","{qthresholds}","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.norm.relaxed.bed"),
        nonStringentBed = join(RESULTSDIR,"peaks","{qthresholds}","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.non.stringent.bed"),
        nonRelaxedBed = join(RESULTSDIR,"peaks","{qthresholds}","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.non.relaxed.bed"),
    params:
        treatment = "{treatment}",
        control = "{control}",
        qthresholds="{qthresholds}",
        outdir = join(RESULTSDIR,"peaks","{qthresholds}","seacr","{treatment}_vs_{control}"),
        dupstatus= "{dupstatus}",
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

        SEACR.sh --bedgraph {input.treatment_bedgraph} \\
            --control {input.control_bedgraph} \\
            --normalize norm \\
            --mode stringent \\
            --threshold {params.qthresholds} \\
            --output {params.treatment}_vs_{params.control}.{params.dupstatus}.norm

        SEACR.sh --bedgraph {input.treatment_bedgraph} \\
            --control {input.control_bedgraph} \\
            --normalize norm \\
            --mode relaxed \\
            --threshold {params.qthresholds} \\
            --output {params.treatment}_vs_{params.control}.{params.dupstatus}.norm

        SEACR.sh --bedgraph {input.treatment_bedgraph} \\
            --control {input.control_bedgraph} \\
            --normalize non \\
            --mode stringent \\
            --threshold {params.qthresholds} \\
            --output {params.treatment}_vs_{params.control}.{params.dupstatus}.non

        SEACR.sh --bedgraph {input.treatment_bedgraph} \\
            --control {input.control_bedgraph} \\
            --normalize non \\
            --mode relaxed \\
            --threshold {params.qthresholds} \\
            --output {params.treatment}_vs_{params.control}.{params.dupstatus}.non
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
        normStringentBed = join(RESULTSDIR,"peaks","{qthresholds}","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.norm.stringent.bigbed"),
        normRelaxedBed= join(RESULTSDIR,"peaks","{qthresholds}","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.norm.relaxed.bigbed"),
        nonStringentBed = join(RESULTSDIR,"peaks","{qthresholds}","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.non.stringent.bigbed"),
        nonRelaxedBed = join(RESULTSDIR,"peaks","{qthresholds}","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.non.relaxed.bigbed"),
    params:
        outdir = join(RESULTSDIR,"peaks","{qthresholds}","seacr","{treatment}_vs_{control}"),
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
        qthresholds = "{qthresholds}",
        prefix = join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","{treatment}_vs_{control}.{dupstatus}")
    threads:
        getthreads("gopeaks")
    output:
        narrowPeaks=join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","{treatment}_vs_{control}.{dupstatus}.narrowGo_peaks.bed"),
        broadPeaks=join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","{treatment}_vs_{control}.{dupstatus}.broadGo_peaks.bed"),
    shell:
        """
            {params.gopeaks} -b {input.treatment_bam} -c {input.control_bam} -p {params.qthresholds} -o {params.prefix}.narrowGo
            {params.gopeaks} -b {input.treatment_bam} -c {input.control_bam} -p {params.qthresholds} -o {params.prefix}.broadGo --broad
        """

localrules: bed2bb_gopeaks

rule bed2bb_gopeaks:
    input:
        narrowPeaks=rules.gopeaks.output.narrowPeaks,
        broadPeaks=rules.gopeaks.output.broadPeaks,
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
    output:
        narrowPeaks=join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","{treatment}_vs_{control}.{dupstatus}.narrowGo_peaks.bigbed"),
        broadPeaks=join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","{treatment}_vs_{control}.{dupstatus}.broadGo_peaks.bigbed"),
    params:
        outdir = join(RESULTSDIR,"peaks","{qthresholds}","gopeaks"),
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