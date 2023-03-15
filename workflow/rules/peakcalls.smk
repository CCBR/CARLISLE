# TSV file should include 6 columns
# 1)condition     2)sample
# 53_H3K4me3	    53_H3K4me3_1         
# 3)bed
# /results/peaks/gopeaks/53_H3K4me3_1_vs_HN6_IgG_rabbit_negative_control_1.dedup.broadGo_peaks.bed
# 4)bedgraph                                      5)scaling factor
# /results/bedgraph/53_H3K4me3_1.dedup.bedgraph	86.32596685082872928176	
# 6)bed
# /results/fragments/53_H3K4me3_1.dedup.fragments.bed

rule macs2_narrow:
    '''
    MACS2 can be run with and without a control. This featured is controlled in the config.yaml filen
    with the flag macs2_control: "N" or macs2_control: "Y"
    '''
    input:
        fragments_bed = expand(join(RESULTSDIR,"fragments","{replicate}.{dupstatus}.fragments.bed"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
    output:
        peak_file = join(RESULTSDIR,"peaks","macs2","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.narrow.peaks.bed"),
        bg_file = join(RESULTSDIR,"peaks","macs2","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.narrow.peaks.bigbed.gz"),
    params:
        frag_bed_path=join(RESULTSDIR,"fragments"),
        tc_file="{treatment_control_list}",
        dupstatus = "{dupstatus}",
        qthresholds = "{qthresholds}",
        control_flag = config["macs2_control"],
        macs2_genome = config["reference"][GENOME]["macs2_g"],
        memG = getmemG("macs2"),
    threads: getthreads("macs2")
    envmodules:
        TOOLS["macs2"],
        TOOLS["ucsc"],
        TOOLS["samtools"]
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
        cd $TMPDIR

        # pull treatment and control ids
        treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
        control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`
         
        # set frag file
        treat_bed={params.frag_bed_path}/${{treatment}}.{params.dupstatus}.fragments.bed
        cntrl_bed={params.frag_bed_path}/${{control}}.{params.dupstatus}.fragments.bed

        # run with or without control
        if [[ {params.control_flag} == "Y" ]]; then 
            file_base="${{treatment}}_vs_${{control}}.{params.dupstatus}"

            macs2 callpeak \\
            -t ${{treat_bed}} \\
            -c ${{cntrl_bed}} \\
            -f BED \\
            -g {params.macs2_genome} \\
            --keep-dup all \\
            -q {params.qthresholds} \\
            -n ${{file_base}} \\
            --SPMR \\
            --shift 0 \\
            --call-summits \\
            --nomodel
        else
            file_base="${{treatment}}_vs_${{control}}.{params.dupstatus}"
                                
            macs2 callpeak \\
            -t ${{treat_bed}} \\
            -f BED \\
            -g {params.macs2_genome} \\
            --keep-dup all \\
            -q {params.qthresholds} \\
            -n ${{file_base}} \\
            --SPMR \\
            --shift 0 \\
            --call-summits \\
            --nomodel
        fi
    
        # mv output and rename for consistency
        mv $TMPDIR/${{file_base}}_peaks.narrowPeak {output.peak_file}
        
        # create bigbed files, zip
        cut -f1-3 {output.peak_file} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/narrow.bed
        bedToBigBed -type=bed3 ${{TMPDIR}}/narrow.bed {input.genome_len} ${{TMPDIR}}/narrow.bigbed
        bgzip -c ${{TMPDIR}}/narrow.bigbed > {output.bg_file}
        """

rule macs2_broad:
    '''
    MACS2 can be run with and without a control. This featured is controlled in the config.yaml filen
    with the flag macs2_control: "N" or macs2_control: "Y"
    '''
    input:
        fragments_bed = expand(join(RESULTSDIR,"fragments","{replicate}.{dupstatus}.fragments.bed"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        align_stats = rules.gather_alignstats.output,
    output:
        peak_file = join(RESULTSDIR,"peaks","macs2","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.broad.peaks.bed"),
        bg_file = join(RESULTSDIR,"peaks","macs2","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.broad.peaks.bigbed.gz"),
    params:
        frag_bed_path=join(RESULTSDIR,"fragments"),
        tc_file="{treatment_control_list}",
        dupstatus = "{dupstatus}",
        qthresholds = "{qthresholds}",
        control_flag = config["macs2_control"],
        broadtreshold = config["macs2_broad_peak_threshold"],
        macs2_genome = config["reference"][GENOME]["macs2_g"],
        memG = getmemG("macs2"),
    threads: getthreads("macs2")
    envmodules:
        TOOLS["macs2"],
        TOOLS["ucsc"],
        TOOLS["samtools"]
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
        cd $TMPDIR

        # pull treatment and control ids
        treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
        control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`
            
        # set frag file
        treat_bed={params.frag_bed_path}/${{treatment}}.{params.dupstatus}.fragments.bed
        cntrl_bed={params.frag_bed_path}/${{control}}.{params.dupstatus}.fragments.bed

        # run with or without control
        if [[ {params.control_flag} == "Y" ]]; then 
            file_base="${{treatment}}_vs_${{control}}.{params.dupstatus}"
                    
            macs2 callpeak \\
            -t ${{treat_bed}} \\
            -c ${{cntrl_bed}} \\
            -f BED \\
            -g {params.macs2_genome} \\
            --keep-dup all \\
            -q {params.qthresholds} \\
            -n ${{file_base}} \\
            --SPMR \\
            --shift 0 \\
            --broad --broad-cutoff {params.broadtreshold} \\
            --nomodel
        else
            file_base="${{treatment}}_vs_nocontrol.{params.dupstatus}"
                    
            macs2 callpeak \\
            -t ${{treat_bed}} \\
            -f BED \\
            -g {params.macs2_genome} \\
            --keep-dup all \\
            -q {params.qthresholds} \\
            -n ${{file_base}} \\
            --SPMR \\
            --shift 0 \\
            --broad --broad-cutoff {params.broadtreshold} \\
            --nomodel
        fi
            
        # mv output and rename for consistency
        mv $TMPDIR/${{file_base}}_peaks.broadPeak {output.peak_file}

        # create bigbed files
        cut -f1-3 {output.peak_file} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/broad.bed
        bedToBigBed -type=bed3 ${{TMPDIR}}/broad.bed {input.genome_len} ${{TMPDIR}}/broad.bigbed
        bgzip -c ${{TMPDIR}}/broad.bigbed > {output.bg_file}
        """

rule seacr_stringent:
    input:
        tc_list=join(RESULTSDIR,"treatment_control_list.txt"),
        bg = expand(join(RESULTSDIR,"bedgraph","{replicate}.{dupstatus}.bedgraph"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        align_stats = rules.gather_alignstats.output,
    output:
        peak_file_norm = join(RESULTSDIR,"peaks","seacr","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.norm_stringent.peaks.bed"),
        peak_file_non = join(RESULTSDIR,"peaks","seacr","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.non_stringent.peaks.bed"),
        bg_file_norm = join(RESULTSDIR,"peaks","seacr","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.norm_stringent.peaks.bigbed.gz"),
        bg_file_non = join(RESULTSDIR,"peaks","seacr","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.non_stringent.peaks.bigbed.gz"),
    params:
        bg_path=join(RESULTSDIR,"bedgraph"),
        frag_bed_path=join(RESULTSDIR,"fragments"),
        tc_file="{treatment_control_list}",
        qthresholds="{qthresholds}",
        dupstatus= "{dupstatus}",
        memG = getmemG("seacr"),
    threads: getthreads("seacr")
    envmodules:
        TOOLS["seacr"],
        TOOLS["ucsc"],
        TOOLS["samtools"]
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
        cd $TMPDIR

        # pull treatment and control ids
        treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
        control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`
         
        # set frag file
        treat_bed={params.bg_path}/${{treatment}}.{params.dupstatus}.bedgraph
        cntrl_bed={params.bg_path}/${{control}}.{params.dupstatus}.bedgraph

        # run norm
        SEACR.sh --bedgraph ${{treat_bed}} \\
            --control ${{cntrl_bed}} \\
            --normalize norm \\
            --mode stringent \\
            --threshold {params.qthresholds} \\
            --output norm

        # run non
        SEACR.sh --bedgraph ${{treat_bed}} \\
            --control ${{cntrl_bed}} \\
            --normalize non \\
            --mode stringent \\
            --threshold {params.qthresholds} \\
            --output non

        # mv output and rename for consistency
        mv $TMPDIR/norm.stringent.bed {output.peak_file_norm}
        mv $TMPDIR/non.stringent.bed {output.peak_file_non}

        # create bigbed files
        cut -f1-3 {output.peak_file_norm} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/norm_stringent.bed
        bedToBigBed -type=bed3 ${{TMPDIR}}/norm_stringent.bed {input.genome_len} ${{TMPDIR}}/norm_stringent.bigbed

        cut -f1-3 {output.peak_file_non} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/non_stringent.bed
        bedToBigBed -type=bed3 ${{TMPDIR}}/non_stringent.bed {input.genome_len} ${{TMPDIR}}/non_stringent.bigbed

        bgzip -c ${{TMPDIR}}/norm_stringent.bigbed > {output.bg_file_norm}
        bgzip -c ${{TMPDIR}}/non_stringent.bigbed > {output.bg_file_non}
        """

rule seacr_relaxed:
    input:
        bg = expand(join(RESULTSDIR,"bedgraph","{replicate}.{dupstatus}.bedgraph"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        align_stats = rules.gather_alignstats.output,
    output:
        peak_file_norm = join(RESULTSDIR,"peaks","seacr","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.norm_relaxed.peaks.bed"),
        peak_file_non = join(RESULTSDIR,"peaks","seacr","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.non_relaxed.peaks.bed"),
        bg_file_norm = join(RESULTSDIR,"peaks","seacr","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.norm_relaxed.peaks.bigbed.gz"),
        bg_file_non = join(RESULTSDIR,"peaks","seacr","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.non_relaxed.peaks.bigbed.gz"),
    params:
        bg_path=join(RESULTSDIR,"bedgraph"),
        frag_bed_path=join(RESULTSDIR,"fragments"),
        tc_file="{treatment_control_list}",
        qthresholds="{qthresholds}",
        dupstatus= "{dupstatus}",
        memG = getmemG("seacr")
    threads: getthreads("seacr")
    envmodules:
        TOOLS["seacr"],
        TOOLS["ucsc"],
        TOOLS["samtools"]
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
            cd $TMPDIR

            # pull treatment and control ids
            treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
            control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`
            
            # set frag file
            treat_bed={params.bg_path}/${{treatment}}.{params.dupstatus}.bedgraph
            cntrl_bed={params.bg_path}/${{control}}.{params.dupstatus}.bedgraph

            # run norm
            SEACR.sh --bedgraph  ${{treat_bed}} \\
                --control ${{cntrl_bed}} \\
                --normalize norm \\
                --mode relaxed \\
                --threshold {params.qthresholds} \\
                --output norm    

            # run non
            SEACR.sh --bedgraph ${{treat_bed}} \\
                --control ${{cntrl_bed}} \\
                --normalize non \\
                --mode relaxed \\
                --threshold {params.qthresholds} \\
                --output non

            # mv output and rename for consistency
            mv $TMPDIR/norm.relaxed.bed {output.peak_file_norm}
            mv $TMPDIR/non.relaxed.bed {output.peak_file_non}

            # create bigbed files
            cut -f1-3 {output.peak_file_norm} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/norm_relaxed.bed
            bedToBigBed -type=bed3 ${{TMPDIR}}/norm_relaxed.bed {input.genome_len} ${{TMPDIR}}/norm_relaxed.bigbed

            cut -f1-3 {output.peak_file_non} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/non_relaxed.bed
            bedToBigBed -type=bed3 ${{TMPDIR}}/non_relaxed.bed {input.genome_len} ${{TMPDIR}}/non_relaxed.bigbed

            bgzip -c ${{TMPDIR}}/norm_relaxed.bigbed > {output.bg_file_norm}
            bgzip -c ${{TMPDIR}}/non_relaxed.bigbed > {output.bg_file_non}

    """

rule gopeaks_narrow:
    input:
        bam = expand(join(RESULTSDIR,"bam","{replicate}.{dupstatus}.bam"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        align_stats = rules.gather_alignstats.output,
    params:
        gopeaks=TOOLS["gopeaks"],
        bam_path=join(RESULTSDIR,"bam"),
        tc_file="{treatment_control_list}",
        dupstatus = "{dupstatus}",
        qthresholds = "{qthresholds}",
        memG = getmemG("gopeaks")
    threads:
        getthreads("gopeaks")
    envmodules:
        TOOLS["ucsc"],
        TOOLS["samtools"]
    output:
        peak_file=join(RESULTSDIR,"peaks","gopeaks","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.narrow.peaks.bed"),
        bg_file=join(RESULTSDIR,"peaks","gopeaks","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.narrow.peaks.bigbed.gz"),
        json=temp(join(RESULTSDIR,"peaks","gopeaks","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.narrow.gopeaks.json")),
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

        # pull treatment and control ids
        treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
        control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`
            
        # set bam file
        treat_bam={params.bam_path}/${{treatment}}.{params.dupstatus}.bam
        cntrl_bam={params.bam_path}/${{control}}.{params.dupstatus}.bam

        # run gopeaks
        {params.gopeaks} -b ${{treat_bam}} -c ${{cntrl_bam}} -p {params.qthresholds} -o ${{TMPDIR}}/narrow

        # mv output and rename for consistency
        mv $TMPDIR/narrow_peaks.bed {output.peak_file}
        mv $TMPDIR/narrow_gopeaks.json {output.json}

        # create bigbed files
        cut -f1-3 {output.peak_file} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/narrow_peaks.bed 
        bedToBigBed -type=bed3 ${{TMPDIR}}/narrow_peaks.bed  {input.genome_len} ${{TMPDIR}}/narrow_peaks.bigbed
        bgzip -c ${{TMPDIR}}/narrow_peaks.bigbed > {output.bg_file}
        """

rule gopeaks_broad:
    input:
        bam = expand(join(RESULTSDIR,"bam","{replicate}.{dupstatus}.bam"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        align_stats = rules.gather_alignstats.output,
    params:
        gopeaks=TOOLS["gopeaks"],
        bam_path=join(RESULTSDIR,"bam"),
        tc_file="{treatment_control_list}",
        dupstatus = "{dupstatus}",
        qthresholds = "{qthresholds}",
        memG = getmemG("gopeaks")
    threads:
        getthreads("gopeaks")
    envmodules:
        TOOLS["ucsc"],
        TOOLS["samtools"]
    output:
        peak_file=join(RESULTSDIR,"peaks","gopeaks","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.broad.peaks.bed"),
        bg_file=join(RESULTSDIR,"peaks","gopeaks","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.broad.peaks.bigbed.gz"),
        json=temp(join(RESULTSDIR,"peaks","gopeaks","{qthresholds}","peak_output","{treatment_control_list}.{dupstatus}.broad.gopeaks.json")),
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

        # pull treatment and control ids
        treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
        control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`
            
        # set bam file
        treat_bam={params.bam_path}/${{treatment}}.{params.dupstatus}.bam
        cntrl_bam={params.bam_path}/${{control}}.{params.dupstatus}.bam

        # run gopeaks
        {params.gopeaks} -b ${{treat_bam}} -c ${{cntrl_bam}} -p {params.qthresholds} -o ${{TMPDIR}}/broad --broad

        # mv output and rename for consistency
        mv $TMPDIR/broad_peaks.bed {output.peak_file}
        mv $TMPDIR/broad_gopeaks.json {output.json}

        # create bigbed files
        cut -f1-3 {output.peak_file} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/broad.bed 
        bedToBigBed -type=bed3 ${{TMPDIR}}/broad.bed  {input.genome_len} ${{TMPDIR}}/broad.bigbed
        bgzip -c ${{TMPDIR}}/broad.bigbed > {output.bg_file}
        """