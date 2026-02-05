# TSV file should include 6 columns
# 1)condition     2)sample
# 53_H3K4me3	    53_H3K4me3_1
# 3)bed
# /results/peaks/gopeaks/53_H3K4me3_1_vs_HN6_IgG_rabbit_negative_control_1.dedup.broadGo_peaks.bed
# 4)bedgraph                                      5)scaling factor
# /results/bedgraph/53_H3K4me3_1.dedup.bedgraph	86.32596685082872928176
# 6)bed
# /results/fragments/53_H3K4me3_1.dedup.fragments.bed

localrules: count_peaks

import re

rule merge_control_bams:
    """
    Merge all replicates of a control sample for pooled control analysis
    """
    input:
        bams = lambda w: expand(join(RESULTSDIR,"bam","{replicate}.{dupstatus}.bam"),
                               replicate=[r for r in REPLICATES if r.startswith(w.control_sample+"_")],
                               dupstatus=w.dupstatus)
    output:
        merged_bam = join(RESULTSDIR,"bam","pooled_controls","{control_sample}.{dupstatus}.merged.bam"),
        merged_bai = join(RESULTSDIR,"bam","pooled_controls","{control_sample}.{dupstatus}.merged.bam.bai")
    params:
        bam_list = lambda w, input: " ".join(input.bams)
    threads: getthreads("samtools")
    envmodules:
        TOOLS["samtools"]
    shell:
        """
        set -exo pipefail
        
        # Merge BAM files
        if [[ {threads} -gt 1 ]]; then
            samtools merge -@ {threads} {output.merged_bam} {params.bam_list}
        else
            samtools merge {output.merged_bam} {params.bam_list}
        fi
        
        # Index merged BAM
        samtools index {output.merged_bam}
        """

rule create_pooled_control_fragments:
    """
    Create fragments.bed file from pooled control BAM
    """
    input:
        merged_bam = join(RESULTSDIR,"bam","pooled_controls","{control_sample}.{dupstatus}.merged.bam")
    output:
        fragments = join(RESULTSDIR,"fragments","pooled_controls","{control_sample}.{dupstatus}.fragments.bed")
    params:
        fragment_len = config["fragment_len_filter"],
        mapping_quality = config["mapping_quality"],
        control_sample = "{control_sample}",
        dupstatus = "{dupstatus}",
        memG = getmemG("create_pooled_control_fragments")
    threads: getthreads("create_pooled_control_fragments")
    envmodules:
        TOOLS["bedtools"],
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
        
        # Name-sort BAM (required for bedtools bamtobed -bedpe)
        samtools sort -n -@ {threads} -T $TMPDIR -o ${{TMPDIR}}/{params.control_sample}.{params.dupstatus}.namesorted.bam {input.merged_bam}
        
        # Convert BAM to fragments bed
        bedtools bamtobed -bedpe -i ${{TMPDIR}}/{params.control_sample}.{params.dupstatus}.namesorted.bam | \
        awk -v OFS='\\t' -v flen={params.fragment_len} -v mapq={params.mapping_quality} \
        '$1 == $4 && $6-$2 < flen && $8 >= mapq {{print $1,$2,$6}}' | \
        LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n > {output.fragments}
        """

rule create_pooled_control_bedgraph:
    """
    Create bedgraph file from pooled control fragments for SEACR and GoPeaks
    """
    input:
        fragments = join(RESULTSDIR,"fragments","pooled_controls","{control_sample}.{dupstatus}.fragments.bed"),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        align_stats = rules.gather_alignstats.output,
    output:
        bedgraph = join(RESULTSDIR,"bedgraph","pooled_controls","{control_sample}.{dupstatus}.bedgraph")
    params:
        norm_method = NORM_METHOD,
        scale = config["spikein_scale"] if NORM_METHOD == "SPIKEIN" else None,
        control_sample = "{control_sample}",
        dupstatus = "{dupstatus}"
    threads: 1
    envmodules:
        TOOLS["bedtools"]
    shell:
        """
        set -exo pipefail
        
        # Get scaling factor from alignment stats
        # For pooled controls, use average of all control replicate scaling factors
        if [[ {params.norm_method} == "LIBRARY" ]]; then
            # Calculate average library scaling factor across control replicates
            scale=$(awk -F'\\t' '$1 ~ /^{params.control_sample}_[0-9]+$/ && $2 == "{params.dupstatus}" {{sum+=$NF; count++}} END {{if(count>0) print sum/count; else print 1}}' {input.align_stats})
        elif [[ {params.norm_method} == "SPIKEIN" ]]; then
            # Calculate average spikein scaling factor across control replicates  
            scale=$(awk -F'\\t' '$1 ~ /^{params.control_sample}_[0-9]+$/ && $2 == "{params.dupstatus}" {{sum+=$(NF-1); count++}} END {{if(count>0) print sum/count; else print 1}}' {input.align_stats})
        else
            scale=1
        fi
        
        # Create bedgraph with scaling
        bedtools genomecov -bg -scale $scale -i {input.fragments} -g {input.genome_len} | \
        sort -k1,1 -k2,2n > {output.bedgraph}
        """

def get_all_peak_files(wildcards):
    files=[]
    pool_controls = config.get("pool_controls", False)

    def control_has_replicate(pair: str) -> bool:
        parts = pair.split("_vs_")
        if len(parts) != 2:
            return False
        ctrl = parts[1]
        return re.search(r"_[0-9]+$", ctrl) is not None

    # Determine control modes to run
    control_modes = ["individual"]
    if pool_controls:
        control_modes.append("pooled")

    for control_mode in control_modes:
        # Use appropriate treatment list based on control mode
        tc_list_m_base = TREATMENT_CONTROL_LIST_POOLED if control_mode == "pooled" else TREATMENT_LIST_M
        tc_list_sg_base = TREATMENT_CONTROL_LIST_POOLED if control_mode == "pooled" else TREATMENT_LIST_SG

        # Filter pairs to enforce correct control suffix per mode
        if control_mode == "pooled":
            tc_list_m = [p for p in tc_list_m_base if not control_has_replicate(p)]
            tc_list_sg = [p for p in tc_list_sg_base if not control_has_replicate(p)]
        else:
            tc_list_m = [p for p in tc_list_m_base if control_has_replicate(p)]
            tc_list_sg = [p for p in tc_list_sg_base if control_has_replicate(p)]

        if "macs2_narrow" in PEAKTYPE:
            n=expand(join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.narrow.peaks.bed"),
                    qthresholds=QTRESHOLDS,treatment_control_list=tc_list_m,dupstatus=DUPSTATUS,control_mode=control_mode)
            files.extend(n)
        if "macs2_broad" in PEAKTYPE:
            b=expand(join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.broad.peaks.bed"),
                    qthresholds=QTRESHOLDS,treatment_control_list=tc_list_m,dupstatus=DUPSTATUS,control_mode=control_mode)
            files.extend(b)
        if "seacr_stringent" in PEAKTYPE:
            s=expand(join(RESULTSDIR,"peaks","{qthresholds}","seacr","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.stringent.peaks.bed"),
                    qthresholds=QTRESHOLDS,treatment_control_list=tc_list_sg,dupstatus=DUPSTATUS,control_mode=control_mode)
            files.extend(s)
        if "seacr_relaxed" in PEAKTYPE:
            r=expand(join(RESULTSDIR,"peaks","{qthresholds}","seacr","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.relaxed.peaks.bed"),
                    qthresholds=QTRESHOLDS,treatment_control_list=tc_list_sg,dupstatus=DUPSTATUS,control_mode=control_mode)
            files.extend(r)
        if "gopeaks_narrow" in PEAKTYPE:
            n=expand(join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.narrow.peaks.bed"),
                    qthresholds=QTRESHOLDS,treatment_control_list=tc_list_sg,dupstatus=DUPSTATUS,control_mode=control_mode)
            files.extend(n)
        if "gopeaks_broad" in PEAKTYPE:
            b=expand(join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.broad.peaks.bed"),
                    qthresholds=QTRESHOLDS,treatment_control_list=tc_list_sg,dupstatus=DUPSTATUS,control_mode=control_mode)
            files.extend(b)

    return files

rule macs2_narrow:
    '''
    MACS2 can be run with and without a control. This featured is controlled in the config.yaml file
    with the flag macs2_control: "N" or macs2_control: "Y"
    The control_mode wildcard determines if individual or pooled controls are used.
    When control_mode="individual", treatment_control_list must match individual patterns (both ends with _N).
    When control_mode="pooled", treatment_control_list must match pooled patterns (control has no _N suffix).
    '''
    input:
        fragments_bed = expand(join(RESULTSDIR,"fragments","{replicate}.{dupstatus}.fragments.bed"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        pooled_controls = expand(join(RESULTSDIR,"fragments","pooled_controls","{control_sample}.{dupstatus}.fragments.bed"),
                                control_sample=CONTROL_SAMPLES,
                                dupstatus=DUPSTATUS) if config.get("pool_controls", False) else [],
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
    output:
        peak_file = join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.narrow.peaks.bed"),
        summit_file = join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.narrow.summits.bed"),
        bg_file = join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.narrow.peaks.bigbed.gz"),
        xls_file = join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.macs2_narrow.peaks.xls"),
    params:
        frag_bed_path=join(RESULTSDIR,"fragments"),
        pooled_frag_path=join(RESULTSDIR,"fragments","pooled_controls"),
        tc_file="{treatment_control_list}",
        control_mode="{control_mode}",
        dupstatus = "{dupstatus}",
        qthresholds = "{qthresholds}",
        control_flag = config["macs2_control"],
        pool_controls = config.get("pool_controls", False),
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

        # Validate control_mode and treatment_control_list combination
        tc_pair="{params.tc_file}"
        control_mode="{params.control_mode}"
        control=`echo "$tc_pair" | awk -F"_vs_" '{{print $NF}}'`
        
        if [[ "$control_mode" == "pooled" ]]; then
            # Pooled mode: control should NOT end with _N pattern
            if [[ "$control" =~ _[0-9]+$ ]]; then
                echo "ERROR: Invalid pooled control pair '$tc_pair' - control has replicate number"
                echo "       Pooled pairs must be like: TREATMENT_X_vs_CONTROL (no _N on control)"
                echo "       Found: $control with _N suffix"
                exit 1
            fi
        else
            # Individual mode: control MUST end with _N pattern
            if [[ ! "$control" =~ _[0-9]+$ ]]; then
                echo "ERROR: Invalid individual control pair '$tc_pair' - control missing replicate number"
                echo "       Individual pairs must be like: TREATMENT_X_vs_CONTROL_X"
                echo "       Found: $control without _N suffix"
                exit 1
            fi
        fi

        # pull treatment and control ids
        treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
        control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`

        # set treatment fragment file
        treat_bed={params.frag_bed_path}/${{treatment}}.{params.dupstatus}.fragments.bed
        
        # set control fragment file based on control_mode wildcard
        if [[ "{params.control_mode}" == "pooled" ]]; then
            # Use pooled control (remove replicate number from control name)
            control_sample=`echo ${{control}} | sed 's/_[0-9]*$//'`
            cntrl_bed={params.pooled_frag_path}/${{control_sample}}.{params.dupstatus}.fragments.bed
        else
            # Use per-replicate control (individual mode)
            control_sample=${{control}}
            cntrl_bed={params.frag_bed_path}/${{control}}.{params.dupstatus}.fragments.bed
        fi

        # run with or without control
        if [[ {params.control_flag} == "Y" ]]; then
            file_base="${{treatment}}_vs_${{control_sample}}.{params.dupstatus}"

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
            file_base="${{treatment}}_vs_${{control_sample}}.{params.dupstatus}"

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
        mv $TMPDIR/${{file_base}}_summits.bed {output.summit_file}
        mv $TMPDIR/${{file_base}}_peaks.xls {output.xls_file}

        # create bigbed files, zip
        cut -f1-3 {output.peak_file} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/narrow.bed
        bedToBigBed -type=bed3 ${{TMPDIR}}/narrow.bed {input.genome_len} ${{TMPDIR}}/narrow.bigbed
        bgzip -c ${{TMPDIR}}/narrow.bigbed > {output.bg_file}
        """

rule macs2_broad:
    '''
    MACS2 can be run with and without a control. This featured is controlled in the config.yaml file
    with the flag macs2_control: "N" or macs2_control: "Y"
    The control_mode wildcard determines if individual or pooled controls are used.
    '''
    input:
        fragments_bed = expand(join(RESULTSDIR,"fragments","{replicate}.{dupstatus}.fragments.bed"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        pooled_controls = expand(join(RESULTSDIR,"fragments","pooled_controls","{control_sample}.{dupstatus}.fragments.bed"),
                                control_sample=CONTROL_SAMPLES,
                                dupstatus=DUPSTATUS) if config.get("pool_controls", False) else [],
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        align_stats = rules.gather_alignstats.output,
    output:
        peak_file = join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.broad.peaks.bed"),
        bg_file = join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.broad.peaks.bigbed.gz"),
        xls_file = join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.macs2_broad.peaks.xls"),
    params:
        frag_bed_path=join(RESULTSDIR,"fragments"),
        pooled_frag_path=join(RESULTSDIR,"fragments","pooled_controls"),
        tc_file="{treatment_control_list}",
        control_mode="{control_mode}",
        dupstatus = "{dupstatus}",
        qthresholds = "{qthresholds}",
        control_flag = config["macs2_control"],
        pool_controls = config.get("pool_controls", False),
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

        # set treatment fragment file
        treat_bed={params.frag_bed_path}/${{treatment}}.{params.dupstatus}.fragments.bed
        
        # set control fragment file based on control_mode wildcard
        if [[ "{params.control_mode}" == "pooled" ]]; then
            # Use pooled control (remove replicate number from control name)
            control_sample=`echo ${{control}} | sed 's/_[0-9]*$//'`
            cntrl_bed={params.pooled_frag_path}/${{control_sample}}.{params.dupstatus}.fragments.bed
        else
            # Use per-replicate control (individual mode)
            control_sample=${{control}}
            cntrl_bed={params.frag_bed_path}/${{control}}.{params.dupstatus}.fragments.bed
        fi

        # run with or without control
        if [[ {params.control_flag} == "Y" ]]; then
            file_base="${{treatment}}_vs_${{control_sample}}.{params.dupstatus}"

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
        mv $TMPDIR/${{file_base}}_peaks.xls {output.xls_file}

        # create bigbed files
        cut -f1-3 {output.peak_file} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/broad.bed
        bedToBigBed -type=bed3 ${{TMPDIR}}/broad.bed {input.genome_len} ${{TMPDIR}}/broad.bigbed
        bgzip -c ${{TMPDIR}}/broad.bigbed > {output.bg_file}
        """

rule seacr_stringent:
    input:
        tc_list=join(RESULTSDIR,"treatment_control_list.txt"),
        bg = expand(join(RESULTSDIR,"bedgraph","{replicate}.{dupstatus}.bedgraph"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        pooled_controls = expand(join(RESULTSDIR,"bedgraph","pooled_controls","{control_sample}.{dupstatus}.bedgraph"),
                                 control_sample=CONTROL_SAMPLES, dupstatus=DUPSTATUS) if config.get("pool_controls", False) else [],
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        align_stats = rules.gather_alignstats.output,
    output:
        peak_file = join(RESULTSDIR,"peaks","{qthresholds}","seacr","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.stringent.peaks.bed"),
        bg_file = join(RESULTSDIR,"peaks","{qthresholds}","seacr","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.stringent.peaks.bigbed.gz"),
    params:
        bg_path=join(RESULTSDIR,"bedgraph"),
        pooled_bg_path=join(RESULTSDIR,"bedgraph","pooled_controls"),
        frag_bed_path=join(RESULTSDIR,"fragments"),
        tc_file="{treatment_control_list}",
        control_mode="{control_mode}",
        qthresholds="{qthresholds}",
        dupstatus= "{dupstatus}",
        memG = getmemG("seacr_stringent"),
        norm_method = NORM_METHOD,
        pool_controls = config.get("pool_controls", False)
    threads: getthreads("seacr_stringent")
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

        # Validate control_mode and treatment_control_list combination
        tc_pair="{params.tc_file}"
        control_mode="{params.control_mode}"
        control=`echo "$tc_pair" | awk -F"_vs_" '{{print $NF}}'`
        
        if [[ "$control_mode" == "pooled" ]]; then
            if [[ "$control" =~ _[0-9]+$ ]]; then
                echo "ERROR: Invalid pooled control pair '$tc_pair' - control has replicate number"
                exit 1
            fi
        else
            if [[ ! "$control" =~ _[0-9]+$ ]]; then
                echo "ERROR: Invalid individual control pair '$tc_pair' - control missing replicate number"
                exit 1
            fi
        fi

        # pull treatment and control ids
        treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
        control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`

        # Check if using pooled controls based on control_mode wildcard
        if [[ "{params.control_mode}" == "pooled" ]]; then
            # Extract control sample name (remove replicate suffix)
            control_sample=$(echo $control | sed 's/_[0-9]*$//')
            cntrl_bed={params.pooled_bg_path}/${{control_sample}}.{params.dupstatus}.bedgraph
        else
            # Use per-replicate control (individual mode)
            cntrl_bed={params.bg_path}/${{control}}.{params.dupstatus}.bedgraph
        fi

        # set treatment bedgraph file
        treat_bed={params.bg_path}/${{treatment}}.{params.dupstatus}.bedgraph

        if [[ {params.norm_method} != "SPIKEIN" ]]; then
            # run norm
            SEACR.sh --bedgraph ${{treat_bed}} \\
                --control ${{cntrl_bed}} \\
                --normalize norm \\
                --mode stringent \\
                --threshold {params.qthresholds} \\
                --output norm

            # mv output and rename for consistency
            mv $TMPDIR/norm.stringent.bed {output.peak_file}
        else
            # run non
            SEACR.sh --bedgraph ${{treat_bed}} \\
                --control ${{cntrl_bed}} \\
                --normalize non \\
                --mode stringent \\
                --threshold {params.qthresholds} \\
                --output non

            # mv output and rename for consistency
            mv $TMPDIR/non.stringent.bed {output.peak_file}
        fi

        # create bigbed files
        cut -f1-3 {output.peak_file} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/stringent.bed
        bedToBigBed -type=bed3 ${{TMPDIR}}/stringent.bed {input.genome_len} ${{TMPDIR}}/stringent.bigbed
        bgzip -c ${{TMPDIR}}/stringent.bigbed > {output.bg_file}
        """

rule seacr_relaxed:
    input:
        bg = expand(join(RESULTSDIR,"bedgraph","{replicate}.{dupstatus}.bedgraph"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        pooled_controls = expand(join(RESULTSDIR,"bedgraph","pooled_controls","{control_sample}.{dupstatus}.bedgraph"),
                                 control_sample=CONTROL_SAMPLES, dupstatus=DUPSTATUS) if config.get("pool_controls", False) else [],
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        align_stats = rules.gather_alignstats.output,
    output:
        peak_file = join(RESULTSDIR,"peaks","{qthresholds}","seacr","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.relaxed.peaks.bed"),
        bg_file = join(RESULTSDIR,"peaks","{qthresholds}","seacr","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.relaxed.peaks.bigbed.gz"),
    params:
        bg_path=join(RESULTSDIR,"bedgraph"),
        pooled_bg_path=join(RESULTSDIR,"bedgraph","pooled_controls"),
        frag_bed_path=join(RESULTSDIR,"fragments"),
        tc_file="{treatment_control_list}",
        control_mode="{control_mode}",
        qthresholds="{qthresholds}",
        dupstatus= "{dupstatus}",
        memG = getmemG("seacr_relaxed"),
        norm_method = NORM_METHOD,
        pool_controls = config.get("pool_controls", False)
    threads: getthreads("seacr_relaxed")
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

            # Validate control_mode and treatment_control_list combination
            tc_pair="{params.tc_file}"
            control_mode="{params.control_mode}"
            control=`echo "$tc_pair" | awk -F"_vs_" '{{print $NF}}'`
            
            if [[ "$control_mode" == "pooled" ]]; then
                if [[ "$control" =~ _[0-9]+$ ]]; then
                    echo "ERROR: Invalid pooled control pair '$tc_pair' - control has replicate number"
                    exit 1
                fi
            else
                if [[ ! "$control" =~ _[0-9]+$ ]]; then
                    echo "ERROR: Invalid individual control pair '$tc_pair' - control missing replicate number"
                    exit 1
                fi
            fi

            # pull treatment and control ids
            treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
            control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`

            # Check if using pooled controls based on control_mode wildcard
            if [[ "{params.control_mode}" == "pooled" ]]; then
                # Extract control sample name (remove replicate suffix)
                control_sample=$(echo $control | sed 's/_[0-9]*$//')
                cntrl_bed={params.pooled_bg_path}/${{control_sample}}.{params.dupstatus}.bedgraph
            else
                # Use per-replicate control (individual mode)
                cntrl_bed={params.bg_path}/${{control}}.{params.dupstatus}.bedgraph
            fi

            # set treatment bedgraph file
            treat_bed={params.bg_path}/${{treatment}}.{params.dupstatus}.bedgraph

            if [[ {params.norm_method} != "SPIKEIN" ]]; then
                # run norm
                SEACR.sh --bedgraph  ${{treat_bed}} \\
                    --control ${{cntrl_bed}} \\
                    --normalize norm \\
                    --mode relaxed \\
                    --threshold {params.qthresholds} \\
                    --output norm

                # mv output and rename for consistency
                mv $TMPDIR/norm.relaxed.bed {output.peak_file}
            else
                # run non
                SEACR.sh --bedgraph ${{treat_bed}} \\
                    --control ${{cntrl_bed}} \\
                    --normalize non \\
                    --mode relaxed \\
                    --threshold {params.qthresholds} \\
                    --output non

                mv $TMPDIR/non.relaxed.bed {output.peak_file}
            fi

            # create bigbed files
            cut -f1-3 {output.peak_file} | LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n | uniq > ${{TMPDIR}}/relaxed.bed
            bedToBigBed -type=bed3 ${{TMPDIR}}/relaxed.bed {input.genome_len} ${{TMPDIR}}/relaxed.bigbed
            bgzip -c ${{TMPDIR}}/relaxed.bigbed > {output.bg_file}
    """

rule gopeaks_narrow:
    input:
        bam = expand(join(RESULTSDIR,"bam","{replicate}.{dupstatus}.bam"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        pooled_controls = expand(join(RESULTSDIR,"bam","pooled_controls","{control_sample}.{dupstatus}.merged.bam"),
                                 control_sample=CONTROL_SAMPLES, dupstatus=DUPSTATUS) if config.get("pool_controls", False) else [],
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        align_stats = rules.gather_alignstats.output,
    params:
        gopeaks=TOOLS["gopeaks"],
        bam_path=join(RESULTSDIR,"bam"),
        pooled_bam_path=join(RESULTSDIR,"bam","pooled_controls"),
        tc_file="{treatment_control_list}",
        control_mode="{control_mode}",
        dupstatus = "{dupstatus}",
        qthresholds = "{qthresholds}",
        memG = getmemG("gopeaks_narrow"),
        pool_controls = config.get("pool_controls", False)
    threads:
        getthreads("gopeaks_narrow")
    envmodules:
        TOOLS["ucsc"],
        TOOLS["samtools"]
    output:
        peak_file=join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.narrow.peaks.bed"),
        bg_file=join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.narrow.peaks.bigbed.gz"),
        json=temp(join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","narrow","{control_mode}","{treatment_control_list}.{dupstatus}.narrow.gopeaks.json")),
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

        # Validate control_mode and treatment_control_list combination
        tc_pair="{params.tc_file}"
        control_mode="{params.control_mode}"
        control=`echo "$tc_pair" | awk -F"_vs_" '{{print $NF}}'`
        
        if [[ "$control_mode" == "pooled" ]]; then
            if [[ "$control" =~ _[0-9]+$ ]]; then
                echo "ERROR: Invalid pooled control pair '$tc_pair' - control has replicate number"
                exit 1
            fi
        else
            if [[ ! "$control" =~ _[0-9]+$ ]]; then
                echo "ERROR: Invalid individual control pair '$tc_pair' - control missing replicate number"
                exit 1
            fi
        fi

        # pull treatment and control ids
        treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
        control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`

        # Check if using pooled controls based on control_mode wildcard
        if [[ "{params.control_mode}" == "pooled" ]]; then
            # Extract control sample name (remove replicate suffix)
            control_sample=$(echo $control | sed 's/_[0-9]*$//')
            cntrl_bam={params.pooled_bam_path}/${{control_sample}}.{params.dupstatus}.merged.bam
        else
            # Use per-replicate control (individual mode)
            cntrl_bam={params.bam_path}/${{control}}.{params.dupstatus}.bam
        fi

        # set treatment bam file
        treat_bam={params.bam_path}/${{treatment}}.{params.dupstatus}.bam

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
        pooled_controls = expand(join(RESULTSDIR,"bam","pooled_controls","{control_sample}.{dupstatus}.merged.bam"),
                                 control_sample=CONTROL_SAMPLES, dupstatus=DUPSTATUS) if config.get("pool_controls", False) else [],
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        align_stats = rules.gather_alignstats.output,
    params:
        gopeaks=TOOLS["gopeaks"],
        bam_path=join(RESULTSDIR,"bam"),
        pooled_bam_path=join(RESULTSDIR,"bam","pooled_controls"),
        tc_file="{treatment_control_list}",
        control_mode="{control_mode}",
        dupstatus = "{dupstatus}",
        qthresholds = "{qthresholds}",
        memG = getmemG("gopeaks_broad"),
        pool_controls = config.get("pool_controls", False)
    threads:
        getthreads("gopeaks_broad")
    envmodules:
        TOOLS["ucsc"],
        TOOLS["samtools"]
    output:
        peak_file=join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.broad.peaks.bed"),
        bg_file=join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.broad.peaks.bigbed.gz"),
        json=temp(join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","broad","{control_mode}","{treatment_control_list}.{dupstatus}.broad.gopeaks.json")),
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

        # Validate control_mode and treatment_control_list combination
        tc_pair="{params.tc_file}"
        control_mode="{params.control_mode}"
        control=`echo "$tc_pair" | awk -F"_vs_" '{{print $NF}}'`
        
        if [[ "$control_mode" == "pooled" ]]; then
            if [[ "$control" =~ _[0-9]+$ ]]; then
                echo "ERROR: Invalid pooled control pair '$tc_pair' - control has replicate number"
                exit 1
            fi
        else
            if [[ ! "$control" =~ _[0-9]+$ ]]; then
                echo "ERROR: Invalid individual control pair '$tc_pair' - control missing replicate number"
                exit 1
            fi
        fi

        # pull treatment and control ids
        treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
        control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`

        # Check if using pooled controls based on control_mode wildcard
        if [[ "{params.control_mode}" == "pooled" ]]; then
            # Extract control sample name (remove replicate suffix)
            control_sample=$(echo $control | sed 's/_[0-9]*$//')
            cntrl_bam={params.pooled_bam_path}/${{control_sample}}.{params.dupstatus}.merged.bam
        else
            # Use per-replicate control (individual mode)
            cntrl_bam={params.bam_path}/${{control}}.{params.dupstatus}.bam
        fi

        # set treatment bam file
        treat_bam={params.bam_path}/${{treatment}}.{params.dupstatus}.bam

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

rule count_peaks:
    input:
        peaks=get_all_peak_files
    output:
        peak_count=join(RESULTSDIR,"peaks","all.peaks.txt"),
        peak_table=join(RESULTSDIR,"peaks","Peak counts.xlsx"),
    params:
        outdir=join(RESULTSDIR,"peaks"),
        rscript=join(SCRIPTSDIR,"_plot_peak_counts.R")
    container: config['containers']['carlisle_r']
    shell:
        """
        wc -l {input.peaks} > {output.peak_count}
        Rscript {params.rscript} {output.peak_count} {params.outdir}
        """
