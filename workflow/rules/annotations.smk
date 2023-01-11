def get_bed_macs2(wildcards):
        if wildcards.macs_types == "narrowPeak":
                bed=join(RESULTSDIR,"peaks","macs2", wildcards.replicate, wildcards.replicate + "." + wildcards.dupstatus + "_peaks.narrowPeak")
        if wildcards.macs_types =="broadPeak":
                bed=join(RESULTSDIR,"peaks","macs2", wildcards.replicate, wildcards.replicate + "." + wildcards.dupstatus + "_peaks.broadPeak")
        return bed

def get_bed_s_and_g(wildcards):
        if wildcards.s_and_g_types =="norm.stringent.bed":
                bed=join(RESULTSDIR,"peaks","seacr", wildcards.t_and_c, wildcards.t_and_c + "." + wildcards.dupstatus + ".norm.stringent.bed")
        if wildcards.s_and_g_types =="norm.relaxed.bed":
                bed=join(RESULTSDIR,"peaks","seacr", wildcards.t_and_c, wildcards.t_and_c + "." + wildcards.dupstatus + ".norm.relaxed.bed")
        if wildcards.s_and_g_types =="narrowGo_peaks.bed":
                bed=join(RESULTSDIR,"peaks","gopeaks", wildcards.t_and_c + "." + wildcards.dupstatus + ".narrowGo_peaks.bed")
        if wildcards.s_and_g_types =="broadGo_peaks.bed":
                bed=join(RESULTSDIR,"peaks","gopeaks", wildcards.t_and_c + "." + wildcards.dupstatus + ".broadGo_peaks.bed")
        return bed

def get_cntrl_bam(wildcards):
        cntrl_sample=TREAT_to_CONTRL_DICT[wildcards.treatment]
        cntrl_file=join(RESULTSDIR,"bam", cntrl_sample + "." + wildcards.dupstatus + ".bam"),
        return cntrl_file

def get_bed_all(wildcards):
        cntrl=TREAT_to_CONTRL_DICT[wildcards.treatment]
        if wildcards.peak_types == "narrowPeak":
                bed=join(RESULTSDIR,"peaks","macs2", wildcards.treatment, wildcards.treatment + "." + wildcards.dupstatus + "_peaks.narrowPeak")
        if wildcards.peak_types =="broadPeak":
                bed=join(RESULTSDIR,"peaks","macs2", wildcards.treatment, wildcards.treatment + "." + wildcards.dupstatus + "_peaks.broadPeak")
        if wildcards.peak_types =="norm.stringent.bed":
                bed=join(RESULTSDIR,"peaks","seacr", wildcards.treatment + "_vs_" + cntrl, wildcards.treatment + "_vs_" + cntrl + "." + wildcards.dupstatus + ".norm.stringent.bed")
        if wildcards.peak_types =="norm.relaxed.bed":
                bed=join(RESULTSDIR,"peaks","seacr", wildcards.treatment + "_vs_" + cntrl, wildcards.treatment + "_vs_" + cntrl + "." + wildcards.dupstatus + ".norm.relaxed.bed")
        if wildcards.peak_types =="narrowGo_peaks.bed":
                bed=join(RESULTSDIR,"peaks","gopeaks", wildcards.treatment + "_vs_" + cntrl + "." + wildcards.dupstatus + ".narrowGo_peaks.bed")
        if wildcards.peak_types =="broadGo_peaks.bed":
                bed=join(RESULTSDIR,"peaks","gopeaks", wildcards.treatment + "_vs_" + cntrl + "." + wildcards.dupstatus + ".broadGo_peaks.bed")
        return bed

rule peakAnnotation_macs2:
        """
        Developed from code: https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk
        """
        input:
                get_bed_macs2,
        output:
                annotation=join(RESULTSDIR,"annotation","{macs_types}","{replicate}","{replicate}.{dupstatus}.annotation.txt"),
                annotation_summary=join(RESULTSDIR,"annotation","{macs_types}","{replicate}","{replicate}.{dupstatus}.annotation.summary")
        envmodules:
                TOOLS["homer"],
        params:
                genome = config["genome"],
        shell:
                """  
                annotatePeaks.pl {input} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
                """

rule peakAnnotation_s_and_g:
        """
        Developed from code: https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk
        """
        input:
                get_bed_s_and_g,
        output:
                annotation=join(RESULTSDIR,"annotation","{s_and_g_types}","{t_and_c}","{t_and_c}.{dupstatus}.annotation.txt"),
                annotation_summary=join(RESULTSDIR,"annotation","{s_and_g_types}","{t_and_c}","{t_and_c}.{dupstatus}.annotation.summary"),
        envmodules:
                TOOLS["homer"],
        params:
                genome = config["genome"],
        shell:
                """  
                annotatePeaks.pl {input} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
                """

rule rose_all:
        """
        Developed from code: https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk
        """
        input:
                bed=get_bed_all,
                bam=join(RESULTSDIR,"bam", "{treatment}.{dupstatus}.bam"),
                cntrl_bam=get_cntrl_bam,
        output:
                no_tss_bed=join(RESULTSDIR,"annotation","{peak_types}","{treatment}","{treatment}.{dupstatus}.no_TSS_{s_dist}.bed"),
                all=join(RESULTSDIR,"annotation","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllEnhancers.table.txt"),
                regular=join(RESULTSDIR,"annotation","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}AllEnhancers.table.regular.bed"),
                super=join(RESULTSDIR,"annotation","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllEnhancers.table.super.bed"),
                regular_great=join(RESULTSDIR,"annotation","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllEnhancers.table.super.GREAT.bed"),
                super_great=join(RESULTSDIR,"annotation","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllEnhancers.table.regular.GREAT.bed"),
                regular_summit=join(RESULTSDIR,"annotation","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllEnhancers.table.regular.summits.bed"),
                super_summit=join(RESULTSDIR,"annotation","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllEnhancers.table.super.summits.bed"),
        envmodules:
                TOOLS["homer"],
                TOOLS["bedtools"],
                TOOLS["rose"],
                TOOLS["python37"],
                TOOLS["R"],
        params:
                genome = config["genome"],
                tss_bed = config["reference"][config["genome"]]["tss_bed"],
                stitch_distance = config["stitch_distance"],
                tss_distance=config["tss_distance"],
                prefix=join(RESULTSDIR,"annotation","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}"),
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
                
                # bedtools
                bedtools intersect -a {input.bed} -b {params.tss_bed} -v > $TMPDIR/tmp.bed
                bedtools merge -i $TMPDIR/tmp.bed -d {params.stitch_distance} -c 4,5,6 -o distinct,sum,distinct > {output.no_tss_bed}

                # rose
                ROSE_main.py \
                        -i {output.no_tss_bed} \
                        -g {params.genome} \
                        -r {input.bam} {input.cntrl_bam} \
                        -t {params.tss_distance} \
                        -s {params.stitch_distance} \
                        -o {params.prefix}
                
                # rose to bed file
                # developed from https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/scripts/roseTable2Bed.sh
                grep -v "^[#|REGION]" {output.all} | awk -v OFS="\\t" -F"\\t" \'$NF==0 {{for(i=2; i<=NF; i++){{printf $i; printf (i<NF?"\\t":"\\n")}}}}\' > $TMPDIR/regular
                bedtools sort -i $TMPDIR/regular > {output.regular}
                
                grep -v "^[#|REGION]" {output.all} | awk -v OFS="\\t" -F"\\t" \'$NF==1 {{for(i=2; i<=NF; i++){{printf $i; printf (i<NF?"\\t":"\\n")}}}}\' > $TMPDIR/super
                bedtools sort -i $TMPDIR/super > {output.super}
                
                # cut rose output files
                cut -f1-3 {output.regular} > {output.regular_great}
                cut -f1-3 {output.super} > {output.super_great}

                # summit files
                bedtools intersect -wa -a {input.bed} -b {output.regular} > {output.regular_summit}
                bedtools intersect -wa -a {input.bed} -b {output.super} > {output.super_summit}
                """

rule findMotif_all:
        """
        Developed from code: https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk
        """
        input:
                bed=get_bed_all,
        output:
                known_html=join(RESULTSDIR,"annotation","{peak_types}","{treatment}","{treatment}.{dupstatus}.motifs","knownResults.html"),
        threads: getthreads("findMotif_all")
        envmodules:
                TOOLS["homer"],
        params:
                genome = config["genome"],
                outDir=join(RESULTSDIR,"annotation","{peak_types}","{treatment}","{treatment}.{dupstatus}.motifs"),
                motif_size = config["motif_size"],
                preparsedDir = config["preparsedDir"],
        shell:
                """
                findMotifsGenome.pl {input.bed} {params.genome} {params.outDir} -size {params.motif_size} -p {threads} -preparsedDir {params.preparsedDir}
                """