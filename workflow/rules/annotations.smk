def get_bed_macs2(wildcards):
        if wildcards.macs_types == "narrowPeak":
                bed=join(RESULTSDIR,"peaks","macs2", wildcards.replicate, wildcards.replicate + "." + wildcards.dupstatus + "_peaks.narrowPeak")
        if wildcards.macs_types =="broadPeak":
                bed=join(RESULTSDIR,"peaks","macs2", wildcards.replicate, wildcards.replicate + "." + wildcards.dupstatus + "_peaks.broadPeak")
        return bed

def get_bed2_macs2(wildcards):
        if wildcards.macs_types == "narrowPeak":
                bed=join(RESULTSDIR,"peaks","macs2", wildcards.treatment, wildcards.treatment + "." + wildcards.dupstatus + "_peaks.narrowPeak")
        if wildcards.macs_types =="broadPeak":
                bed=join(RESULTSDIR,"peaks","macs2", wildcards.treatment, wildcards.treatment + "." + wildcards.dupstatus + "_peaks.broadPeak")
        return bed

def get_cntrl_bam(wildcards):
        cntrl_sample=TREAT_to_CONTRL_DICT[wildcards.treatment]
        cntrl_file=join(RESULTSDIR,"bam", cntrl_sample + "." + wildcards.dupstatus + ".bam"),
        return cntrl_file

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
        input:
                get_bed_macs2,
        output:
                annotation=join(RESULTSDIR,"annotation","{macs_types}","{replicate}.{dupstatus}_peaks.annotation.txt"),
                annotation_summary=join(RESULTSDIR,"annotation","{macs_types}","{replicate}.{dupstatus}_peaks.annotation.summary")
        envmodules:
                TOOLS["homer"],
        params:
                genome = config["genome"],
        shell:
                """  
                annotatePeaks.pl {input} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
                """

rule peakAnnotation_s_and_g:
        input:
                get_bed_s_and_g,
        output:
                annotation=join(RESULTSDIR,"annotation","{s_and_g_types}","{t_and_c}.{dupstatus}_peaks.annotation.txt"),
                annotation_summary=join(RESULTSDIR,"annotation","{s_and_g_types}","{t_and_c}.{dupstatus}_peaks.annotation.summary"),
        envmodules:
                TOOLS["homer"],
        params:
                genome = config["genome"],
        shell:
                """  
                annotatePeaks.pl {input} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
                """

rule rose_all:
        input:
                bed=get_bed_all,
                bam=join(RESULTSDIR,"bam", "{treatment}.{dupstatus}.bam"),
                cntrl_bam=get_cntrl_bam,
        output:
                no_tss_bed=join(RESULTSDIR,"annotation","{peak_types}","{treatment}.{dupstatus}_peaks.no_TSS_{s_dist}.bed"),
                all=join(RESULTSDIR,"annotation","{peak_types}","{treatment}.{dupstatus}.{s_dist},{treatment}_peaks_AllEnhancers.table.txt"),
                regular=join(RESULTSDIR,"annotation","{peak_types}","{treatment}.{dupstatus}.{s_dist},{treatment}_peaks_AllEnhancers.table.regular.bed"),
                super=join(RESULTSDIR,"annotation","{peak_types}","{treatment}.{dupstatus}.{s_dist},{treatment}_peaks_AllEnhancers.table.super.bed"),
                regular_great=join(RESULTSDIR,"annotation","{peak_types}","{treatment}.{dupstatus}.{s_dist},{treatment}_peaks_AllEnhancers.table.super.GREAT.bed"),
                super_great=join(RESULTSDIR,"annotation","{peak_types}","{treatment}.{dupstatus}.{s_dist},{treatment}_peaks_peaks_AllEnhancers.table.regular.GREAT.bed"),
                regular_summit=join(RESULTSDIR,"annotation","{peak_types}","{treatment}.{dupstatus}.{s_dist},{treatment}_peaks_AllEnhancers.table.regular.summits.bed"),
                super_summit=join(RESULTSDIR,"annotation","{peak_types}","{treatment}.{dupstatus}.{s_dist},{treatment}_peaks_AllEnhancers.table.super.summits.bed"),
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
                prefix=join(RESULTSDIR,"annotation","{peak_types}","{treatment}.{dupstatus}.{s_dist}"),
                rose_script=join(SCRIPTSDIR,"roseTable2Bed.sh")
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
                {params.rose_script} {output.all} {output.regular} {output.super}
                
                # cut rose output files
                cut -f1-3 {output.regular} > {output.regular_great}
                cut -f1-3 {output.super} > {output.super_great}

                # summit files
                bedtools intersect -wa -a {input.bed} -b {output.regular} > {output.regular_summit}
                bedtools intersect -wa -a {input.bed} -b {output.super} > {output.super_summit}
                """

rule findMotif_all:
        input:
                bed=get_bed_all,
        output:
                known_html=join(RESULTSDIR,"annotation","{peak_types}","{treatment}.{dupstatus}.motifs","knownResults.html"),
        threads: getthreads("findMotif_all")
        envmodules:
                TOOLS["homer"],
        params:
                genome = config["genome"],
                outDir=join(RESULTSDIR,"annotation","{peak_types}","{treatment}.{dupstatus}.motifs"),
                motif_size = config["motif_size"],
                preparsedDir = config["preparsedDir"],
        shell:
                """
                findMotifsGenome.pl {input.bed} {params.genome} {params.outDir} -size {params.motif_size} -p {threads} -preparsedDir {params.preparsedDir}
                """