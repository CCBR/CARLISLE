def get_peakAnnotation_macs2(wildcards):
        if wildcards.macs_types == "narrowPeak":
                bed=join(RESULTSDIR,"peaks","macs2", wildcards.replicate, wildcards.replicate + "." + wildcards.dupstatus + "_peaks.narrowPeak")
        if wildcards.macs_types =="broadPeak":
                bed=join(RESULTSDIR,"peaks","macs2", wildcards.replicate, wildcards.replicate + "." + wildcards.dupstatus + "_peaks.broadPeak")
        return bed

rule peakAnnotation_macs2:
        input:
                get_peakAnnotation_macs2,
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

def get_peakAnnotation_s_and_g(wildcards):
        if wildcards.s_and_g_types =="norm.stringent.bed":
                bed=join(RESULTSDIR,"peaks","seacr", wildcards.t_and_c, wildcards.t_and_c + "." + wildcards.dupstatus + ".norm.stringent.bed")
        if wildcards.s_and_g_types =="norm.relaxed.bed":
                bed=join(RESULTSDIR,"peaks","seacr", wildcards.t_and_c, wildcards.t_and_c + "." + wildcards.dupstatus + ".norm.relaxed.bed")
        if wildcards.s_and_g_types =="narrowGo_peaks.bed":
                bed=join(RESULTSDIR,"peaks","gopeaks", wildcards.t_and_c + "." + wildcards.dupstatus + ".narrowGo_peaks.bed")
        if wildcards.s_and_g_types =="broadGo_peaks.bed":
                bed=join(RESULTSDIR,"peaks","gopeaks", wildcards.t_and_c + "." + wildcards.dupstatus + ".broadGo_peaks.bed")
        return bed

rule peakAnnotation_s_and_g:
    input:
        get_peakAnnotation_s_and_g,
    output:
        annotation=join(RESULTSDIR,"annotation","{s_and_g_types}","{t_and_c}.{dupstatus}_peaks.annotation.txt"),
        annotation_summary=join(RESULTSDIR,"annotation","{s_and_g_types}","{t_and_c}.{dupstatus}_peaks.annotation.summary")
        f=expand(join(RESULTSDIR,"annotation","{s_and_g_types}","{t_and_c}.{dupstatus}_peaks.annotation.txt"),s_and_g_types=S_AND_G_TYPES,t_and_c=TREATMENT_CONTROL_LIST,dupstatus=DUPSTATUS),
    envmodules:
        TOOLS["homer"],
    params:
        genome = config["genome"],
    shell:
        """  
        annotatePeaks.pl {input} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
        """