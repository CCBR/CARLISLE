def get_peakAnnotation_input(wildcards):
    files=[]
    if "narrowPeak" in config["peaktype"]:
        n=expand(join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.narrowPeak"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        files.extend(n)
    if "broadPeak" in config["peaktype"]:
        b=expand(join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.broadPeak"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        files.extend(b)
    if "norm.stringent.bed" in config["peaktype"]:
        s=expand([join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.norm.stringent.bed")],zip,treatment=TREATMENTS,control=CONTROLS,dupstatus=DUPSTATUS),
        files.extend(s)
    if "norm.relaxed.bed" in config["peaktype"]:
        r=expand([join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.{dupstatus}.norm.relaxed.bed")],zip,treatment=TREATMENTS,control=CONTROLS,dupstatus=DUPSTATUS),
        files.extend(r)
    if "narrowGo_peaks.bed" in config["peaktype"]:
        n=expand([join(RESULTSDIR,"peaks","gopeaks","{treatment}_vs_{control}.dedup.narrowGo_peaks.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        files.extend(n)
    if "broadGo_peaks.bed" in config["peaktype"]:
        b=expand([join(RESULTSDIR,"peaks","gopeaks","{treatment}_vs_{control}.dedup.broadGo_peaks.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        files.extend(b)
    return files

print("Im here")
rule peakAnnotation:
    input:
        get_peakAnnotation_input,
    output:
        annotation=join(RESULTSDIR,"annotation","{pt}","{replicate}.{dupstatus}_peaks.annotation.txt"),pt=PT,replicate=REPLICATES,dupstatus=DUPSTATUS),
        annotation=join(RESULTSDIR,"annotation","{pt}","{replicate}.{dupstatus}_peaks.annotation.summary"),pt=PT,replicate=REPLICATES,dupstatus=DUPSTATUS)
    envmodules:
        TOOLS["homer"]
    params:
        genome = config["genome"],
        pl_script = join(SCRIPTSDIR,"_annotatePeaks.pl"),  
    shell:
        """  
        {params.pl_script} {input} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
        """

rule rose:
    input:
            bed=lambda wildcards: wildcards.sample + "/MACS_Out_" + wildcards.cutoff + "/" + wildcards.sample + "_peaks." + samples[wildcards.sample]["PeakCalling"] + "Peak.nobl.bed",
            #bed="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.no_TSS.bed",
            bam="{sample}/{sample}.bam"
    output:
            rose_out="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.txt",
            se_table="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.super.bed",
            re_table="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.regular.bed",
            se_great="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.super.GREAT.bed",
            re_great="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.regular.GREAT.bed",
            #no_tss="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.no_TSS.bed",
            #no_tss_stitched="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.no_TSS_{stitch_distance}.bed",
    version:
            config["version"]["rose"],            
    benchmark:
            "{sample}/benchmark/rose.{sample}.{cutoff}.{stitch_distance}benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_rose"],
            pipeline_home=config["pipeline_home"],
            peak_type = lambda wildcards: samples[wildcards.sample]["PeakCalling"],
            tss_distance=config["rose"]["tss_distance"],
            tss_bed=lambda wildcards: config["pipeline_home"] + "/" + config[samples[wildcards.sample]["Genome"]]["tss_bed"],
            genome = lambda wildcards: samples[wildcards.sample]["Genome"].upper(),
            annotation = lambda wildcards: config[samples[wildcards.sample]["Genome"]]["rose"],
            input_control = lambda wildcards: " -c " + config["work_dir"] + "/" + samples[wildcards.sample]["Matched normal"] + "/" + samples[wildcards.sample]["Matched normal"] + ".bam" if samples[wildcards.sample]["Matched normal"] != "." and samples[wildcards.sample]["Matched normal"] != "" else "",
            rulename = "rose",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            version_R=config["version_common"]["R"],
            version_python=config["version_common"]["python2"],
            version_bedtools=config["version_common"]["bedtools"],
    shell:
            """
            module load bedtools/{params.version_bedtools}
            bedtools intersect -a {input.bed} -b {params.tss_bed} -v > {wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS.bed
            bedtools merge -i {wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS.bed -d {wildcards.stitch_distance} -c 4,5,6 -o distinct,sum,distinct > {wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS_{wildcards.stitch_distance}.bed
            module load rose/{version}
            #module load bamliquidator
            module load python/{params.version_python}
            module load R/{params.version_R}
            #cd /usr/local/apps/bamliquidator/pipeline
            cd $(dirname `which ROSE_main.py`)
            #cd {params.pipeline_home}/apps/rose
            #export PATH=$PATH:{params.pipeline_home}/apps/rose
            ./ROSE_main.py -i {params.work_dir}/{wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS_{wildcards.stitch_distance}.bed -g {params.genome} -r {params.work_dir}/{input.bam} {params.input_control} -t {params.tss_distance} -s {wildcards.stitch_distance} -o {params.work_dir}/{wildcards.sample}/MACS_Out_{wildcards.cutoff}/ROSE_out_{wildcards.stitch_distance}
            {params.pipeline_home}/scripts/roseTable2Bed.sh {params.work_dir}/{output.rose_out} {params.work_dir}/{output.re_table} {params.work_dir}/{output.se_table}
            cd {params.work_dir}
            cut -f1-3 {output.se_table} > {output.se_great}
            cut -f1-3 {output.re_table} > {output.re_great}
            """

rule findMotif:
    input:
            "{sample}/{folders}/{sample}.{peak_type}_summits.bed",
    output:
            "{sample}/{folders}/motif_{peak_type}/knownResults.html",
    wildcard_constraints:
            folders = "MACS.+"
    version:
            config["version"]["homer"]
    benchmark:
            "{sample}/benchmark/findMotif.{sample}.{folders}.{peak_type}.benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_motif"],
            pipeline_home=config["pipeline_home"],
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            rulename = "findMotif",
            motif_size = str(config["homer"]["motif_size"]),
            log_dir = lambda wildcards: wildcards.sample + '/log',
    shell:
            """  
            module load homer/{version}
            out_dir=$(dirname {output})
            if [ -s {input} ];then
                findMotifsGenome.pl {input} {params.genome} ${{out_dir}} -size {params.motif_size} -p ${{THREADS}} -preparsedDir {params.pipeline_home}/ref/preparsedDir
            else
                touch {output}
            fi
            """

rule prepareRoseSummit:
    input:
            summit=lambda wildcards: wildcards.sample + "/MACS_Out_" + wildcards.cutoff + "/" + wildcards.sample + "." + samples[wildcards.sample]["PeakCalling"] + "_summits.bed",
            se_table="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.super.bed",
            re_table="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.regular.bed",
    output:
            se_summit="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}.super_summits.bed",
            re_summit="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}.regular_summits.bed",
    version:
            config["version_common"]["bedtools"]
    params:
            work_dir = config["work_dir"],
            pipeline_home=config["pipeline_home"],
            rulename = "prepareRoseSummit",
            log_dir = lambda wildcards: wildcards.sample + '/log',
    shell:
            """  
            module load bedtools/{version}
            bedtools intersect -wa -a {input.summit} -b {input.se_table} > {output.se_summit}
            bedtools intersect -wa -a {input.summit} -b {input.re_table} > {output.re_summit}
            """

    input:
            bed=lambda wildcards: wildcards.sample + "/MACS_Out_" + wildcards.cutoff + "/" + wildcards.sample + "_peaks." + samples[wildcards.sample]["PeakCalling"] + "Peak.nobl.bed",
            #bed="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.no_TSS.bed",
            bam="{sample}/{sample}.bam"
    output:
            rose_out="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.txt",
            se_table="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.super.bed",
            re_table="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.regular.bed",
            se_great="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.super.GREAT.bed",
            re_great="{sample}/MACS_Out_{cutoff}/ROSE_out_{stitch_distance}/{sample}_peaks_AllEnhancers.table.regular.GREAT.bed",
            #no_tss="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.no_TSS.bed",
            #no_tss_stitched="{sample}/MACS_Out_{cutoff}/{sample}_peaks.{peak_type}Peak.nobl.no_TSS_{stitch_distance}.bed",
    version:
            config["version"]["rose"],            
    benchmark:
            "{sample}/benchmark/rose.{sample}.{cutoff}.{stitch_distance}benchmark.txt"
    params:
            work_dir = config["work_dir"],
            batch    = config["cluster"]["job_rose"],
            pipeline_home=config["pipeline_home"],
            peak_type = lambda wildcards: samples[wildcards.sample]["PeakCalling"],
            tss_distance=config["rose"]["tss_distance"],
            tss_bed=lambda wildcards: config["pipeline_home"] + "/" + config[samples[wildcards.sample]["Genome"]]["tss_bed"],
            genome = lambda wildcards: samples[wildcards.sample]["Genome"].upper(),
            annotation = lambda wildcards: config[samples[wildcards.sample]["Genome"]]["rose"],
            input_control = lambda wildcards: " -c " + config["work_dir"] + "/" + samples[wildcards.sample]["Matched normal"] + "/" + samples[wildcards.sample]["Matched normal"] + ".bam" if samples[wildcards.sample]["Matched normal"] != "." and samples[wildcards.sample]["Matched normal"] != "" else "",
            rulename = "rose",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            version_R=config["version_common"]["R"],
            version_python=config["version_common"]["python2"],
            version_bedtools=config["version_common"]["bedtools"],
    shell:
            """
            module load bedtools/{params.version_bedtools}
            bedtools intersect -a {input.bed} -b {params.tss_bed} -v > {wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS.bed
            bedtools merge -i {wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS.bed -d {wildcards.stitch_distance} -c 4,5,6 -o distinct,sum,distinct > {wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS_{wildcards.stitch_distance}.bed
            module load rose/{version}
            #module load bamliquidator
            module load python/{params.version_python}
            module load R/{params.version_R}
            #cd /usr/local/apps/bamliquidator/pipeline
            cd $(dirname `which ROSE_main.py`)
            #cd {params.pipeline_home}/apps/rose
            #export PATH=$PATH:{params.pipeline_home}/apps/rose
            ./ROSE_main.py -i {params.work_dir}/{wildcards.sample}/MACS_Out_{wildcards.cutoff}/{wildcards.sample}_peaks.{params.peak_type}Peak.nobl.no_TSS_{wildcards.stitch_distance}.bed -g {params.genome} -r {params.work_dir}/{input.bam} {params.input_control} -t {params.tss_distance} -s {wildcards.stitch_distance} -o {params.work_dir}/{wildcards.sample}/MACS_Out_{wildcards.cutoff}/ROSE_out_{wildcards.stitch_distance}
            {params.pipeline_home}/scripts/roseTable2Bed.sh {params.work_dir}/{output.rose_out} {params.work_dir}/{output.re_table} {params.work_dir}/{output.se_table}
            cd {params.work_dir}
            cut -f1-3 {output.se_table} > {output.se_great}
            cut -f1-3 {output.re_table} > {output.re_great}
            """
