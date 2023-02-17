def get_bed_macs2(wildcards):
        if wildcards.macs_types == "narrowPeak":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"macs2", wildcards.replicate, wildcards.replicate + "." + wildcards.dupstatus + "_peaks.narrowPeak")
        if wildcards.macs_types =="broadPeak":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"macs2", wildcards.replicate, wildcards.replicate + "." + wildcards.dupstatus + "_peaks.broadPeak")
        return bed

rule peakAnnotation_macs2:
        """
        Developed from code: https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk
        """
        input:
                get_bed_macs2,
        output:
                annotation=join(RESULTSDIR,"annotation","{qthresholds}","{macs_types}","{replicate}","{replicate}.{dupstatus}.annotation.txt"),
                annotation_summary=join(RESULTSDIR,"annotation","{qthresholds}","{macs_types}","{replicate}","{replicate}.{dupstatus}.annotation.summary")
        envmodules:
                TOOLS["homer"],
        params:
                genome = config["genome"],
                fa=config["reference"][config["genome"]]["fa"],
                gtf = config["reference"][config["genome"]]["gtf"],
        shell:
                """  
                if [[ {params.genome} == "hs1" ]]; then
                        annotatePeaks.pl {input} {params.fa} -annStats {output.annotation_summary} -gtf {params.gtf} > {output.annotation}
                else
                        annotatePeaks.pl {input} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
                fi
                """

def get_bed_s_and_g(wildcards):
        # SEACR OPTIONS
        if wildcards.s_and_g_types =="norm.stringent.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr", wildcards.t_and_c, wildcards.t_and_c + "." + wildcards.dupstatus + ".norm.stringent.bed")
        if wildcards.s_and_g_types =="norm.relaxed.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr", wildcards.t_and_c, wildcards.t_and_c + "." + wildcards.dupstatus + ".norm.relaxed.bed")
        if wildcards.s_and_g_types =="non.stringent.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr", wildcards.t_and_c, wildcards.t_and_c + "." + wildcards.dupstatus + ".non.stringent.bed")
        if wildcards.s_and_g_types =="non.relaxed.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr", wildcards.t_and_c, wildcards.t_and_c + "." + wildcards.dupstatus + ".non.relaxed.bed")

        #GOPEAKS OPTIONS
        if wildcards.s_and_g_types =="narrowGo_peaks.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"gopeaks", wildcards.t_and_c + "." + wildcards.dupstatus + ".narrowGo_peaks.bed")
        if wildcards.s_and_g_types =="broadGo_peaks.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"gopeaks", wildcards.t_and_c + "." + wildcards.dupstatus + ".broadGo_peaks.bed")
        return bed

rule peakAnnotation_s_and_g:
        """
        Developed from code: https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk
        """
        input:
                get_bed_s_and_g,
        output:
                annotation=join(RESULTSDIR,"annotation","{qthresholds}","{s_and_g_types}","{t_and_c}","{t_and_c}.{dupstatus}.annotation.txt"),
                annotation_summary=join(RESULTSDIR,"annotation","{qthresholds}","{s_and_g_types}","{t_and_c}","{t_and_c}.{dupstatus}.annotation.summary"),
        envmodules:
                TOOLS["homer"],
        params:
                genome = config["genome"],
                fa=config["reference"][config["genome"]]["fa"],
                gtf = config["reference"][config["genome"]]["gtf"],
        shell:
                """  
                if [[ {params.genome} == "hs1" ]]; then
                        annotatePeaks.pl {input} {params.fa} -annStats {output.annotation_summary} -gtf {params.gtf} > {output.annotation}
                else
                        annotatePeaks.pl {input} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
                fi
                """

def get_cntrl_bam(wildcards):
        cntrl_sample=TREAT_to_CONTRL_DICT[wildcards.treatment]
        cntrl_file=join(RESULTSDIR, "bam", cntrl_sample + "." + wildcards.dupstatus + ".bam"),
        return cntrl_file

def get_bed_all(wildcards):
        cntrl=TREAT_to_CONTRL_DICT[wildcards.treatment]
        #MACS2 OPTIONS
        if wildcards.peak_types == "narrowPeak":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"macs2", wildcards.treatment, wildcards.treatment + "." + wildcards.dupstatus + "_peaks.narrowPeak")
        if wildcards.peak_types =="broadPeak":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"macs2", wildcards.treatment, wildcards.treatment + "." + wildcards.dupstatus + "_peaks.broadPeak")
        
        #SEACR OPTIONS
        if wildcards.peak_types =="norm.stringent.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr", wildcards.treatment + "_vs_" + cntrl, wildcards.treatment + "_vs_" + cntrl + "." + wildcards.dupstatus + ".norm.stringent.bed")
        if wildcards.peak_types =="norm.relaxed.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr", wildcards.treatment + "_vs_" + cntrl, wildcards.treatment + "_vs_" + cntrl + "." + wildcards.dupstatus + ".norm.relaxed.bed")
        if wildcards.peak_types =="non.stringent.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr", wildcards.treatment + "_vs_" + cntrl, wildcards.treatment + "_vs_" + cntrl + "." + wildcards.dupstatus + ".non.stringent.bed")
        if wildcards.peak_types =="non.relaxed.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr", wildcards.treatment + "_vs_" + cntrl, wildcards.treatment + "_vs_" + cntrl + "." + wildcards.dupstatus + ".non.relaxed.bed")
        
        #GOPEAKS OPTIONS
        if wildcards.peak_types =="narrowGo_peaks.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"gopeaks", wildcards.treatment + "_vs_" + cntrl + "." + wildcards.dupstatus + ".narrowGo_peaks.bed")
        if wildcards.peak_types =="broadGo_peaks.bed":
                bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"gopeaks", wildcards.treatment + "_vs_" + cntrl + "." + wildcards.dupstatus + ".broadGo_peaks.bed")
        return bed

rule rose:
        """
        Developed from code: https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk

        # SEACR bed file output format
        <1>     <2>     <3>     <4>             <5>             <6>
        <chr>   <start> <end>   <total signal>  <max signal>	<max signal region>

        # GOPEAKS bed file output format
        <1>     <2>     <3>
        <chr>   <start> <end>
        chr1	29107	29364	42861.3	222.849	chr1:29151-29291


        # MACS2 bed file output format
        https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md
        <1>     <2>     <3>     <4>     <5>                             <6>     <7>             <8>             <9>             <10>
        <chr>   <start> <end>   <name>  <integer score for display>     <empty> <fold-change> <-log10pvalue>    <-log10qvalue>  <relative summit position to peak start>
        ## integer score for display: It's calculated as int(-10*log10pvalue) or int(-10*log10qvalue) depending on whether -p (pvalue) or -q (qvalue) is used as score cutoff. 
        ### Please note that currently this value might be out of the [0-1000] range defined in UCSC ENCODE narrowPeak format. You can let the value saturated at 1000 (i.e. p/q-value = 10^-100) 
        ### by using the following 1-liner awk: awk -v OFS="\t" '{$5=$5>1000?1000:$5} {print}' NAME_peaks.narrowPeak
        ## broadPeak does not have 10th column
        ### Since in the broad peak calling mode, the peak summit won't be called, the values in the 5th, and 7-9th columns are the mean value across all positions in the peak region

        # rose input (differs from documentation)
        column 1: chromosome (chr#)
        column 2: unique ID for each constituent enhancer region
        column 3: start of constituent
        column 4: end of constituent
        Column 5: ignored
        column 6: strand (+,-,.)

        """
        input:
                bed=get_bed_all,
                bam=join(RESULTSDIR,"bam", "{treatment}.{dupstatus}.bam"),
                cntrl_bam=get_cntrl_bam,
        output:
                no_tss_bed=join(RESULTSDIR,"annotation","{qthresholds}","{peak_types}","{treatment}","{treatment}.{dupstatus}.no_TSS_{s_dist}.bed"),
                regular_summit=join(RESULTSDIR,"annotation","{qthresholds}","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllStitched.table.regular.summits.bed"),
                super_summit=join(RESULTSDIR,"annotation","{qthresholds}","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllStitched.table.super.summits.bed"),
        envmodules:
                TOOLS["bedtools"],
                TOOLS["rose"],
                TOOLS["python37"],
                TOOLS["samtools"],
        threads: getthreads("rose")
        params:
                genome = config["genome"],
                regions=config["reference"][config["genome"]]["regions"],
                tss_bed = config["reference"][config["genome"]]["tss_bed"],
                stitch_distance = config["stitch_distance"],
                tss_distance=config["tss_distance"],
                peak_type="{peak_types}",
                workdir=join(WORKDIR),
                refseq=config["reference"][config["genome"]]["rose"],
                sampleID="{treatment}",
                prefix=join(RESULTSDIR,"annotation","{qthresholds}","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}"),
                all=join(RESULTSDIR,"annotation","{qthresholds}","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllStitched.table.txt"),
                regular=join(RESULTSDIR,"annotation","{qthresholds}","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}AllEnhancers.table.regular.bed"),
                super=join(RESULTSDIR,"annotation","{qthresholds}","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllStitched.table.super.bed"),
                regular_great=join(RESULTSDIR,"annotation","{qthresholds}","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllStitched.table.super.GREAT.bed"),
                super_great=join(RESULTSDIR,"annotation","{qthresholds}","{peak_types}","{treatment}","{treatment}.{dupstatus}.{s_dist}","{treatment}_AllStitched.table.regular.GREAT.bed"),
        shell:
                """
                # set tmp
                set -exo pipefail
                if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then 
                        TMPDIR="/lscratch/$SLURM_JOB_ID"
                else
                        dirname=$(basename $(mktemp))
                        TMPDIR="/dev/shm/$dirname"
                        mkdir -p $TMPDIR
                fi
                
                # set ROSE specific paths
                PATHTO=/usr/local/apps/ROSE/1.3.1
                PYTHONPATH=/usr/local/apps/ROSE/1.3.1/src/lib
                export PYTHONPATH
                export PATH=$PATH:$PATHTO/bin
                
                # remove NC from bams and beds
                echo "Cleaning"
                samtools view -b {input.bam} {params.regions} > $TMPDIR/subset.bam
                samtools index $TMPDIR/subset.bam
                grep -v "NC_" {input.bed} > $TMPDIR/subset.bed
                      
                # prep for ROSE
                # macs2 output is prepared for ROSE formatting
                # seacr and gopeaks must be edited to correct for formatting
                ## correct GOPEAKS
                ### original: <chr>   <start> <end>
                ### output: <chr>   <start> <end> <$sampleid_uniquenumber> <0> <.>
                echo "prep Rose"
                if [[ {params.peak_type} == "narrowGo_peaks.bed" ]] || [[ {params.peak_type} == "broadGo_peaks.bed" ]]; then
                        cp $TMPDIR/subset.bed $TMPDIR/save.bed
                        nl --number-format=rz --number-width=3 $TMPDIR/subset.bed | awk -v sample_id="{params.sampleID}_" \'{{print sample_id$1"\\t0\\t."}}\' > $TMPDIR/col.txt
                        paste -d "\t" $TMPDIR/save.bed $TMPDIR/col.txt > $TMPDIR/subset.bed
                fi
                ## correct SEACR
                ### original: <chr>   <start> <end>   <total signal>  <max signal>	<max signal region>
                ### output: <chr>   <start> <end> <$sampleid_uniquenumber> <total signal> <.>
                if [[ {params.peak_type} == "norm.stringent.bed" ]] || [[ {params.peak_type} == "norm.relaxed.bed" ]] || [[ {params.peak_type} == "non.relaxed.bed" ]] || [[ {params.peak_type} == "non.relaxed.bed" ]]; then
                        cp $TMPDIR/subset.bed $TMPDIR/save.bed
                        awk -v sample_id="{params.sampleID}_" \'{{print $1"\\t"$2"\\t"$3"\\t"sample_id$1"\\t"$4"\\t."}}\' $TMPDIR/subset.bed > $TMPDIR/col.txt
                        paste -d "\t" $TMPDIR/save.bed $TMPDIR/col.txt > $TMPDIR/subset.bed
                fi

                # bedtools
                echo "intersects"
                bedtools intersect -a $TMPDIR/subset.bed -b {params.tss_bed} -v > $TMPDIR/tmp.bed
                bedtools merge -i $TMPDIR/tmp.bed -d {params.stitch_distance} -c 4,5,6 -o distinct,sum,distinct > {output.no_tss_bed}

                # if there are less than 5 peaks, annotation will fail
                num_of_peaks=`cat {output.no_tss_bed} | wc -l`
                if [[ $num_of_peaks -gt 5 ]]; then
                        # rose
                        # developed from outdated verison of rose: https://github.com/younglab/ROSE
                        echo "rose"
                        cd {params.workdir}
                        
                        if [[ {params.genome} == "hs1" ]]; then
                                ROSE_main.py \
                                        -i {output.no_tss_bed} \
                                        --custom={params.refseq} \
                                        -r $TMPDIR/subset.bam {input.cntrl_bam} \
                                        -t {params.tss_distance} \
                                        -s {params.stitch_distance} \
                                        -o {params.prefix}                        
                        else
                                ROSE_main.py \
                                        -i {output.no_tss_bed} \
                                        -g {params.genome} \
                                        -r $TMPDIR/subset.bam {input.cntrl_bam} \
                                        -t {params.tss_distance} \
                                        -s {params.stitch_distance} \
                                        -o {params.prefix}
                        fi
                        
                        # rose to bed file
                        echo "convert bed"
                        # developed from https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/scripts/roseTable2Bed.sh
                        grep -v "^[#|REGION]" {params.all} | awk -v OFS="\\t" -F"\\t" \'$NF==0 {{for(i=2; i<=NF; i++){{printf $i; printf (i<NF?"\\t":"\\n")}}}}\' > $TMPDIR/regular
                        bedtools sort -i $TMPDIR/regular > {params.regular}
                        
                        grep -v "^[#|REGION]" {params.all} | awk -v OFS="\\t" -F"\\t" \'$NF==1 {{for(i=2; i<=NF; i++){{printf $i; printf (i<NF?"\\t":"\\n")}}}}\' > $TMPDIR/super
                        bedtools sort -i $TMPDIR/super > {params.super}
                        
                        # cut rose output files, create summits
                        echo "cut and summit"
                        cut -f1-3 {params.regular} > {params.regular_great}
                        cut -f1-3 {params.super} > {params.super_great}
                        bedtools intersect -wa -a {input.bed} -b {params.regular} > {output.regular_summit}
                        bedtools intersect -wa -a {input.bed} -b {params.super} > {output.super_summit}
                else
                        echo "Less than 5 usable peaks detected ($num_of_peaks)" > {output.regular_summit}
                        echo "Less than 5 usable peaks detected ($num_of_peaks)" > {output.super_summit}
                fi
                """


rule findMotif:
        """
        Developed from code: https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk

        Notes on using alternative genomes now in config
        - http://homer.ucsd.edu/homer/introduction/update.html
        - Config is located on Biowulf here: /usr/local/apps/homer/4.11.1/.//config.txt

        """
        input:
                bed=get_bed_all,
        output:
                known_html=join(RESULTSDIR,"annotation","{qthresholds}","{peak_types}","{treatment}","{treatment}.{dupstatus}.motifs","knownResults.html"),
        threads: getthreads("findMotif")
        envmodules:
                TOOLS["homer"],
        params:
                genome = config["genome"],
                outDir=join(RESULTSDIR,"annotation","{qthresholds}","{peak_types}","{treatment}","{treatment}.{dupstatus}.motifs"),
                motif_size = config["motif_size"],
                preparsedDir = config["preparsedDir"],
                fa=config["reference"][config["genome"]]["fa"]
        shell:
                """
                # hs1 is not part of HOMER's config genome db. Must add it as a separate param
                if [[ {params.genome} == "hs1" ]]; then
                        findMotifsGenome.pl {input.bed} {params.fa} {params.outDir} -size {params.motif_size} -p {threads} -preparsedDir {params.preparsedDir}
                else
                        findMotifsGenome.pl {input.bed} {params.genome} {params.outDir} -size {params.motif_size} -p {threads} -preparsedDir {params.preparsedDir}
                fi
                """