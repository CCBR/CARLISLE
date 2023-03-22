def get_peak_file(wildcards):
    # MACS2 OPTIONS
    if wildcards.peak_caller_type == "macs2_narrow":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"macs2","peak_output", wildcards.treatment_control_list + "." + wildcards.dupstatus + ".narrow.peaks.bed")
    if wildcards.peak_caller_type == "macs2_broad":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"macs2","peak_output",wildcards.treatment_control_list + "." + wildcards.dupstatus + ".broad.peaks.bed")

    # SEACR OPTIONS
    if wildcards.peak_caller_type =="seacr_norm_stringent":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr","peak_output",wildcards.treatment_control_list + "." + wildcards.dupstatus + ".norm_stringent.peaks.bed")
    if wildcards.peak_caller_type =="seacr_norm_relaxed":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr","peak_output",wildcards.treatment_control_list + "." + wildcards.dupstatus + ".norm_relaxed.peaks.bed")
    if wildcards.peak_caller_type =="seacr_non_stringent":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr","peak_output",wildcards.treatment_control_list + "." + wildcards.dupstatus + ".non_stringent.peaks.bed")
    if wildcards.peak_caller_type =="seacr_non_relaxed":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr","peak_output",wildcards.treatment_control_list + "." + wildcards.dupstatus + ".non_relaxed.peaks.bed")

    #GOPEAKS OPTIONS
    if wildcards.peak_caller_type =="gopeaks_narrow":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"gopeaks","peak_output",wildcards.treatment_control_list + "." + wildcards.dupstatus + ".narrow.peaks.bed")
    if wildcards.peak_caller_type =="gopeaks_broad":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"gopeaks","peak_output",wildcards.treatment_control_list + "." + wildcards.dupstatus + ".broad.peaks.bed")
    return bed

# localrules: create_contrast_peakcaller_files
rule findMotif:
    """
    Developed from code: https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk

    Notes on using alternative genomes now in config
    - http://homer.ucsd.edu/homer/introduction/update.html
    - Config is located on Biowulf here: /usr/local/apps/homer/4.11.1/.//config.txt
    """
    input:
        peak_file=get_peak_file,
    output:
        annotation=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{treatment_control_list}.{dupstatus}.{peak_caller_type}.annotation.txt"),
        annotation_summary=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{treatment_control_list}.{dupstatus}.{peak_caller_type}.annotation.summary"),
        known_html=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{treatment_control_list}.{dupstatus}.{peak_caller_type}.motifs","knownResults.html"),
    threads: getthreads("findMotif")
    envmodules:
        TOOLS["homer"],
    params:
        genome = config["genome"],
        fa=config["reference"][config["genome"]]["fa"],
        gtf = config["reference"][config["genome"]]["gtf"],
        outDir=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{treatment_control_list}.{dupstatus}.{peak_caller_type}.motifs"),
        motif_size = config["motif_size"],
        preparsedDir = config["preparsedDir"],
    shell:
        """
        # run homer
        if [[ {params.genome} == "hs1" ]]; then
            annotatePeaks.pl {input.peak_file} {params.fa} -annStats {output.annotation_summary} -gtf {params.gtf} > {output.annotation}
        else
            annotatePeaks.pl {input.peak_file} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
        fi
        
        # hs1 is not part of HOMER's config genome db. Must add it as a separate param
        if [[ {params.genome} == "hs1" ]]; then
            findMotifsGenome.pl {input.peak_file} {params.fa} {params.outDir} -size {params.motif_size} -p {threads} -preparsedDir {params.preparsedDir}
        else
            findMotifsGenome.pl {input.peak_file} {params.genome} {params.outDir} -size {params.motif_size} -p {threads} -preparsedDir {params.preparsedDir}
        fi
        """

rule rose:
    """
    Developed from code: 
    https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk
    outdated verison of rose: https://github.com/younglab/ROSE
            
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
        peak_file=get_peak_file,
        bam = expand(join(RESULTSDIR,"bam","{replicate}.{dupstatus}.bam"),replicate=REPLICATES,dupstatus=DUPSTATUS),
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
        refseq=config["reference"][config["genome"]]["rose"],
        stitch_distance = config["stitch_distance"],
        tss_distance=config["tss_distance"],
        tc_file="{treatment_control_list}",
        peak_caller_type="{peak_caller_type}",
        dupstatus = "{dupstatus}",
        bam_path=join(RESULTSDIR,"bam"),
        workdir=join(WORKDIR),
        file_base=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}"),
        control_flag = config["macs2_control"],
        mapped_gff_dir=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","mappedGFF"),
        gff_dir=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","gff"),
    output:
        no_tss_bed=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{treatment_control_list}.{dupstatus}.{peak_caller_type}.no_TSS_{s_dist}.bed"),
        all=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","{treatment_control_list}_AllStitched.table.txt"),
        regular=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","{treatment_control_list}_AllEnhancers.table.regular.bed"),
        super=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","{treatment_control_list}_AllStitched.table.super.bed"),
        regular_great=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","{treatment_control_list}_AllStitched.table.super.GREAT.bed"),
        super_great=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","{treatment_control_list}_AllStitched.table.regular.GREAT.bed"),
        regular_summit=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","{treatment_control_list}_AllStitched.table.regular.summits.bed"),
        super_summit=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","{treatment_control_list}_AllStitched.table.super.summits.bed"),
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
                
        # pull treatment and control ids
        treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
        control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`

        # set bam file
        treat_bam={params.bam_path}/${{treatment}}.{params.dupstatus}.bam
        cntrl_bam={params.bam_path}/${{control}}.{params.dupstatus}.bam

        # remove NC from bams and beds
        echo "## Cleaning"
        samtools view -b ${{treat_bam}} {params.regions} > $TMPDIR/subset.bam
        samtools index $TMPDIR/subset.bam
        grep -v "NC_" {input.peak_file} > $TMPDIR/subset.bed
                      
        # prep for ROSE
        # macs2 output is prepared for ROSE formatting
        # seacr and gopeaks must be edited to correct for formatting
        
        ## correct GOPEAKS
        ### original: <chr>   <start> <end>
        ### output: <chr>   <start> <end> <$sampleid_uniquenumber> <0> <.>
        echo "## Prep Rose"
        if [[ {params.peak_caller_type} == "gopeaks_narrow" ]] || [[ {params.peak_caller_type} == "gopeaks_broad" ]]; then
            echo "#### Fixing GoPeaks"
            cp $TMPDIR/subset.bed $TMPDIR/save.bed
            nl --number-format=rz --number-width=3 $TMPDIR/subset.bed | awk -v sample_id="${{treatment}}_" \'{{print sample_id$1"\\t0\\t."}}\' > $TMPDIR/col.txt
            paste -d "\t" $TMPDIR/save.bed $TMPDIR/col.txt > $TMPDIR/subset.bed
        fi
        
        ## correct SEACR
        ### original: <chr>   <start> <end>   <total signal>  <max signal>	<max signal region>
        ### output: <chr>   <start> <end> <$sampleid_uniquenumber> <total signal> <.>
        if [[ {params.peak_caller_type} == "seacr_norm_stringent" ]] || [[ {params.peak_caller_type} == "seacr_norm_relaxed" ]] || [[ {params.peak_caller_type} == "seacr_non_relaxed" ]] || [[ {params.peak_caller_type} == "seacr_non_relaxed" ]]; then
            echo "#### Fixing SECAR"
            cp $TMPDIR/subset.bed $TMPDIR/save.bed
            awk -v sample_id="${{treatment}}_" \'{{print $1"\\t"$2"\\t"$3"\\t"sample_id$1"\\t"$4"\\t."}}\' $TMPDIR/subset.bed > $TMPDIR/col.txt
            paste -d "\t" $TMPDIR/save.bed $TMPDIR/col.txt > $TMPDIR/subset.bed
        fi

        # bedtools
        echo "## Intersecting"
        bedtools intersect -a $TMPDIR/subset.bed -b {params.tss_bed} -v > $TMPDIR/tmp.bed
        bedtools merge -i $TMPDIR/tmp.bed -d {params.stitch_distance} -c 4,5,6 -o distinct,sum,distinct > {output.no_tss_bed}

        # if there are less than 5 peaks, annotation will fail
        # if there are more, run ROSE
        num_of_peaks=`cat {output.no_tss_bed} | wc -l`
        if [[ ${{num_of_peaks}} -gt 5 ]]; then
            echo "## More than 5 usable peaks detected ${{num_of_peaks}} - Running rose"
            cd {params.workdir}
            
            # if macs2 control is off, there will be no macs2 control to annotate
            if [[ {params.control_flag} == "N" ]] & [[ {params.peak_caller_type} == "macs2_narrow" ]] || [[ {params.peak_caller_type} == "macs2_broad" ]]; then
                echo "#### No control was used"
                rose_files="$TMPDIR/subset.bam"  
            else
                echo "#### A control used ${{cntrl_bam}}"
                rose_files="$TMPDIR/subset.bam ${{cntrl_bam}}"
            fi

            if [[ {params.genome} == "hs1" ]]; then
                ROSE_main.py \
                    -i {output.no_tss_bed} \
                    --custom={params.refseq} \
                    -r ${{rose_files}} \
                    -t {params.tss_distance} \
                    -s {params.stitch_distance} \
                    -o {params.file_base}                        
            else
                ROSE_main.py \
                    -i {output.no_tss_bed} \
                    -g {params.genome} \
                    -r ${{rose_files}} \
                    -t {params.tss_distance} \
                    -s {params.stitch_distance} \
                    -o {params.file_base}
            fi
                        
            # rose to bed file
            echo "## Convert bed"
            # developed from https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/scripts/roseTable2Bed.sh
            grep -v "^[#|REGION]" {output.all} | awk -v OFS="\\t" -F"\\t" \'$NF==0 {{for(i=2; i<=NF; i++){{printf $i; printf (i<NF?"\\t":"\\n")}}}}\' > $TMPDIR/regular
            bedtools sort -i $TMPDIR/regular > {output.regular}
                        
            grep -v "^[#|REGION]" {output.all} | awk -v OFS="\\t" -F"\\t" \'$NF==1 {{for(i=2; i<=NF; i++){{printf $i; printf (i<NF?"\\t":"\\n")}}}}\' > $TMPDIR/super
            bedtools sort -i $TMPDIR/super > {output.super}
                        
            # cut rose output files, create summits
            echo "## Cut and summit"
            cut -f1-3 {output.regular} > {output.regular_great}
            cut -f1-3 {output.super} > {output.super_great}
            bedtools intersect -wa -a {input.peak_file} -b {output.regular} > {output.regular_summit}
            bedtools intersect -wa -a {input.peak_file} -b {output.super} > {output.super_summit}

            # cleanup
            echo "## Cleaning up"
            rm -r {params.mapped_gff_dir}
            rm -r {params.gff_dir}
        else
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})"
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.no_tss_bed}
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.all}
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.regular}
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.super}
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.regular_great}
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.super_great}
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.regular_summit}
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.super_summit}
        fi
    """

rule create_contrast_peakcaller_files:
    """
    Reads in all of the output from Rules create_contrast_data_files which match the same peaktype and merges them together
    """
    input:
        contrast_files=expand(join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.txt"),qthresholds=QTRESHOLDS, contrast_list=CONTRAST_LIST,dupstatus=DUPSTATUS,peak_caller_type=PEAKTYPE)
    params:
        qthresholds = "{qthresholds}",
        contrast_list = "{contrast_list}",
        dupstatus = "{dupstatus}",
        peak_caller = "{peak_caller}"
    output:
        peak_contrast_files=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","go_enrichment","{contrast_list}.{dupstatus}.txt")
    shell:
        """
        set -exo pipefail
        if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then 
            TMPDIR="/lscratch/$SLURM_JOB_ID"
        else
            TMPDIR="/dev/shm"
        fi

        # for each of the file, find matches to the peak_type
        if [[ -f $$TMPDIR/merge.txt ]]; then rm $TMPDIR/merge.txt; fi

        for f in {input.contrast_files}; do
            touch $TMPDIR/merge.txt

            # pull peaktype: macs2, seacr, gopeaks
            # /data/sevillas2/carlisle/v2.0/results/peaks/0.05/contrasts/53_H3K4me3_vs_HN6_H3K4me3.dedup/53_H3K4me3_vs_HN6_H3K4me3.dedup.seacr_norm_stringent.txt
            qthresholds=`echo $f | cut -f10 -d"/"`
            contrast_list=`echo $f | cut -f12 -d"/" | cut -f1 -d"."`
            dupstatus=`echo $f | cut -f12 -d"/" | cut -f2 -d"."`
            peak_caller=`echo $f | cut -f13 -d"/" | cut -f3 -d"." | cut -f1 -d"_"`

            if [[ $qthresholds == {params.qthresholds} ]] && [[ $contrast_list == {params.contrast_list} ]] && [[ $dupstatus == {params.dupstatus} ]] && [[ $peak_caller == {params.peak_caller} ]]; then
                cat $f >> $TMPDIR/merge.txt
            fi
        done

        # save to output
        cat $TMPDIR/merge.txt | sort | uniq > {output.peak_contrast_files}
        rm $TMPDIR/merge.txt
        """

rule go_enrichment:
    """
    https://bioconductor.org/packages/devel/bioc/vignettes/chipenrich/inst/doc/chipenrich-vignette.html#peak-distance-to-tss-distribution
    """
    input:
        contrast_file=rules.create_contrast_peakcaller_files.output.peak_contrast_files
    output:
        html=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","go_enrichment","{contrast_list}.{dupstatus}.go_enrichment.html"),
    params:
        rscript_wrapper=join(SCRIPTSDIR,"_go_enrichment_wrapper.R"),
        rmd=join(SCRIPTSDIR,"_go_enrichment.Rmd"),
        rscript_functions=join(SCRIPTSDIR,"_carlisle_functions.R"),
        output_dir = join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","go_enrichment"),
        species = config["genome"],
        geneset_id = GENESET_ID
    envmodules:
        TOOLS["R"]
    shell:
        """
        set -exo pipefail
        if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then 
            TMPDIR="/lscratch/$SLURM_JOB_ID"
        else
            TMPDIR="/dev/shm"
        fi

        # get sample list
        sample_list=`awk '{{print $3}}' {input.contrast_file}`

        # rum script       
        Rscript {params.rscript_wrapper} \\
            --rmd {params.rmd} \\
            --sourcefile {params.rscript_functions} \\
            --output_dir {params.output_dir} \\
            --tmpdir $TMPDIR \\
            --report {output.html} \\
            --peak_list "$sample_list" \\
            --species {params.species} \\
            --geneset_id {params.geneset_id}
        """
