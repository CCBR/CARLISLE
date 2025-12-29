def get_peak_file(wildcards):
    # MACS2 OPTIONS
    if wildcards.peak_caller_type == "macs2_narrow":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"macs2","peak_output",wildcards.control_mode, wildcards.treatment_control_list + "." + wildcards.dupstatus + ".narrow.peaks.bed")
    if wildcards.peak_caller_type == "macs2_broad":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"macs2","peak_output",wildcards.control_mode,wildcards.treatment_control_list + "." + wildcards.dupstatus + ".broad.peaks.bed")

    # SEACR OPTIONS
    if wildcards.peak_caller_type =="seacr_stringent":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr","peak_output",wildcards.control_mode,wildcards.treatment_control_list + "." + wildcards.dupstatus + ".stringent.peaks.bed")
    if wildcards.peak_caller_type =="seacr_relaxed":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"seacr","peak_output",wildcards.control_mode,wildcards.treatment_control_list + "." + wildcards.dupstatus + ".relaxed.peaks.bed")

    #GOPEAKS OPTIONS
    if wildcards.peak_caller_type =="gopeaks_narrow":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"gopeaks","peak_output",wildcards.control_mode,wildcards.treatment_control_list + "." + wildcards.dupstatus + ".narrow.peaks.bed")
    if wildcards.peak_caller_type =="gopeaks_broad":
        bed=join(RESULTSDIR,"peaks",wildcards.qthresholds,"gopeaks","peak_output",wildcards.control_mode,wildcards.treatment_control_list + "." + wildcards.dupstatus + ".broad.peaks.bed")
    return bed

localrules: create_contrast_peakcaller_files, homer_annotations, combine_homer
rule homer_motif:
    """
    HOMER peak annotation and motif discovery
    
    Developed from code: https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk

    Notes on using alternative genomes now in config
    - http://homer.ucsd.edu/homer/introduction/update.html
    - Config is located on Biowulf here: /usr/local/apps/homer/4.11.1/.//config.txt
    """
    input:
        peak_file=get_peak_file,
    output:
        annotation=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.annotation.txt"),
        annotation_summary=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.annotation.summary"),
        known_html=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.motifs","knownResults.html"),
        target_fasta=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.motifs","target.fa"),
        background_fasta=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.motifs","background.fa"),
    threads: getthreads("homer_motif")
    envmodules:
        TOOLS["homer"],
    params:
        genome = config["genome"],
        fa=config["reference"][config["genome"]]["fa"],
        gtf = config["reference"][config["genome"]]["gtf"],
        outDir=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.motifs"),
        hocomoco_motif = config["hocomoco_motifs"]
    shell:
        """
        set -euo pipefail
        
        echo "=========================================="
        echo "DEBUG: Starting HOMER motif analysis"
        echo "DEBUG: Peak file: {input.peak_file}"
        echo "DEBUG: Genome: {params.genome}"
        echo "DEBUG: Output directory: {params.outDir}"
        echo "DEBUG: Threads: {threads}"
        echo "=========================================="
        
        # Check if peak file is empty or has no peaks
        num_peaks=$(wc -l < {input.peak_file} || echo 0)
        echo "DEBUG: Number of peaks detected: $num_peaks"
        
        if [[ $num_peaks -lt 5 ]]; then
            echo "WARNING: Only $num_peaks peaks found in {input.peak_file}"
            echo "INFO: Skipping HOMER analysis and creating empty output files"
            
            # Create empty annotation files
            echo "# No peaks found for HOMER annotation" > {output.annotation}
            echo -e "Annotation\\tDistance to TSS\\tNumber of Peaks\\t% of Peaks\\tTotal size (bp)\\tLog10 p-value\\tLog2 Ratio (vs. Genome)\\tLogP enrichment (+values depleted)" > {output.annotation_summary}
            
            # Create minimal motif output directory and files
            mkdir -p {params.outDir}
            echo "<html><body><h1>No peaks available for motif analysis</h1><p>Peak file contained fewer than 5 peaks.</p></body></html>" > {output.known_html}
            touch {output.target_fasta}
            touch {output.background_fasta}
            echo "DEBUG: Empty output files created successfully"
        else
            echo "INFO: Found $num_peaks peaks, proceeding with HOMER analysis..."
            
            # ============================================
            # STEP 1: HOMER Peak Annotation
            # ============================================
            echo "DEBUG: STEP 1 - Running HOMER peak annotation"
            if [[ {params.genome} == "hs1" ]]; then
                echo "DEBUG: Using hs1 genome with custom FA and GTF"
                echo "DEBUG: FA file: {params.fa}"
                echo "DEBUG: GTF file: {params.gtf}"
                annotatePeaks.pl {input.peak_file} {params.fa} -annStats {output.annotation_summary} -gtf {params.gtf} > {output.annotation}
            else
                echo "DEBUG: Using standard genome: {params.genome}"
                annotatePeaks.pl {input.peak_file} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
            fi
            echo "DEBUG: HOMER annotation completed"
            echo "DEBUG: Annotation file: {output.annotation}"
            echo "DEBUG: Annotation summary: {output.annotation_summary}"

            # ============================================
            # STEP 2: HOMER Motif Discovery
            # ============================================
            echo "DEBUG: STEP 2 - Running findMotifsGenome.pl"
            echo "DEBUG: Motif size: given (use peak widths)"
            echo "DEBUG: HOCOMOCO motif file: {params.hocomoco_motif}"
            echo "DEBUG: Checking if HOCOMOCO motif file exists..."
            if [[ -f {params.hocomoco_motif} ]]; then
                echo "DEBUG: HOCOMOCO motif file found"
            else
                echo "ERROR: HOCOMOCO motif file NOT found at {params.hocomoco_motif}"
            fi
            
            if [[ {params.genome} == "hs1" ]]; then
                echo "DEBUG: Running findMotifsGenome.pl with hs1 genome"
                findMotifsGenome.pl {input.peak_file} {params.fa} {params.outDir} \
                    -nomotif \
                    -size given \
                    -mknown {params.hocomoco_motif} \\
                    -p {threads} \\
                    -dumpFasta -cpg -maxN 0.1 -len 10 \\
                    
            else
                echo "DEBUG: Running findMotifsGenome.pl with standard genome"
                findMotifsGenome.pl {input.peak_file} {params.genome} {params.outDir} \
                    -nomotif \
                    -size given \
                    -mknown {params.hocomoco_motif} \\
                    -p {threads} \\
                    -dumpFasta -cpg -maxN 0.1 -len 10 \\
                    
            fi
            echo "DEBUG: findMotifsGenome.pl completed"
        fi
        
        echo "=========================================="
        echo "DEBUG: HOMER motif analysis completed successfully"
        echo "=========================================="
        """

rule ame_motif_enrichment:
    """
    AME (Analysis of Motif Enrichment) using HOCOMOCO v14 CORE database
    """
    input:
        target_fasta=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.motifs","target.fa"),
        background_fasta=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.motifs","background.fa"),
    output:
        ame_results=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.motifs","ame_results.txt"),
    threads: getthreads("ame_motif_enrichment")
    envmodules:
        TOOLS["parallel"],
        TOOLS["meme"],
    params:
        outDir=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.motifs"),
        hocomoco_memes_tar = config["hocomoco_memes_targz"],
        python_script=join(SCRIPTSDIR,"_parse_ame_output.py")
    shell:
        """
        set -euo pipefail
        
        echo "=========================================="
        echo "DEBUG: Starting AME motif enrichment analysis"
        echo "DEBUG: Target FASTA: {input.target_fasta}"
        echo "DEBUG: Background FASTA: {input.background_fasta}"
        echo "DEBUG: Output directory: {params.outDir}"
        echo "DEBUG: Threads: {threads}"
        echo "=========================================="
        
        # Check if FASTA files are empty
        target_lines=$(wc -l < {input.target_fasta} 2>/dev/null || echo 0)
        background_lines=$(wc -l < {input.background_fasta} 2>/dev/null || echo 0)
        
        if [[ $target_lines -eq 0 ]] || [[ $background_lines -eq 0 ]]; then
            echo "WARNING: Empty FASTA files detected (target: $target_lines, background: $background_lines)"
            echo "INFO: Skipping AME analysis and creating empty results file"
            echo "# No sequences for AME analysis" > {output.ame_results}
        else
            echo "INFO: FASTA files ready (target: $target_lines lines, background: $background_lines lines)"
            
            # ============================================
            # STEP 1: Prepare FASTA files for AME
            # ============================================
            echo "DEBUG: STEP 1 - Preparing FASTA files for AME analysis"
            cd {params.outDir}
            echo "DEBUG: Changed directory to {params.outDir}"
            echo "DEBUG: Current working directory: $(pwd)"

            # Clean and create tmpdir
            if [[ -d tmpdir ]] ; then
                echo "DEBUG: Removing existing tmpdir"
                rm -rf tmpdir
            fi
            mkdir -p tmpdir
            echo "DEBUG: Created tmpdir"

            
            # Check if HOMER generated fasta files
            if [[ -f {params.outDir}/target.fa ]]; then
                echo "DEBUG: target.fa exists ($(wc -l < {params.outDir}/target.fa) lines)"
            cp {params.outDir}/target.fa tmpdir/target.fa 
            else
                echo "WARNING: target.fa NOT found"
touch {params.outDir}/target.fa
            fi
            
            if [[ -f {params.outDir}/background.fa ]]; then
                echo "DEBUG: background.fa exists ($(wc -l < {params.outDir}/background.fa) lines)"
            cp {params.outDir}/background.fa tmpdir/background.fa
            else
                echo "WARNING: background.fa NOT found"
touch {params.outDir}/background.fa
            fi

            # Copy fasta files to tmpdir
            echo "DEBUG: Copied fasta files to tmpdir"

            # ============================================
            # STEP 2: AME Motif Enrichment Analysis
            # ============================================
            echo "DEBUG: STEP 2 - Running AME motif enrichment analysis"
            cd tmpdir
            echo "DEBUG: Changed to tmpdir: $(pwd)"
            
            # Initialize AME results file with header
            printf '%b\n' 'rank\tmotif_DB\tmotif_ID\tmotif_ALT_ID\tconsensus\tp-value\tadjusted-p-value\tE-value\ttests\tFAMP\tn_sequences\tTP\t%TP\tFP\t%FP' > {output.ame_results}
            echo "DEBUG: Created AME results file with header"

            # Extract and process HOCOMOCO meme files
            echo "DEBUG: Checking for HOCOMOCO meme tar.gz: {params.hocomoco_memes_tar}"
            if [[ -f {params.hocomoco_memes_tar} ]]; then
                echo "DEBUG: HOCOMOCO meme tar.gz found, extracting..."
                cp {params.hocomoco_memes_tar} .
                tar xzf $(basename {params.hocomoco_memes_tar})
                echo "DEBUG: Extraction complete"
                
                # Create list of meme files
                echo "DEBUG: Creating list of meme files..."
                ls *.meme 2>/dev/null | sort > memes || touch memes
                num_meme_files=$(wc -l < memes 2>/dev/null || echo 0)
                echo "DEBUG: Found $num_meme_files meme files"
                
                # Check fasta file availability
                echo "DEBUG: Checking fasta files in tmpdir..."
                if [[ -f target.fa ]]; then
                    echo "DEBUG: target.fa found in tmpdir ($(wc -l < target.fa) lines)"
                else
                    echo "WARNING: target.fa NOT found in tmpdir"
                fi
                
                if [[ -f background.fa ]]; then
                    echo "DEBUG: background.fa found in tmpdir ($(wc -l < background.fa) lines)"
                else
                    echo "WARNING: background.fa NOT found in tmpdir"
                fi
                
                if [[ -s memes ]] && [[ -f target.fa ]] && [[ -f background.fa ]]; then
                    echo "DEBUG: All prerequisites met, generating AME commands..."
                    
                    # Generate AME commands
                    while read a; do
                        echo "ame --o ${{a}}_ame_out --noseq --control background.fa --seed 12345 --verbose 3 target.fa ${{a}}"
                    done < memes > do_memes
                    
                    num_commands=$(wc -l < do_memes)
                    echo "DEBUG: Generated $num_commands AME commands"
                    
                    # Run AME in parallel
                    echo "DEBUG: Running AME in parallel with {threads} threads..."
                    parallel -j {threads} < do_memes
                    echo "DEBUG: AME parallel execution completed"
                    
                    # Collect and process AME results
                    echo "DEBUG: Collecting AME results..."
                    find . -name 'ame.tsv' -exec cat {{}} \\; | \\
                    grep -A1 ^rank | \\
                    grep -v '^--$' | \\
                    grep -v ^rank | \\
                    sort | \\
                    uniq | \\
                    sort -k7,7g | \\
                    python {params.python_script} >> {output.ame_results}
                    
                    echo "DEBUG: Processed results from AME outputs"
                    final_lines=$(wc -l < {output.ame_results})
                    echo "DEBUG: Final AME results file has $final_lines lines"
                    
                else
                    echo "WARNING: Prerequisites not met for AME analysis"
                    echo "DEBUG: memes file size: $(wc -l < memes 2>/dev/null || echo 0)"
                    echo "# No meme files or fasta files found for AME analysis" > {output.ame_results}
                fi
            else
                echo "ERROR: HOCOMOCO meme files not found at {params.hocomoco_memes_tar}"
                echo "# HOCOMOCO meme files not found at {params.hocomoco_memes_tar}" > {output.ame_results}
            fi
            
            cd {params.outDir}
            echo "DEBUG: Returned to {params.outDir}"
            rm -rf tmpdir
            echo "DEBUG: Deleted tmpdir"
        fi
        
        echo "=========================================="
        echo "DEBUG: AME motif enrichment analysis completed successfully"
        echo "=========================================="
        """

def get_annotation_files(wildcards):
    """
    treatment_control_list depends on the peak caller
    """
    return expand(join(RESULTSDIR,"peaks", wildcards.qthresholds, wildcards.peak_caller, "annotation","homer", wildcards.control_mode,
           "{treatment_control_list}" + "." + wildcards.dupstatus + "." + wildcards.peak_caller_type + ".annotation.summary"),
           treatment_control_list = TREATMENT_LIST_M if wildcards.peak_caller.startswith('macs2') else TREATMENT_LIST_SG)


rule homer_annotations:
    """
    Plot enrichment over genic features
    """
    input:
        annotation_summary=get_annotation_files
    output:
        enrich_png=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","enrichment.{dupstatus}.{peak_caller_type}.png")
    params:
        annotation_dir=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}"),
        peak_mode="{peak_caller_type}",
        dupstatus="{dupstatus}",
        rscript=join(SCRIPTSDIR,"_plot_feature_enrichment.R")
    container: config['containers']['carlisle_r']
    shell:
        """
        Rscript {params.rscript} {params.annotation_dir} {params.peak_mode} {params.dupstatus} {output.enrich_png}
        """

rule combine_homer:
    """
    Add MACS2 q-value and FC to HOMER peak annotation
    """
    input:
        annotation=join(RESULTSDIR,"peaks","{qthresholds}","macs2","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.annotation.txt"),
        peaks_file=join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.peaks.xls")
    output:
        combined_tsv=join(RESULTSDIR,"peaks","{qthresholds}","macs2","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.annotation_qvalue.tsv"),
        combined_xlsx=join(RESULTSDIR,"peaks","{qthresholds}","macs2","annotation","homer","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.annotation_qvalue.xlsx")
    container: config['containers']['carlisle_r']
    params:
        rscript=join(SCRIPTSDIR,"_combine_macs2_homer.R")
    shell:
        """
        Rscript {params.rscript} {input.peaks_file} {input.annotation} {output.combined_tsv} {output.combined_xlsx}
        """

rule rose:
    """
    Developed from code:
    https://github.com/CCRGeneticsBranch/khanlab_pipeline/blob/master/rules/pipeline.chipseq.smk
    outdated version of rose: https://github.com/younglab/ROSE

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
        TOOLS["python3"],
        TOOLS["samtools"],
        TOOLS["R"]
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
        if [[ {params.peak_caller_type} == "seacr_stringent" ]] || [[ {params.peak_caller_type} == "seacr_relaxed" ]]; then
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
if config["run_contrasts"]:
    rule create_contrast_peakcaller_files:
        """
        Reads in all of the output from Rules create_contrast_data_files which match the same peaktype and merges them together
        """
        input:
            contrast_files=expand(join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.txt"),qthresholds=QTRESHOLDS, contrast_list=CONTRAST_LIST,dupstatus=DUPSTATUS,peak_caller_type=PEAKTYPE,control_mode=CONTROL_MODES)
        params:
            qthresholds = "{qthresholds}",
            contrast_list = "{contrast_list}",
            dupstatus = "{dupstatus}",
            peak_caller = "{peak_caller}",
            control_mode = "{control_mode}",
            search_dir=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}")
        output:
            peak_contrast_files=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","go_enrichment","{control_mode}","{contrast_list}.{dupstatus}.txt")
        shell:
            """
            cat {input.contrast_files} | sort | uniq > {output.peak_contrast_files}
            """

    rule go_enrichment:
        """
        https://bioconductor.org/packages/devel/bioc/vignettes/chipenrich/inst/doc/chipenrich-vignette.html#peak-distance-to-tss-distribution
        """
        input:
            contrast_file=rules.create_contrast_peakcaller_files.output.peak_contrast_files
        params:
            rscript_wrapper=join(SCRIPTSDIR,"_go_enrichment_wrapper.R"),
            rmd=join(SCRIPTSDIR,"_go_enrichment.Rmd"),
            carlisle_functions=join(SCRIPTSDIR,"_carlisle_functions.R"),
            rscript_diff=join(SCRIPTSDIR,"_diff_markdown_wrapper.R"),
            rscript_functions=join(SCRIPTSDIR,"_carlisle_functions.R"),
            output_dir = join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","go_enrichment","{control_mode}"),
            species = config["genome"],
            geneset_id = GENESET_ID,
            dedup_status =  "{dupstatus}"
        output:
            html=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","go_enrichment","{control_mode}","{contrast_list}.{dupstatus}.go_enrichment.html"),
        container: config['containers']['carlisle_r']
        shell:
            """
            set -exo pipefail

            # get sample list
            sample_list=`awk '{{print $3}}' {input.contrast_file}`
            clean_sample_list=`echo $sample_list | sed "s/\s/xxx/g"`

            # rum script
            Rscript {params.rscript_wrapper} \\
                --rmd {params.rmd} \\
                --carlisle_functions {params.carlisle_functions} \\
                --output_dir {params.output_dir} \\
                --report {output.html} \\
                --peak_list "$clean_sample_list" \\
                --species {params.species} \\
                --geneset_id {params.geneset_id} \\
                --dedup_status {params.dedup_status}
            """
rule motif_enrichment:
    """
    Perform motif enrichment analysis if enabled in the configuration.
    """
    input:
        annotation_summary=get_annotation_files
    output:
        enrich_png=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}","enrichment.{dupstatus}.{peak_caller_type}.png")
    params:
        annotation_dir=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","homer","{control_mode}"),
        peak_mode="{peak_caller_type}",
        dupstatus="{dupstatus}",
        rscript=join(SCRIPTSDIR,"_plot_feature_enrichment.R")
    container: config['containers']['carlisle_r']
    shell:
        """
        Rscript {params.rscript} {params.annotation_dir} {params.peak_mode} {params.dupstatus} {output.enrich_png}
        """
