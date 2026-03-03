import re

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
def get_deg_bed(wildcards):
    # DEG-based peak sets produced by diffbb
    # method: AUCbased | fragmentsbased
    # group: up_group1 | up_group2
    bed = join(
        RESULTSDIR,
        "peaks",
        wildcards.qthresholds,
        "contrasts",
        wildcards.control_mode,
        wildcards.contrast_list + "." + wildcards.dupstatus,
        wildcards.contrast_list
        + "."
        + wildcards.dupstatus
        + "."
        + wildcards.peak_caller_type
        + "."
        + wildcards.method
        + "_"
        + wildcards.group
        + ".bed",
    )
    return bed


localrules: homer_annotations, combine_homer
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
        # set tmp - create directory directly in target location
        if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then
            TMPDIR=$(mktemp -d "/lscratch/$SLURM_JOB_ID/tmp.XXXXXX")
        else
            TMPDIR=$(mktemp -d "/dev/shm/tmp.XXXXXX")
        fi

        preparsedDir="$TMPDIR/preparsedDir"
        mkdir -p $preparsedDir        
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
                findMotifsGenome.pl {input.peak_file} {params.fa} {params.outDir} \\
                    -nomotif \\
                    -size given \\
                    -mknown {params.hocomoco_motif} \\
                    -p {threads} \\
                    -dumpFasta -cpg -maxN 0.1 -len 10 \\
                    -preparsedDir $preparsedDir 
            else
                echo "DEBUG: Running findMotifsGenome.pl with standard genome"
                findMotifsGenome.pl {input.peak_file} {params.genome} {params.outDir} \\
                    -nomotif \\
                    -size given \\
                    -mknown {params.hocomoco_motif} \\
                    -p {threads} \\
                    -dumpFasta -cpg -maxN 0.1 -len 10 \\
                    -preparsedDir $preparsedDir 
            fi
            echo "DEBUG: findMotifsGenome.pl completed"
        fi
        
        echo "=========================================="
        echo "DEBUG: HOMER motif analysis completed successfully"
        echo "=========================================="
        """

rule homer_motif_deg:
    """
    HOMER motif discovery on DEG-based peak sets (AUC/FRAG up in group1/group2)
    """
    input:
        deg_peak_file=get_deg_bed,
    output:
        annotation=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.{method,(AUCbased|fragmentsbased)}_{group,(up_group1|up_group2)}.motifs","annotation.txt"),
        annotation_summary=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.{method,(AUCbased|fragmentsbased)}_{group,(up_group1|up_group2)}.motifs","annotation.summary"),
        known_html=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.{method,(AUCbased|fragmentsbased)}_{group,(up_group1|up_group2)}.motifs","knownResults.html"),
        target_fasta=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.{method,(AUCbased|fragmentsbased)}_{group,(up_group1|up_group2)}.motifs","target.fa"),
        background_fasta=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.{method,(AUCbased|fragmentsbased)}_{group,(up_group1|up_group2)}.motifs","background.fa"),
    threads: getthreads("homer_motif_deg")
    envmodules:
        TOOLS["homer"],
    params:
        genome = config["genome"],
        fa=config["reference"][config["genome"]]["fa"],
        gtf = config["reference"][config["genome"]]["gtf"],
        outDir=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.{method}_{group}.motifs"),
        hocomoco_motif = config["hocomoco_motifs"]
    shell:
        """
        set -euo pipefail

        # set tmp - create directory directly in target location
        if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then
            TMPDIR=$(mktemp -d "/lscratch/$SLURM_JOB_ID/tmp.XXXXXX")
        else
            TMPDIR=$(mktemp -d "/dev/shm/tmp.XXXXXX")
        fi

        preparsedDir="$TMPDIR/preparsedDir"
        mkdir -p $preparsedDir

        echo "=========================================="
        echo "DEBUG: Starting HOMER motif analysis (DEG)"
        echo "DEBUG: DEG Peak file: {input.deg_peak_file}"
        echo "DEBUG: Genome: {params.genome}"
        echo "DEBUG: Output directory: {params.outDir}"
        echo "DEBUG: Threads: {threads}"
        echo "=========================================="

        num_peaks=$(wc -l < {input.deg_peak_file} || echo 0)
        echo "DEBUG: Number of DEG peaks detected: $num_peaks"

        if [[ $num_peaks -lt 5 ]]; then
            echo "WARNING: Only $num_peaks peaks found in {input.deg_peak_file}"
            echo "INFO: Skipping HOMER analysis and creating empty output files"

            echo "# No peaks found for HOMER annotation" > {output.annotation}
            echo -e "Annotation\tDistance to TSS\tNumber of Peaks\t% of Peaks\tTotal size (bp)\tLog10 p-value\tLog2 Ratio (vs. Genome)\tLogP enrichment (+values depleted)" > {output.annotation_summary}

            mkdir -p {params.outDir}
            echo "<html><body><h1>No peaks available for motif analysis</h1><p>DEG peak file contained fewer than 5 peaks.</p></body></html>" > {output.known_html}
            touch {output.target_fasta}
            touch {output.background_fasta}
            echo "DEBUG: Empty output files created successfully"
        else
            echo "INFO: Found $num_peaks DEG peaks, proceeding with HOMER analysis..."

            echo "DEBUG: STEP 1 - Running HOMER peak annotation (DEG)"
            if [[ {params.genome} == "hs1" ]]; then
                echo "DEBUG: Using hs1 genome with custom FA and GTF"
                echo "DEBUG: FA file: {params.fa}"
                echo "DEBUG: GTF file: {params.gtf}"
                annotatePeaks.pl {input.deg_peak_file} {params.fa} -annStats {output.annotation_summary} -gtf {params.gtf} > {output.annotation}
            else
                echo "DEBUG: Using standard genome: {params.genome}"
                annotatePeaks.pl {input.deg_peak_file} {params.genome} -annStats {output.annotation_summary} > {output.annotation}
            fi
            echo "DEBUG: HOMER annotation (DEG) completed"

            echo "DEBUG: STEP 2 - Running findMotifsGenome.pl (DEG)"
            echo "DEBUG: Motif size: given (use peak widths)"
            echo "DEBUG: HOCOMOCO motif file: {params.hocomoco_motif}"
            if [[ -f {params.hocomoco_motif} ]]; then
                echo "DEBUG: HOCOMOCO motif file found"
            else
                echo "ERROR: HOCOMOCO motif file NOT found at {params.hocomoco_motif}"
            fi

            if [[ {params.genome} == "hs1" ]]; then
                findMotifsGenome.pl {input.deg_peak_file} {params.fa} {params.outDir} \\
                    -nomotif \\
                    -size given \\
                    -mknown {params.hocomoco_motif} \\\
                    -p {threads} \\
                    -dumpFasta -cpg -maxN 0.1 -len 10 \\
                    -preparsedDir $preparsedDir 
            else
                findMotifsGenome.pl {input.deg_peak_file} {params.genome} {params.outDir} \\
                    -nomotif \\
                    -size given \\
                    -mknown {params.hocomoco_motif} \\
                    -p {threads} \\
                    -dumpFasta -cpg -maxN 0.1 -len 10 \\
                    -preparsedDir $preparsedDir
            fi
            echo "DEBUG: findMotifsGenome.pl (DEG) completed"
        fi

        echo "=========================================="
        echo "DEBUG: HOMER motif analysis (DEG) completed successfully"
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

        # set tmp
            dirname=$(basename $(mktemp))
        if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then
            TMPDIR="/lscratch/$SLURM_JOB_ID/$dirname"
        else

            TMPDIR="/dev/shm/$dirname"
        fi
        mkdir -p $TMPDIR
        
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

rule ame_motif_enrichment_deg:
    """
    AME (Analysis of Motif Enrichment) on DEG-based motif FASTAs using HOCOMOCO v14 CORE
    """
    input:
        target_fasta=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.{method}_{group}.motifs","target.fa"),
        background_fasta=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.{method}_{group}.motifs","background.fa"),
    output:
        ame_results=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.{method,(AUCbased|fragmentsbased)}_{group,(up_group1|up_group2)}.motifs","ame_results.txt"),
    threads: getthreads("ame_motif_enrichment_deg")
    wildcard_constraints:
        method="AUCbased|fragmentsbased",
        group="up_group1|up_group2"
    envmodules:
        TOOLS["parallel"],
        TOOLS["meme"],
    params:
        outDir=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.{method}_{group}.motifs"),
        hocomoco_memes_tar = config["hocomoco_memes_targz"],
        python_script=join(SCRIPTSDIR,"_parse_ame_output.py")
    shell:
        """
        set -euo pipefail

        echo "=========================================="
        echo "DEBUG: Starting AME motif enrichment analysis (DEG)"
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
            echo "DEBUG: STEP 1 - Preparing FASTA files for AME analysis (DEG)"
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
            echo "DEBUG: STEP 2 - Running AME motif enrichment analysis (DEG)"
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
                    find . -name 'ame.tsv' -exec cat {{}} \; | \
                    grep -A1 ^rank | \
                    grep -v '^--$' | \
                    grep -v ^rank | \
                    sort | \
                    uniq | \
                    sort -k7,7g | \
                    python {params.python_script} >> {output.ame_results}
                    
                    echo "DEBUG: Processed results from AME outputs"
                    final_lines=$(wc -l < {output.ame_results})
                    echo "DEBUG: Final AME results file has $final_lines lines"
                    
                else
                    echo "WARNING: Prerequisites not met for AME analysis (DEG)"
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
        echo "DEBUG: AME motif enrichment analysis (DEG) completed successfully"
        echo "=========================================="
        """

def get_annotation_files(wildcards):
    """
    treatment_control_list depends on the peak caller
    """
    def control_has_replicate(pair):
        parts = pair.split("_vs_")
        if len(parts) != 2:
            return False
        return re.search(r"_[0-9]+$", parts[1]) is not None

    base_list = TREATMENT_LIST_M if wildcards.peak_caller.startswith('macs2') else TREATMENT_LIST_SG
    if wildcards.control_mode == "pooled":
        tc_list = [p for p in base_list if not control_has_replicate(p)]
    else:
        tc_list = [p for p in base_list if control_has_replicate(p)]

    return expand(join(RESULTSDIR,"peaks", wildcards.qthresholds, wildcards.peak_caller, "annotation","homer", wildcards.control_mode,
           "{treatment_control_list}" + "." + wildcards.dupstatus + "." + wildcards.peak_caller_type + ".annotation.summary"),
           treatment_control_list = tc_list)


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
    Run ROSE with a containerized two-step flow:
    1) `run-prep-rose` normalizes peak BED inputs (MACS2/SEACR/GoPeaks),
       removes TSS-overlapping regions, stitches peaks, and writes a stitched GFF.
    2) `ROSE_main.py` is executed from `/opt/ROSE` with:
       `-g <genome> -i <stitched.gff> -r <treatment.bam> [-c <control.bam>] -s -t -o`.

    Control BAM behavior:
    - If `macs2_control == "N"` and peak caller is MACS2, ROSE runs without `-c`.
    - Otherwise ROSE runs with control BAM (`-c`), including pooled controls.

    Required outputs:
    - `rose_input_Enhancers_withSuper.bed`
    - `rose_input_Gateway_Enhancers.bed`
    - `rose_input_Gateway_SuperEnhancers.bed`
    - `rose_input_AllEnhancers.table.txt`
    - Cleans up large ROSE intermediates (`gff/`, `mappedGFF/`) after success.
    """
    input:
        peak_file=get_peak_file,
        bam = expand(join(RESULTSDIR,"bam","{replicate}.{dupstatus}.bam"),replicate=REPLICATES,dupstatus=DUPSTATUS),
    threads: getthreads("rose")
    params:
        genome = config["genome"],
        tss_bed = config["reference"][config["genome"]]["tss_bed"],
        stitch_distance = config["stitch_distance"],
        tss_distance=config["tss_distance"],
        tc_file="{treatment_control_list}",
        peak_caller_type="{peak_caller_type}",
        dupstatus = "{dupstatus}",
        bam_path=join(RESULTSDIR,"bam"),
        file_base=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}"),
        control_flag = config["macs2_control"],
        rose_root="/opt/ROSE",
        rose_python="/opt/conda/envs/rose/bin/python",
        prep_bed_name="rose_input.prepared.stitched.bed",
        prep_gff_name="rose_input.prepared.stitched.gff",
    output:
        enh_with_super=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","rose_input_Enhancers_withSuper.bed"),
        gateway_enhancers=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","rose_input_Gateway_Enhancers.bed"),
        gateway_super=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","rose_input_Gateway_SuperEnhancers.bed"),
        all_enhancers_table=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","rose","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_caller_type}.{s_dist}","rose_input_AllEnhancers.table.txt"),
    container: config["containers"].get("rose", "docker://nciccbr/ccbr_rose:v1")
    shell:
        """
        set -euo pipefail
        if [[ "{params.genome}" != "hg19" && "{params.genome}" != "hg38" && "{params.genome}" != "mm10" ]]; then
            echo "ERROR: rule rose supports only hg19, hg38, and mm10. Found genome={params.genome}"
            exit 1
        fi

        # pull treatment and control ids
        treatment=`echo {params.tc_file} | awk -F"_vs_" '{{print $1}}'`
        control=`echo {params.tc_file} | awk -F"_vs_" '{{print $2}}'`

        mkdir -p {params.file_base}
        # set bam files
        treat_bam={params.bam_path}/${{treatment}}.{params.dupstatus}.bam
        if [[ "{wildcards.control_mode}" == "pooled" ]]; then
            cntrl_bam={params.bam_path}/pooled_controls/${{control}}.{params.dupstatus}.merged.bam
        else
            cntrl_bam={params.bam_path}/${{control}}.{params.dupstatus}.bam
        fi

        # Set control bam usage for ROSE based on MACS2 control mode
        control_arg=""
        if [[ "{params.control_flag}" == "N" ]] && [[ "{params.peak_caller_type}" == "macs2_narrow" || "{params.peak_caller_type}" == "macs2_broad" ]]; then
            echo "ROSE control BAM omitted for MACS2 with macs2_control=N"
        else
            control_arg="--control-bam ${{cntrl_bam}}"
        fi

        # STEP 1: prep ROSE input (BED/GFF), keeping intermediates for no_tss output.
        run-prep-rose \
            --peak-file {input.peak_file} \
            --peak-format auto \
            --sample-id {wildcards.treatment_control_list} \
            --treatment-bam ${{treat_bam}} \
            ${{control_arg}} \
            --genome {params.genome} \
            --tss-bed {params.tss_bed} \
            --stitch-distance {params.stitch_distance} \
            --tss-distance {params.tss_distance} \
            --output-dir {params.file_base} \
            --prepared-bed-name {params.prep_bed_name} \
            --prepared-gff-name {params.prep_gff_name}

        # If there are more than 5 peaks, run ROSE.
        prep_bed={params.file_base}/{params.prep_bed_name}
        prep_gff={params.file_base}/{params.prep_gff_name}
        num_of_peaks=`cat "$prep_bed" | wc -l`
        if [[ ${{num_of_peaks}} -gt 5 ]]; then
            echo "More than 5 usable peaks detected (${{num_of_peaks}}). Running ROSE."
            cd {params.rose_root}
            if [[ -n "$control_arg" ]]; then
                {params.rose_python} ROSE_main.py \
                    -g {params.genome} \
                    -i "$prep_gff" \
                    -r ${{treat_bam}} \
                    -c ${{cntrl_bam}} \
                    -s {params.stitch_distance} \
                    -t {params.tss_distance} \
                    -o {params.file_base}
            else
                {params.rose_python} ROSE_main.py \
                    -g {params.genome} \
                    -i "$prep_gff" \
                    -r ${{treat_bam}} \
                    -s {params.stitch_distance} \
                    -t {params.tss_distance} \
                    -o {params.file_base}
            fi

            # Cleanup large ROSE intermediates not needed downstream.
            rm -rf {params.file_base}/gff {params.file_base}/mappedGFF
        else
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})"
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.enh_with_super}
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.gateway_enhancers}
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.gateway_super}
            echo "Less than 5 usable peaks detected (N=${{num_of_peaks}})" > {output.all_enhancers_table}
        fi
    """
if config["run_go_enrichment"]:
    rule go_enrichment_peaks:
        """
        Run GO enrichment on all peak BED files.
        """
        wildcard_constraints:
            peak_type="narrow|broad|stringent|relaxed",
            geneset_id="[^/]+",
            method="[^/]+",
        input:
            peaks_bed=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","peak_output","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_type}.peaks.bed")
        params:
            rscript=join(SCRIPTSDIR,"_get_enrichment.R"),
            genome=config["genome"],
            output_tsv_base=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","go_enrichment","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_type}.go_enrichment.tsv")
        output:
            tsv=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","go_enrichment","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_type}.go_enrichment.{geneset_id}.{method}.tsv")
        log:
            runlog=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","go_enrichment","{control_mode}","logs","{treatment_control_list}.{dupstatus}.{peak_type}.go_enrichment.{geneset_id}.{method}.log")
        threads: getthreads("go_enrichment_peaks")
        container: config['containers']['go_enrichment']
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname "{log.runlog}")"
            {{
                echo "[$(date -Is)] START go_enrichment_peaks"
                echo "host=$(hostname) slurm_job_id=${{SLURM_JOB_ID:-NA}}"
                echo "input_peaks={input.peaks_bed}"
                echo "output_expected={output.tsv}"
                echo "output_base={params.output_tsv_base}"
                echo "geneset_id={wildcards.geneset_id} method={wildcards.method} n_cores={threads}"

                if ! Rscript {params.rscript} \\
                    --peaks_bed {input.peaks_bed} \\
                    --geneset_id {wildcards.geneset_id} \\
                    --genome {params.genome} \\
                    --output_tsv {params.output_tsv_base} \\
                    --methods {wildcards.method} \\
                    --n_cores {threads}; then
                    echo "[$(date -Is)] ERROR: Rscript failed"
                    exit 1
                fi

                outdir="$(dirname "{params.output_tsv_base}")"
                outprefix="$(basename "{params.output_tsv_base}" .tsv)"
                echo "Generated TSV candidates:"
                ls -lh "$outdir"/"$outprefix".*.tsv 2>/dev/null || echo "(none found)"

                if [[ ! -f "{output.tsv}" ]]; then
                    echo "[$(date -Is)] ERROR: Expected output missing: {output.tsv}"
                    exit 1
                fi

                echo "[$(date -Is)] DONE go_enrichment_peaks"
            }} > >(tee -a {log.runlog}) 2>&1
            """

if config["run_contrasts"] and config["run_go_enrichment"]:
    rule go_enrichment_diffbed:
        """
        Run GO enrichment on differential peak BED files.
        """
        wildcard_constraints:
            diff_type="AUCbased_up_group1|AUCbased_up_group2|fragmentsbased_up_group1|fragmentsbased_up_group2",
            peak_caller_type="[^/]+",
            geneset_id="[^/]+",
            method="[^/]+",
        input:
            peaks_bed=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.{diff_type}.bed")
        params:
            rscript=join(SCRIPTSDIR,"_get_enrichment.R"),
            genome=config["genome"],
            output_tsv_base=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","go_enrichment","{contrast_list}.{dupstatus}.{peak_caller_type}.{diff_type}.go_enrichment.tsv")
        output:
            tsv=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","go_enrichment","{contrast_list}.{dupstatus}.{peak_caller_type}.{diff_type}.go_enrichment.{geneset_id}.{method}.tsv")
        log:
            runlog=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","go_enrichment","logs","{contrast_list}.{dupstatus}.{peak_caller_type}.{diff_type}.go_enrichment.{geneset_id}.{method}.log")
        threads: getthreads("go_enrichment_diffbed")
        container: config['containers']['go_enrichment']
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname "{log.runlog}")"
            {{
                echo "[$(date -Is)] START go_enrichment_diffbed"
                echo "host=$(hostname) slurm_job_id=${{SLURM_JOB_ID:-NA}}"
                echo "input_peaks={input.peaks_bed}"
                echo "output_expected={output.tsv}"
                echo "output_base={params.output_tsv_base}"
                echo "geneset_id={wildcards.geneset_id} method={wildcards.method} n_cores={threads}"

                if ! Rscript {params.rscript} \\
                    --peaks_bed {input.peaks_bed} \\
                    --geneset_id {wildcards.geneset_id} \\
                    --genome {params.genome} \\
                    --output_tsv {params.output_tsv_base} \\
                    --methods {wildcards.method} \\
                    --n_cores {threads}; then
                    echo "[$(date -Is)] ERROR: Rscript failed"
                    exit 1
                fi

                outdir="$(dirname "{params.output_tsv_base}")"
                outprefix="$(basename "{params.output_tsv_base}" .tsv)"
                echo "Generated TSV candidates:"
                ls -lh "$outdir"/"$outprefix".*.tsv 2>/dev/null || echo "(none found)"

                if [[ ! -f "{output.tsv}" ]]; then
                    echo "[$(date -Is)] ERROR: Expected output missing: {output.tsv}"
                    exit 1
                fi

                echo "[$(date -Is)] DONE go_enrichment_diffbed"
            }} > >(tee -a {log.runlog}) 2>&1
            """

if config["run_go_enrichment"]:
    rule go_enrichment_dotplot_peaks:
        """
        Generate one GO enrichment dot plot PNG from one GO enrichment TSV (peak-call outputs).
        """
        wildcard_constraints:
            peak_type="narrow|broad|stringent|relaxed",
            geneset_id="[^/]+",
            method="[^/]+",
        input:
            tsv=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","go_enrichment","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_type}.go_enrichment.{geneset_id}.{method}.tsv")
        output:
            png=join(RESULTSDIR,"peaks","{qthresholds}","{peak_caller}","annotation","go_enrichment","{control_mode}","{treatment_control_list}.{dupstatus}.{peak_type}.go_enrichment.{geneset_id}.{method}.png")
        params:
            rscript=join(SCRIPTSDIR,"_dotplot_enrichment.R")
        container: config['containers']['go_enrichment']
        shell:
            """
            set -euo pipefail
            Rscript {params.rscript} --input {input.tsv} --output {output.png}
            """

if config["run_contrasts"] and config["run_go_enrichment"]:
    rule go_enrichment_dotplot_diffbed:
        """
        Generate one GO enrichment dot plot PNG from one GO enrichment TSV (differential peak outputs).
        """
        wildcard_constraints:
            diff_type="AUCbased_up_group1|AUCbased_up_group2|fragmentsbased_up_group1|fragmentsbased_up_group2",
            peak_caller_type="[^/]+",
            geneset_id="[^/]+",
            method="[^/]+",
        input:
            tsv=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","go_enrichment","{contrast_list}.{dupstatus}.{peak_caller_type}.{diff_type}.go_enrichment.{geneset_id}.{method}.tsv")
        output:
            png=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","go_enrichment","{contrast_list}.{dupstatus}.{peak_caller_type}.{diff_type}.go_enrichment.{geneset_id}.{method}.png")
        log:
            runlog=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{control_mode}","go_enrichment","logs","{contrast_list}.{dupstatus}.{peak_caller_type}.{diff_type}.go_enrichment.{geneset_id}.{method}.dotplot.log")
        params:
            rscript=join(SCRIPTSDIR,"_dotplot_enrichment.R")
        container: config['containers']['go_enrichment']
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname "{log.runlog}")"
            {{
                echo "[$(date -Is)] START go_enrichment_dotplot_diffbed"
                echo "host=$(hostname) slurm_job_id=${{SLURM_JOB_ID:-NA}}"
                echo "input_tsv={input.tsv}"
                echo "output_png={output.png}"

                tmp_err="$(mktemp)"
                if ! Rscript {params.rscript} --input {input.tsv} --output {output.png} 2>"$tmp_err"; then
                    cat "$tmp_err"
                    if grep -Eq "No plottable rows after filtering|No enriched pathways found|No enriched pathways with valid p-values|No rows with positive geneset size|No rows available in enrichment TSV|Input TSV is empty" "$tmp_err"; then
                        echo "[$(date -Is)] WARN: Non-plottable enrichment table; creating placeholder PNG."
                        Rscript -e "png(filename='{output.png}', width=3000, height=2100, res=300, bg='white'); plot.new(); title('No GO Enrichment Results'); text(0.5, 0.45, 'No plottable rows after filtering', cex=1.2); dev.off()"
                    else
                        echo "[$(date -Is)] ERROR: dotplot Rscript failed with unexpected error"
                        rm -f "$tmp_err"
                        exit 1
                    fi
                fi
                rm -f "$tmp_err"

                if [[ ! -f "{output.png}" ]]; then
                    echo "[$(date -Is)] ERROR: Expected output missing: {output.png}"
                    exit 1
                fi

                echo "[$(date -Is)] DONE go_enrichment_dotplot_diffbed"
            }} > >(tee -a {log.runlog}) 2>&1
            """
