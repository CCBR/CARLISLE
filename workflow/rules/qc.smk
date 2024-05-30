localrules: spikein_assessment
rule qc_fastqc:
    """
    1) Runs FastQC report on each sample after adaptors have been removed

    2) Runs FASTQ SCREEN & VALIDATOR
    - this will align first to human, mouse, bacteria
    - fastq validator
    Quality-control step to ensure the input FastQC files are not corrupted or
    incomplete prior to running the entire workflow.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        Log file containing any warnings or errors on file
    """
    input:
        R1 = rules.trim.output.R1,
        R2 = rules.trim.output.R2,
    params:
        base = join(RESULTSDIR, 'qc', 'fastqc_raw'),
        conf_species = join(WORKDIR,'config','fqscreen_config.conf'),
        base_screen = join(RESULTSDIR, 'qc', 'fqscreen_raw'),
        R1="{replicate}_R1.fastq",
        R2="{replicate}_R2.fastq",
        base_val = join(RESULTSDIR,'qc'),
        fastq_val=TOOLS["fastq_val"],
    envmodules:
        TOOLS['fastqc'],
        TOOLS['bowtie2'],
        TOOLS['perl'],
        TOOLS['fastq_screen'],
    output:
        htmlR1 = temp(join(RESULTSDIR, 'qc','fastqc_raw','{replicate}.R1.trim_fastqc.html')),
        htmlR2 = temp(join(RESULTSDIR, 'qc','fastqc_raw','{replicate}.R2.trim_fastqc.html')),
        speciesR1 = temp(join(RESULTSDIR, 'qc', 'fqscreen_raw','{replicate}_R1_screen.txt')),
        speciesR2 = temp(join(RESULTSDIR, 'qc', 'fqscreen_raw','{replicate}_R2_screen.txt')),
        logR1 = join(RESULTSDIR,'qc','{replicate}.validated.fastqR1.log'),
        logR2 = join(RESULTSDIR,'qc','{replicate}.validated.fastqR2.log'),
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

        # run FASTQC
        fastqc {input.R1} {input.R2} -o {params.base}
 
        # Gzip input files
        gunzip -c {input.R1} > ${{TMPDIR}}/{params.R1}
        gunzip -c {input.R2} > ${{TMPDIR}}/{params.R2}

        # Run FastQ Screen
        fastq_screen ${{TMPDIR}}/{params.R1} ${{TMPDIR}}/{params.R2} \\
            --conf {params.conf_species} \\
            --outdir {params.base_screen} \\
            --threads {threads} \\
            --subset 1000000 \\
            --aligner bowtie2 \\
            --force
        
        # Run FastQ Validator
        mkdir -p {params.base_val}
        {params.fastq_val} \\
            --disableSeqIDCheck \\
            --noeof \\
            --printableErrors 100000000 \\
            --baseComposition \\
            --avgQual \\
            --file {input.R1} > {output.logR1};
        {params.fastq_val} \\
            --disableSeqIDCheck \\
            --noeof \\
            --printableErrors 100000000 \\
            --baseComposition \\
            --avgQual \\
            --file {input.R2} > {output.logR2};
        """

rule spikein_assessment:
    """
    Runs assessment on bam.idxstats files to ensure that spikein control amounts are uniform amongst replicates

    @Input:
        Bamidxstats files
    @Output:
        R HTML Report of spike in control values samples
    """
    input:
        bams = expand(rules.align.output.bamidxstats,replicate=REPLICATES,dupstatus=DUPSTATUS),
    output:
        html=join(RESULTSDIR,'qc',"spikein_qc_report.html"),
    params:
        rscript_wrapper=join(SCRIPTSDIR,"_generate_spikein_wrapper.R"),
        rmd=join(SCRIPTSDIR,"_generate_spikein_plots.Rmd"),
        carlisle_functions=join(SCRIPTSDIR,"_carlisle_functions.R"),
        rscript_diff=join(SCRIPTSDIR,"_diff_markdown_wrapper.R"),
        spikein=config["spikein_genome"],
    container: config['containers']['carlisle_r']
    shell:
        """
        if [[ {params.spikein} == "ecoli" ]]; then species_name="NC_000913.3"; else species_name=""; fi
        
        # get sample list
        sample_list="{input.bams}"
        clean_sample_list=`echo $sample_list | sed "s/\s/xxx/g"`

        # rum script       
        Rscript {params.rscript_wrapper} \\
            --rmd {params.rmd} \\
            --carlisle_functions {params.carlisle_functions} \\
            --report {output.html} \\
            --bam_list "$clean_sample_list" \\
            --spikein_control $species_name
        """

if ("gopeaks_narrow" in PEAKTYPE) or ("gopeaks_broad" in PEAKTYPE):
    rule multiqc:
        """
        merges FastQC reports for pre/post trimmed fastq files into MultiQC report
        https://multiqc.info/docs/#running-multiqc
        """
        input:
            fqR1=expand(rules.qc_fastqc.output.htmlR1,replicate=REPLICATES),
            fqR2=expand(rules.qc_fastqc.output.htmlR2,replicate=REPLICATES),
            screenR1=expand(rules.qc_fastqc.output.speciesR1,replicate=REPLICATES),
            screenR2=expand(rules.qc_fastqc.output.speciesR2,replicate=REPLICATES),
            flagstat=expand(rules.align.output.bamflagstat,replicate=REPLICATES),
            idxstat=expand(rules.align.output.bamidxstats,replicate=REPLICATES),
            gopeaks_broad=expand(rules.gopeaks_broad.output.json,qthresholds=QTRESHOLDS,treatment_control_list=TREATMENT_LIST_SG,dupstatus=DUPSTATUS),
            gopeaks_narrow=expand(rules.gopeaks_narrow.output.json,qthresholds=QTRESHOLDS,treatment_control_list=TREATMENT_LIST_SG,dupstatus=DUPSTATUS)
        params:
            qc_config = join(WORKDIR,'config','multiqc_config.yaml'),
            dir_fqc = join(RESULTSDIR, 'qc', 'fastqc_raw'),
            dir_fqscreen = join(RESULTSDIR, 'qc', 'fqscreen_raw'),
            dir_samtools = join(RESULTSDIR,"bam","raw"),
            dir_gopeaks = join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","peak_output"),
        envmodules:
            TOOLS['multiqc']
        output:
            report = join(RESULTSDIR,'qc','multiqc_report_{qthresholds}.html')
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

            multiqc -f -v \\
                -c {params.qc_config} \\
                -d -dd 1 \\
                {params.dir_fqc} \\
                {params.dir_fqscreen} \\
                {params.dir_samtools} \\
                {params.dir_gopeaks} \\
                -o $TMPDIR/qc

            mv $TMPDIR/qc/multiqc_report.html {output.report}
            """
else:
    rule multiqc:
        """
        merges FastQC reports for pre/post trimmed fastq files into MultiQC report
        https://multiqc.info/docs/#running-multiqc
        """
        input:
            fqR1=expand(rules.qc_fastqc.output.htmlR1,replicate=REPLICATES),
            fqR2=expand(rules.qc_fastqc.output.htmlR2,replicate=REPLICATES),
            screenR1=expand(rules.qc_fastqc.output.speciesR1,replicate=REPLICATES),
            screenR2=expand(rules.qc_fastqc.output.speciesR2,replicate=REPLICATES),
            flagstat=expand(rules.align.output.bamflagstat,replicate=REPLICATES),
            idxstat=expand(rules.align.output.bamidxstats,replicate=REPLICATES),
        params:
            qc_config = join(WORKDIR,'config','multiqc_config.yaml'),
            dir_fqc = join(RESULTSDIR, 'qc', 'fastqc_raw'),
            dir_fqscreen = join(RESULTSDIR, 'qc', 'fqscreen_raw'),
            dir_samtools = join(RESULTSDIR,"bam","raw"),
        envmodules:
            TOOLS['multiqc']
        output:
            report = join(RESULTSDIR,'qc','multiqc_report.html')
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

            multiqc -f -v \\
                -c {params.qc_config} \\
                -d -dd 1 \\
                {params.dir_fqc} \\
                {params.dir_fqscreen} \\
                {params.dir_samtools} \\
                -o $TMPDIR/qc

            mv $TMPDIR/qc/multiqc_report.html {output.report}
            """