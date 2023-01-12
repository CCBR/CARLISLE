rule qc_fastqc:
    """
    Runs FastQC report on each sample after adaptors have been removed
    """
    input:
        R1 = rules.trim.output.R1,
        R2 = rules.trim.output.R2,
    params:
        base = join(RESULTSDIR, 'qc', 'fastqc_raw'),
    envmodules:
        TOOLS['fastqc']
    output:
        htmlR1 = temp(join(RESULTSDIR, 'qc','fastqc_raw','{replicate}.R1.trim_fastqc.html')),
        htmlR2 = temp(join(RESULTSDIR, 'qc','fastqc_raw','{replicate}.R2.trim_fastqc.html'))
    shell:
        """
        set -exo pipefail
        # run FASTQC
        fastqc {input.R1} {input.R2} -o {params.base}
        """


rule qc_fastq_screen_validator:
    """
    #fastq screen
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
        conf_species = join(WORKDIR,'config','fqscreen_config.conf'),
        base_screen = join(RESULTSDIR, 'qc', 'fqscreen_raw'),
        R1="{replicate}_R1.fastq",
        R2="{replicate}_R2.fastq",
        base_val = join(RESULTSDIR,'qc'),
        fastq_val=TOOLS["fastq_val"],
    threads: getthreads("qc_fastq_screen_validator")        
    envmodules:
        TOOLS['bowtie2'],
        TOOLS['perl'],
        TOOLS['fastq_screen'],
    output:
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

rule multiqc:
    """
    merges FastQC reports for pre/post trimmed fastq files into MultiQC report
    https://multiqc.info/docs/#running-multiqc
    """
    input:
        f1=expand(rules.qc_fastqc.output.htmlR1,replicate=REPLICATES),
        f2=expand(rules.qc_fastqc.output.htmlR2,replicate=REPLICATES),
        f3=expand(rules.qc_fastq_screen_validator.output.speciesR1,replicate=REPLICATES),
        f4=expand(rules.qc_fastq_screen_validator.output.speciesR2,replicate=REPLICATES),
        #f3=expand(rules.align.output.bamflagstat,replicate=REPLICATES),
        #f4=expand(rules.align.output.bamidxstats,replicate=REPLICATES)
    params:
        qc_config = join(WORKDIR,'config','multiqc_config.yaml'),
        dir_fqc = join(RESULTSDIR, 'qc', 'fastqc_raw'),
        dir_fqscreen = join(RESULTSDIR, 'qc', 'fqscreen_raw'),
        dir_samtools = join(RESULTSDIR,"bam","raw"),
        outDir = join(RESULTSDIR,'qc'),
    envmodules:
        TOOLS['multiqc']
    output:
        o1 = join(RESULTSDIR,'qc', 'multiqc_report.html')
    shell:
        """
        set -exo pipefail
        multiqc -f -v \\
            -c {params.qc_config} \\
            -d -dd 1 \\
            {params.dir_fqc} \\
            {params.dir_fqscreen} \\
            -o {params.outDir}
        """
