localrules:contrast_init

rule contrast_init:
    input:
        expand(join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.narrowPeak"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        expand(join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.broadPeak"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        expand(join(RESULTSDIR,"fragments","{replicate}.{dupstatus}.fragments.bed"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        expand([join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.dedup.norm.stringent.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        expand([join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.dedup.norm.relaxed.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        expand([join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.no_dedup.norm.stringent.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        expand([join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.no_dedup.norm.relaxed.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        expand([join(RESULTSDIR,"peaks","gopeaks","{treatment}_vs_{control}.dedup.narrowGo_peaks.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        expand([join(RESULTSDIR,"peaks","gopeaks","{treatment}_vs_{control}.dedup.broadGo_peaks.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        join(RESULTSDIR,"replicate_sample.tsv"),
        expand(join(RESULTSDIR,"bedgraph","{replicate}.{dupstatus}.sf.yaml"),replicate=REPLICATES,dupstatus=DUPSTATUS)
    output:
        outtsv=join(RESULTSDIR,"peaks","contrasts","bed_bedgraph_paths.tsv"),
    params:
        resultsdir = RESULTSDIR,
        bedgraphdir = join(RESULTSDIR,"bedgraph"),
        dedup_list=DUPSTATUS,
        peak_list=PEAKTYPE
    shell:
        """
        while read replicate sample;do
            for dupstatus in {params.dedup_list};do
                bedgraph=$(find {params.resultsdir} -name "*${{replicate}}.${{dupstatus}}.bedgraph")
                fragment_bed=$(find {params.resultsdir} -name "*${{replicate}}.${{dupstatus}}.fragments.bed")
                sf=$(cat {params.bedgraphdir}/${{replicate}}.${{dupstatus}}.sf.yaml | awk -F"=" '{{print $NF}}')
                for peaktype in {params.peak_list}; do
                    if [[ $dupstatus == "dedup" ]];then
                        bed=$(find {params.resultsdir} -name "*${{replicate}}*${{dupstatus}}*${{peaktype}}" |grep -v no_dedup)
                    else
                        bed=$(find {params.resultsdir} -name "*${{replicate}}*${{dupstatus}}*${{peaktype}}")
                    fi
                    echo -ne "$replicate\\t$sample\\t$dupstatus\\t$peaktype\\t$bed\\t$bedgraph\\t$sf\\t$fragment_bed\\n"
                done
            done
        done < {params.resultsdir}/replicate_sample.tsv > {output.outtsv}    
        """

localrules:make_inputs

rule make_inputs:
    input:
        join(RESULTSDIR,"peaks","contrasts","bed_bedgraph_paths.tsv"),
    output:
        inputs=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_inputs.txt"),
    params:
        condition1 = "{c1}",
        condition2 = "{c2}",
        ds = "{ds}",
        pt = "{pt}",
    shell:
        """
        for condition in {params.condition1} {params.condition2}
        do
            while read s g d p bed bg sf fbed;do
                if [[ "$g" == "$condition" ]];then
                    if [[ "$d" == "{params.ds}" ]];then
                        if [[ "$p" == "{params.pt}" ]];then
                            echo -ne "$g\\t$s\\t$bed\\t$bg\\t$sf\\t$fbed\\n"
                        fi
                    fi
                fi
            done < {input}
        done > {output.inputs}
        """

rule make_counts_matrix:
    input:
        inputs=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_inputs.txt"),
    output:
        cm=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_countsmatrix.txt"),
        fcm=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_fragmentscountsmatrix.txt"),
        si=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_sampleinfo.txt"),
    params:
        pyscript=join(SCRIPTSDIR,"_make_counts_matrix.py"),
        condition1 = "{c1}",
        condition2 = "{c2}",
        ds = "{ds}",
        pt = "{pt}",
    envmodules:
        TOOLS["python37"],
        TOOLS["bedtools"],
        TOOLS["bedops"]
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
        python {params.pyscript} --bedbedgraph {input.inputs} --tmpdir $TMPDIR --countsmatrix {output.cm} --fragmentscountsmatrix {output.fcm} --sampleinfo {output.si}
        """

rule DESeq:
    input:
        bbpaths=join(RESULTSDIR,"peaks","contrasts","bed_bedgraph_paths.tsv"), # this has the scaling factors
        cm=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_countsmatrix.txt"),
        si=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_sampleinfo.txt"),
    output:
        results=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_AUCbased_diffresults.txt"),
        html=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_AUCbased_diffanalysis.html"),
        elbowlimits=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_AUCbased_diffanalysis_elbowlimits.yaml"),
    params:
        rscript=join(SCRIPTSDIR,"_diff_markdown_wrapper.R"),
        rmd=join(SCRIPTSDIR,"_diff_markdown.Rmd"),
        condition1 = "{c1}",
        condition2 = "{c2}",
        ds = "{ds}",
        pt = "{pt}",
        spiked = SPIKED, # "Y" for spiked
        fdr_cutoff = FDRCUTOFF,
        log2fc_cutoff = LFCCUTOFF,
        species = config["genome"],
    envmodules:
        TOOLS["R"]
    shell:
        """
        set -exo pipefail
        dirname=$(basename $(mktemp))
        if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then 
            TMPDIR="/lscratch/$SLURM_JOB_ID/$dirname"
        else
            TMPDIR="/dev/shm/$dirname"
        fi
        mkdir -p $TMPDIR
        mkdir -p ${{TMPDIR}}/intermediates_dir
        mkdir -p ${{TMPDIR}}/knit_root_dir
        cd $TMPDIR
        Rscript {params.rscript} \\
            --rmd {params.rmd} \\
            --countsmatrix {input.cm} \\
            --sampleinfo {input.si} \\
            --dupstatus {params.ds} \\
            --condition1 {params.condition1} \\
            --condition2 {params.condition2} \\
            --fdr_cutoff {params.fdr_cutoff} \\
            --log2fc_cutoff {params.log2fc_cutoff} \\
            --results {output.results} \\
            --report {output.html} \\
            --elbowlimits {output.elbowlimits} \\
            --spiked {params.spiked} \\
            --rawcountsprescaled \\
            --scalesfbymean \\
            --bbpaths {input.bbpaths} \\
            --tmpdir $TMPDIR \\
            --species {params.species}
        """

rule DESeq2:
    input:
        bbpaths=join(RESULTSDIR,"peaks","contrasts","bed_bedgraph_paths.tsv"), # this has the scaling factors
        cm=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_fragmentscountsmatrix.txt"),
        si=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_sampleinfo.txt"),
    output:
        results=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_fragmentsbased_diffresults.txt"),
        html=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_fragmentsbased_diffanalysis.html"),
        elbowlimits=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_fragmentsbased_diffanalysis_elbowlimits.yaml"),
    params:
        rscript=join(SCRIPTSDIR,"_diff_markdown_wrapper.R"),
        rmd=join(SCRIPTSDIR,"_diff_markdown.Rmd"),
        condition1 = "{c1}",
        condition2 = "{c2}",
        ds = "{ds}",
        pt = "{pt}",
        spiked = SPIKED, # "Y" for spiked
        fdr_cutoff = FDRCUTOFF,
        log2fc_cutoff = LFCCUTOFF,
        species = config["genome"],
    envmodules:
        TOOLS["R"]
    shell:
        """
        set -exo pipefail
        dirname=$(basename $(mktemp))
        if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then 
            TMPDIR="/lscratch/$SLURM_JOB_ID/$dirname"
        else
            TMPDIR="/dev/shm/$dirname"
        fi
        mkdir -p $TMPDIR
        mkdir -p ${{TMPDIR}}/intermediates_dir
        mkdir -p ${{TMPDIR}}/knit_root_dir
        cd $TMPDIR
        # Do not use --rawcountsprescaled as these counts are not prescaled!
        Rscript {params.rscript} \\
            --rmd {params.rmd} \\
            --countsmatrix {input.cm} \\
            --sampleinfo {input.si} \\
            --dupstatus {params.ds} \\
            --condition1 {params.condition1} \\
            --condition2 {params.condition2} \\
            --fdr_cutoff {params.fdr_cutoff} \\
            --log2fc_cutoff {params.log2fc_cutoff} \\
            --results {output.results} \\
            --report {output.html} \\
            --elbowlimits {output.elbowlimits} \\
            --spiked {params.spiked} \\
            --scalesfbymean \\
            --bbpaths {input.bbpaths} \\
            --tmpdir $TMPDIR \\
            --species {params.species}
        """

rule diffbb:
    input:
        results=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_AUCbased_diffresults.txt"),
        elbowlimits=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_AUCbased_diffanalysis_elbowlimits.yaml"),
        fresults=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_fragmentsbased_diffresults.txt"),
        felbowlimits=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_fragmentsbased_diffanalysis_elbowlimits.yaml"),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
    output:
        bed=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_AUCbased_diffresults.bed"),
        bigbed=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_AUCbased_diffresults.bigbed"),
        fbed=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_fragmentsbased_diffresults.bed"),
        fbigbed=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_fragmentsbased_diffresults.bigbed"),
    params:
        fdr=FDRCUTOFF,
        lfc=LFCCUTOFF,
        randstr = RANDOMSTR,
        script = join(SCRIPTSDIR,"_make_results_bed.py"),
    envmodules:
        TOOLS["ucsc"],
        TOOLS["python37"]
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

        python {params.script} \\
            --results {input.results} \\
            --fdr_cutoff {params.fdr} \\
            --log2FC_cutoff {params.lfc} \\
            --elbowyaml {input.elbowlimits} \\
            --bed {output.bed}

        bedToBigBed -type=bed9 {output.bed} {input.genome_len} {output.bigbed}

        python {params.script} \\
            --results {input.fresults} \\
            --fdr_cutoff {params.fdr} \\
            --log2FC_cutoff {params.lfc} \\
            --elbowyaml {input.felbowlimits} \\
            --bed {output.fbed}

        bedToBigBed -type=bed9 {output.fbed} {input.genome_len} {output.fbigbed}
        """

localrules: venn

rule venn:
    input:
        results=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_AUCbased_diffresults.txt"),
        fresults=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_fragmentsbased_diffresults.txt"),
    output:
        pdf=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_venn.pdf")
    params:
        rscript=join(SCRIPTSDIR,"_plot_results_venn.R"),
        condition1 = "{c1}",
        condition2 = "{c2}",
        ds = "{ds}",
        pt = "{pt}",
        fdr_cutoff = FDRCUTOFF,
        log2fc_cutoff = LFCCUTOFF,
    envmodules:
        TOOLS["R"]
    shell:
        """
        set -exo pipefail
        dirname=$(basename $(mktemp))
        if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then 
            TMPDIR="/lscratch/$SLURM_JOB_ID/$dirname"
        else
            TMPDIR="/dev/shm/$dirname"
        fi
        mkdir -p $TMPDIR
        cd $TMPDIR
        Rscript {params.rscript} \\
            --aucresults {input.results} \\
            --fragmentsresults {input.fresults} \\
            --pdf {output.pdf} \\
            --title "{params.condition1}_vs_{params.condition2}__{params.ds}__{params.pt}"
        """