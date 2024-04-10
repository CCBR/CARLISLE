def get_all_peak_files(wildcards):
    files=[]
    if "macs2_narrow" in PEAKTYPE:
        n=expand(join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{treatment_control_list}.{dupstatus}.narrow.peaks.bed"),qthresholds=QTRESHOLDS,treatment_control_list=TREATMENT_LIST_M,dupstatus=DUPSTATUS),
        files.extend(n)
    if "macs2_broad" in PEAKTYPE:
        b=expand(join(RESULTSDIR,"peaks","{qthresholds}","macs2","peak_output","{treatment_control_list}.{dupstatus}.broad.peaks.bed"),qthresholds=QTRESHOLDS,treatment_control_list=TREATMENT_LIST_M,dupstatus=DUPSTATUS),
        files.extend(b)
    if "seacr_stringent" in PEAKTYPE:
        s=expand(join(RESULTSDIR,"peaks","{qthresholds}","seacr","peak_output","{treatment_control_list}.{dupstatus}.stringent.peaks.bed"),qthresholds=QTRESHOLDS,treatment_control_list=TREATMENT_LIST_SG,dupstatus=DUPSTATUS),
        files.extend(s)
    if "seacr_relaxed" in PEAKTYPE:
        r=expand(join(RESULTSDIR,"peaks","{qthresholds}","seacr","peak_output","{treatment_control_list}.{dupstatus}.relaxed.peaks.bed"),qthresholds=QTRESHOLDS,treatment_control_list=TREATMENT_LIST_SG,dupstatus=DUPSTATUS),
        files.extend(r)
    if "gopeaks_narrow" in PEAKTYPE:
        n=expand(join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","peak_output","{treatment_control_list}.{dupstatus}.narrow.peaks.bed"),qthresholds=QTRESHOLDS,treatment_control_list=TREATMENT_LIST_SG,dupstatus=DUPSTATUS),
        files.extend(n)
    if "gopeaks_broad" in PEAKTYPE:
        b=expand(join(RESULTSDIR,"peaks","{qthresholds}","gopeaks","peak_output","{treatment_control_list}.{dupstatus}.broad.peaks.bed"),qthresholds=QTRESHOLDS,treatment_control_list=TREATMENT_LIST_SG,dupstatus=DUPSTATUS),
        files.extend(b)

    files_list=list(itertools.chain.from_iterable(files))
    return files_list

localrules: create_contrast_data_files,make_counts_matrix

rule create_contrast_data_files:
    """
    Output file should include 6 columns
    1)condition     2)sample
    53_H3K4me3	    53_H3K4me3_1 
    
    3)bed
    /results/peaks/gopeaks/53_H3K4me3_1_vs_HN6_IgG_rabbit_negative_control_1.dedup.broadGo_peaks.bed
    
    4)bedgraph                                      5)scaling factor
    /results/bedgraph/53_H3K4me3_1.dedup.bedgraph	86.32596685082872928176	
    
    6)bed
    /results/fragments/53_H3K4me3_1.dedup.fragments.bed
    """
    input:
        replicate_tsv = join(RESULTSDIR,"replicate_sample.tsv"),
        align_stats = rules.gather_alignstats.output.table,
        peaks = get_all_peak_files
    params:
        t_control_list = join(RESULTSDIR,"treatment_control_list.txt"),
        contrast_list="{contrast_list}",
        peak_dir=join(RESULTSDIR,"peaks"),
        bedgraph_dir=join(RESULTSDIR,"bedgraph"),
        frag_bed_path=join(RESULTSDIR,"fragments"),
        qthresholds="{qthresholds}",
        dupstatus="{dupstatus}",
        peak_caller_type="{peak_caller_type}",
        control_flag = config["macs2_control"],
    output:
        contrast_data=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.txt")
    shell:
        """
        # pull conditions
        condition1=`echo {params.contrast_list} | awk -F"_vs_" '{{print $1}}'`
        condition2=`echo {params.contrast_list} | awk -F"_vs_" '{{print $2}}'`

        # set replicates per condition
        replicates1=`cat {input.replicate_tsv} | grep "${{condition1}}" | awk '{{print $1}}'`
        replicates2=`cat {input.replicate_tsv} | grep "${{condition2}}" | awk '{{print $1}}'`

        # identify peak callers
        peak_caller=`echo {params.peak_caller_type} | cut -d"_" -f1`
        peak_type=`echo {params.peak_caller_type} | awk -F"_" '{{print $2}}'`

        # set up output
        if [[ -f {output.contrast_data} ]]; then rm {output.contrast_data}; fi
        touch {output.contrast_data}

        # for each replicate in contrast 1
        for rep in ${{replicates1[@]}}; do
            # set variables
            if [[ {params.control_flag} == "N" ]] && [[ ${{peak_caller}} == "macs2" ]]; then
                treatment_control="${{rep}}_vs_nocontrol"
                treatment=`echo ${{treatment_control}} | awk -F"_vs_" '{{print $1}}'`
            else
                treatment_control=`cat {params.t_control_list} | grep "${{rep}}" | awk \'{{print $1\"_vs_\"$2}}\'`
                treatment=`echo ${{treatment_control}} | awk -F"_vs_" '{{print $1}}'`
            fi
            peakfile="{params.peak_dir}/{params.qthresholds}/${{peak_caller}}/peak_output/${{treatment_control}}.{params.dupstatus}.${{peak_type}}.peaks.bed"
            bedgraphfile="{params.bedgraph_dir}/${{rep}}.{params.dupstatus}.bedgraph"
            sf=`cat {input.align_stats} | grep "${{treatment}}" | awk '{{print $5}}'`
            fragmentbed="{params.frag_bed_path}/${{treatment}}.{params.dupstatus}.fragments.bed"

            # create output file needed for counts matrix
            echo -ne "${{condition1}}\\t${{rep}}\\t${{peakfile}}\\t${{bedgraphfile}}\\t${{sf}}\\t${{fragmentbed}}\\n" >> {output.contrast_data}
        done

        # for each replicate in contrast 2
        for rep in ${{replicates2[@]}}; do
            # set variables
            if [[ {params.control_flag} == "N" ]] && [[ ${{peak_caller}} == "macs2" ]]; then
                treatment_control="${{rep}}_vs_nocontrol"
                treatment=`echo ${{treatment_control}} | awk -F"_vs_" '{{print $1}}'`
            else
                treatment_control=`cat {params.t_control_list} | grep "${{rep}}" | awk \'{{print $1\"_vs_\"$2}}\'`
                treatment=`echo ${{treatment_control}} | awk -F"_vs_" '{{print $1}}'`
            fi
            peakfile="{params.peak_dir}/{params.qthresholds}/${{peak_caller}}/peak_output/${{treatment_control}}.{params.dupstatus}.${{peak_type}}.peaks.bed"
            bedgraphfile="{params.bedgraph_dir}/${{rep}}.{params.dupstatus}.bedgraph"
            sf=`cat {input.align_stats} | grep "${{treatment}}" | awk '{{print $5}}'`
            fragmentbed="{params.frag_bed_path}/${{treatment}}.{params.dupstatus}.fragments.bed"

            # create output file needed for counts matrix
            echo -ne "${{condition2}}\\t${{rep}}\\t${{peakfile}}\\t${{bedgraphfile}}\\t${{sf}}\\t${{fragmentbed}}\\n" >> {output.contrast_data}
        done
        """

rule make_counts_matrix:
    """
    """
    input:
        contrast_data=rules.create_contrast_data_files.output.contrast_data
    output:
        cm=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.countsmatrix.csv"),
        fcm=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.fragmentscountsmatrix.csv"),
        si=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.sampleinfo.csv"),
    params:
        pyscript=join(SCRIPTSDIR,"_make_counts_matrix.py"),
    envmodules:
        TOOLS["python3"],
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
        python {params.pyscript} --bedbedgraph {input.contrast_data} --tmpdir $TMPDIR --countsmatrix {output.cm} --fragmentscountsmatrix {output.fcm} --sampleinfo {output.si}
        """

rule DESeq:
    input:
        contrast_data=rules.create_contrast_data_files.output.contrast_data,
        cm_auc=rules.make_counts_matrix.output.cm,
        cm_frag=rules.make_counts_matrix.output.fcm,
        si=rules.make_counts_matrix.output.si
    output:
        results_auc=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.AUCbased_diffresults.csv"),
        html_auc=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.AUCbased_diffanalysis.html"),
        elbowlimits_auc=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.AUCbased_diffanalysis_elbowlimits.yaml"),
        results_frag=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.fragmentsbased_diffresults.csv"),
        html_frag=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.fragmentsbased_diffanalysis.html"),
        elbowlimits_frag=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.fragmentsbased_diffanalysis_elbowlimits.yaml"),
        pdf=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.venn.pdf")
    params:
        rmd=join(SCRIPTSDIR,"_diff_markdown.Rmd"),
        carlisle_functions=join(SCRIPTSDIR,"_carlisle_functions.R"),
        Rlib_dir=config["Rlib_dir"],
        Rpkg_config=config["Rpkg_config"],
        rscript_diff=join(SCRIPTSDIR,"_diff_markdown_wrapper.R"),
        rscript_venn=join(SCRIPTSDIR,"_plot_results_venn.R"),
        contrast_list="{contrast_list}",
        dupstatus = "{dupstatus}",
        peak_caller_type = "{peak_caller_type}",
        fdr_cutoff = FDRCUTOFF,
        log2fc_cutoff = LFCCUTOFF,
        spiked = NORM_METHOD,
        species = config["genome"],
        gtf=config["reference"][config["genome"]]["gtf"]
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
        TMPDIR_AUC=$TMPDIR/AUC
        TMPDIR_FRAG=$TMPDIR/FRAG
        TMPDIR_VENN=$TMPDIR/VENN

        # pull conditions
        condition1=`echo {params.contrast_list} | awk -F"_vs_" '{{print $1}}'`
        condition2=`echo {params.contrast_list} | awk -F"_vs_" '{{print $2}}'`

        ## Run AUC method
        mkdir -p ${{TMPDIR_AUC}}
        mkdir -p ${{TMPDIR_AUC}}/intermediates_dir
        mkdir -p ${{TMPDIR_AUC}}/knit_root_dir
        cd $TMPDIR_AUC

        Rscript {params.rscript_diff} \\
            --rmd {params.rmd} \\
            --carlisle_functions {params.carlisle_functions} \\
            --Rlib_dir {params.Rlib_dir} \\
            --Rpkg_config {params.Rpkg_config} \\
            --countsmatrix {input.cm_auc} \\
            --sampleinfo {input.si} \\
            --dupstatus {params.dupstatus} \\
            --condition1 ${{condition1}} \\
            --condition2 ${{condition2}} \\
            --fdr_cutoff {params.fdr_cutoff} \\
            --log2fc_cutoff {params.log2fc_cutoff} \\
            --results {output.results_auc} \\
            --report {output.html_auc} \\
            --elbowlimits {output.elbowlimits_auc} \\
            --spiked {params.spiked} \\
            --rawcountsprescaled \\
            --scalesfbymean \\
            --contrast_data {input.contrast_data} \\
            --tmpdir $TMPDIR_AUC \\
            --species {params.species} \\
            --gtf {params.gtf}

        # change elbow limits to provided log2fc if limit is set to .na.real
        sed -i "s/low_limit: .na.real/low_limit: -{params.log2fc_cutoff}/" {output.elbowlimits_auc}
        sed -i "s/up_limit: .na.real/up_limit: {params.log2fc_cutoff}/g" {output.elbowlimits_auc}

        ## Run FRAG method
        mkdir -p ${{TMPDIR_FRAG}}
        mkdir -p ${{TMPDIR_FRAG}}/intermediates_dir
        mkdir -p ${{TMPDIR_FRAG}}/knit_root_dir
        cd $TMPDIR_FRAG

        # Do not use --rawcountsprescaled as these counts are not prescaled!
        Rscript {params.rscript_diff} \\
            --rmd {params.rmd} \\
            --carlisle_functions {params.carlisle_functions} \\
            --Rlib_dir {params.Rlib_dir} \\
            --Rpkg_config {params.Rpkg_config} \\
            --countsmatrix {input.cm_frag} \\
            --sampleinfo {input.si} \\
            --dupstatus {params.dupstatus} \\
            --condition1 ${{condition1}} \\
            --condition2 ${{condition2}} \\
            --fdr_cutoff {params.fdr_cutoff} \\
            --log2fc_cutoff {params.log2fc_cutoff} \\
            --results {output.results_frag} \\
            --report {output.html_frag} \\
            --elbowlimits {output.elbowlimits_frag} \\
            --spiked {params.spiked} \\
            --scalesfbymean \\
            --contrast_data {input.contrast_data} \\
            --tmpdir $TMPDIR_FRAG \\
            --species {params.species} \\
            --gtf {params.gtf}

        # change elbow limits to provided log2fc if limit is set to .na.real
        sed -i "s/low_limit: .na.real/low_limit: -{params.log2fc_cutoff}/" {output.elbowlimits_frag}
        sed -i "s/up_limit: .na.real/up_limit: {params.log2fc_cutoff}/g" {output.elbowlimits_frag}

        ## Run venns
        mkdir -p ${{TMPDIR_VENN}}
        mkdir -p ${{TMPDIR_VENN}}/intermediates_dir
        mkdir -p ${{TMPDIR_VENN}}/knit_root_dir
        cd $TMPDIR_VENN
        Rscript {params.rscript_venn} \\
            --aucresults {output.results_auc} \\
            --fragmentsresults {output.results_frag} \\
            --pdf {output.pdf} \\
            --title "${{condition1}}_vs_${{condition2}}__{params.dupstatus}__{params.peak_caller_type}"
    fi
    """

rule diffbb:
    input:
        results=rules.DESeq.output.results_auc,
        elbowlimits=rules.DESeq.output.elbowlimits_auc,
        fresults=rules.DESeq.output.results_frag,
        felbowlimits=rules.DESeq.output.elbowlimits_frag,
        genome_len=join(BOWTIE2_INDEX,"genome.len"),
    output:
        bed=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.AUCbased_diffresults.bed"),
        bigbed=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.AUCbased_diffresults.bigbed"),
        fbed=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.fragmentsbased_diffresults.bed"),
        fbigbed=join(RESULTSDIR,"peaks","{qthresholds}","contrasts","{contrast_list}.{dupstatus}","{contrast_list}.{dupstatus}.{peak_caller_type}.fragmentsbased_diffresults.bigbed"),
    params:
        script = join(SCRIPTSDIR,"_make_results_bed.py"),
        fdr=FDRCUTOFF,
        lfc=LFCCUTOFF,
    envmodules:
        TOOLS["ucsc"],
        TOOLS["python3"]
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

        ## add check to ensure that DESEQ2 ran to completion
        ## mainly used in tinytest scenarios, but also used if 
        ## Nsamples/group is 1
        check=`wc -c {input.results} | cut -f1 -d" "`
        if [[ $check == 0 ]]; then
            echo "There is only 1 sample per group - this is not allowed in DESEQ2 and leads to incomplete DIFFBB results"
            touch {output.bed}
            touch {output.bigbed}
            touch {output.fbed}
            touch {output.fbigbed}
        else
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
        fi
        """
