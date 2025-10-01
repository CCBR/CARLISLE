
def get_input_fastqs(wildcards):
    d = dict()
    d["R1"] = replicateName2R1[wildcards.replicate]
    d["R2"] = replicateName2R2[wildcards.replicate]
    return d

def get_library_input(wildcards):
    if (NORM_METHOD=="LIBRARY"):
        stats_file=join(RESULTSDIR,"alignment_stats","library_scale.tsv")
    else:
        stats_file=join(RESULTSDIR,"alignment_stats","alignment_stats.tsv")
    return(stats_file)


# check adapters
check_readaccess(config["adapters"])

localrules: gather_alignstats, create_library_norm_scales

rule trim:
    """
    Remove adapters using cutadapt:
    * min read length is 35
    * min avg. sliding window quality score of 10 per 10 bp window is required
    """
    input:
        unpack(get_input_fastqs)
    output:
        R1 = temp(join(RESULTSDIR,"trim","{replicate}.R1.trim.fastq.gz")),
        R2 = temp(join(RESULTSDIR,"trim","{replicate}.R2.trim.fastq.gz")),
    params:
        adapters = config["adapters"],
    threads: getthreads("trim")
    envmodules:
        TOOLS["cutadapt"],
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
        cutadapt \\
        --pair-filter=any \\
        --nextseq-trim=2 \\
        --trim-n \\
        -n 5 -O 5 \\
        -q 10,10 \\
        -m 35:35 \\
        -b file:{params.adapters} \\
        -B file:{params.adapters} \\
        -j {threads} \\
        -o {output.R1} \\
        -p {output.R2} \\
        {input.R1} {input.R2}
        """

rule align:
    """
    Align using bowtie:
    * use --dovetail option via "bowtie2_parameters" in config.yaml. This is recommended for
    Cut and Run where overlapping R1 and R2 alignments are expected
    BAM is sorted and indexed, and stats are collected using flagstat and idxstats
    """
    input:
        R1 = rules.trim.output.R1,
        R2 = rules.trim.output.R2,
        bt2 = join(BOWTIE2_INDEX,"ref.1.bt2")
    output:
        bam=temp(join(RESULTSDIR,"bam","raw","{replicate}.bam")),
        bai=temp(join(RESULTSDIR,"bam","raw","{replicate}.bam.bai")),
        bamflagstat=temp(join(RESULTSDIR,"bam","raw","{replicate}.bam.flagstat")),
        bamidxstats=temp(join(RESULTSDIR,"bam","raw","{replicate}.bam.idxstats")),
    params:
        replicate = "{replicate}",
        bowtie2_parameters = config["bowtie2_parameters"],
        bt2_base = join(BOWTIE2_INDEX,"ref"),
        pyscript = join(SCRIPTSDIR,"_filter_bam.py"),
        qfilter = config["mapping_quality"],
    threads: getthreads("align")
    envmodules:
        TOOLS["bowtie2"],
        TOOLS["samtools"],
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
        bowtie2 \\
            -p {threads} \\
            {params.bowtie2_parameters} \\
            -x {params.bt2_base}  \\
            -1 {input.R1} -2 {input.R2}  | \\
            samtools view -f 0x2 -q {params.qfilter} -bS - |  \\
            samtools sort -T ${{TMPDIR}} -@{threads} -o {output.bam}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.bamflagstat}
        samtools idxstats {output.bam} > {output.bamidxstats}
        """

rule filter:
    """
    Raw alignment BAMs are filtered for:
    * no duplicates removed for no_dedup DUPSTATUS
    * all linear duplicates (except spikein regions) removed for all regions for dedup DUPSTATUS
    * alignments which are not proper pairs are removed.
    * alignments with fragment length larger than "fragment_len_filter" from config.yaml are removed.
    """
    input:
        bam=rules.align.output.bam,
        bai=rules.align.output.bai,
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        spikein_len = join(BOWTIE2_INDEX,"spikein.len"),
    output:
        bam=join(RESULTSDIR,"bam","{replicate}.{dupstatus}.bam"),
        bai=join(RESULTSDIR,"bam","{replicate}.{dupstatus}.bam.bai"),
        bamflagstat=join(RESULTSDIR,"bam","{replicate}.{dupstatus}.bam.flagstat"),
        bamidxstats=join(RESULTSDIR,"bam","{replicate}.{dupstatus}.bam.idxstats"),
    params:
        replicate = "{replicate}",
        dupstatus = "{dupstatus}",  # can be "dedup" or "no_dedup"
        fragment_len_filter = config["fragment_len_filter"],
        pyscript = join(SCRIPTSDIR,"_filter_bam.py"),
        memg = getmemg("filter")
    threads: getthreads("filter")
    envmodules:
        TOOLS["bowtie2"],
        TOOLS["samtools"],
        TOOLS["python3"],
        TOOLS["picard"],
        TOOLS["ucsc"]
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
        if [[ "{params.dupstatus}" == "dedup" ]];then

            # deduplicate only the spikeins
            genome_regions=$(cut -f1 {input.genome_len} | tr '\\n' ' ')
            samtools view -@{threads} -b {input.bam} $genome_regions | \\
                samtools sort -@{threads} -T $TMPDIR -o ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.genome.bam -

            spikein_regions=$(cut -f1 {input.spikein_len} | tr '\\n' ' ')
            samtools view -@{threads} -b {input.bam} $spikein_regions | \\
                samtools sort -@{threads} -T $TMPDIR -o ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.spikein.bam

            samtools index ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.genome.bam
            samtools index ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.spikein.bam

            # remove linear duplicates from genome only
            python {params.pyscript} --bam_in ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.genome.bam \\
            --bam_out ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.genome.filtered.bam \\
            --metrics_path ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.genome.metrics.txt \\
            --fraglen {params.fragment_len_filter} \\
            --minmapq 20

            # only filter spikein bam
            python {params.pyscript} --bam_in ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.spikein.bam \\
            --bam_out ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.spikein.filtered.bam \\
            --metrics_path ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.spikein.metrics.txt \\
            --fraglen {params.fragment_len_filter} \\
            --minmapq 20 --no-linear-dedup          

            samtools merge -@{threads} -O BAM -o ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.bam \\
                ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.genome.filtered.bam \\
                ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.spikein.filtered.bam
        else
            # no_dedup
            python {params.pyscript} --bam_in {input.bam} \\
            --bam_out ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.bam \\
            --metrics_path ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.metrics.txt \\
            --fraglen {params.fragment_len_filter} \\
            --minmapq 20 --no-linear-dedup
        fi

        samtools sort -T ${{TMPDIR}} -@{threads} -o {output.bam} ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.bam
        samtools index -@{threads} {output.bam}
        samtools flagstat {output.bam} > {output.bamflagstat}
        samtools idxstats {output.bam} > {output.bamidxstats}
        """

rule alignstats:
    """
    Number of reads stats collected in per_replicate YAML file:
    * nreads in fastq
    * nreads aligned (to genome and spikein)
    * nreads aligned (to genome and spikein) after fragment_len_filter filtering no_dedup DUPSTATUS
    * nreads aligned (to genome and spikein) after fragment_len_filter and duplicate filtering dedup DUPSTATUS
    """
    input:
        R1 = rules.trim.output.R1,
        raw_alignment_idxstats = rules.align.output.bamidxstats,
        filtered_no_dedup_alignment_idxstats = join(RESULTSDIR,"bam","{replicate}.no_dedup.bam.idxstats"),
        filtered_dedup_alignment_idxstats = join(RESULTSDIR,"bam","{replicate}.dedup.bam.idxstats"),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        spikein_len = join(BOWTIE2_INDEX,"spikein.len"),
    output:
        outyaml = join(RESULTSDIR,"alignment_stats","{replicate}.alignment_stats.yaml")
    params:
        replicate = "{replicate}",
        pyscript = join(SCRIPTSDIR,"_get_nreads_stats.py"),
    threads: getthreads("alignstats")
    envmodules:
        TOOLS["python3"],
    shell:
        """
        set -exo pipefail

        # _get_nreads_stats parameters:
        # @param1 = R1 fastq file
        # @param2 = spike-in len file
        # @param3 = genome len file
        # @param4 = raw idx stats file
        # @param5 = filtered idx stats file - no_dedup
        # @param6 = filtered idx stats file - dedup
        # @param7 = out yaml file

        python {params.pyscript} \\
            {input.R1} \\
            {input.spikein_len} \\
            {input.genome_len} \\
            {input.raw_alignment_idxstats} \\
            {input.filtered_no_dedup_alignment_idxstats} \\
            {input.filtered_dedup_alignment_idxstats} \\
            {output.outyaml}
        """

rule gather_alignstats:
    input:
        stats=expand(join(RESULTSDIR,"alignment_stats","{replicate}.alignment_stats.yaml"),replicate=REPLICATES),
    output:
        table=join(RESULTSDIR,"alignment_stats","alignment_stats.tsv")
    params:
        rscript = join(SCRIPTSDIR,"_make_alignment_stats_table.R"),
        spikein_scale = config["spikein_scale"],
    envmodules:
        TOOLS["R"]
    shell:
        """
        set -exo pipefail
        file1=$(echo {input.stats} | awk '{{print $1}}')
        dir=$(dirname $file1)
        Rscript {params.rscript} \\
        --yamlDir $dir \\
        --excludeFromName ".alignment_stats.yaml" \\
        --scaleConstant {params.spikein_scale} \\
        --outTable {output.table}
        """

rule create_library_norm_scales:
    input:
        stats=expand(join(RESULTSDIR,"alignment_stats","{replicate}.alignment_stats.yaml"),replicate=REPLICATES),
    output:
        scalefile=join(RESULTSDIR,"alignment_stats","library_scale.tsv")
    params:
        rscript = join(SCRIPTSDIR,"_make_library_norm_table.R"),
    envmodules:
        TOOLS["R"]
    shell:
        """
        set -exo pipefail
        file1=$(echo {input.stats} | awk '{{print $1}}')
        dir=$(dirname $file1)
        Rscript {params.rscript} \\
        --yamlDir $dir \\
        --excludeFromName ".alignment_stats.yaml" \\
        --outTable {output.scalefile}
        """

rule bam2bg:
    """
    Converted filtered BAM files to bedgraph and bigwig formats. SEACR needs bedgraph files as input.
    sf = Constant / [ Nreads aligning to spikein (deduped)] where Constant is defined as "spikein_scale" in config.yaml
    The above sf (scaling factor) is used to scale the bedgraph file. Scaled bedgraph is then converted to bigwig.
    """
    input:
        bam = rules.filter.output.bam,
        bai = rules.filter.output.bai,
        bamidxstats = rules.filter.output.bamidxstats,
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
        spikein_len = join(BOWTIE2_INDEX,"spikein.len"),
        library_file = get_library_input
    output:
        fragments_bed = join(RESULTSDIR,"fragments","{replicate}.{dupstatus}.fragments.bed"),
        bg=join(RESULTSDIR,"bedgraph","{replicate}.{dupstatus}.bedgraph"),
        bw=join(RESULTSDIR,"bigwig","{replicate}.{dupstatus}.bigwig"),
        sf_yaml=join(RESULTSDIR,"bedgraph","{replicate}.{dupstatus}.sf.yaml")
    params:
        spikein = NORM_METHOD,
        replicate = "{replicate}",
        dupstatus = "{dupstatus}",
        fragment_len_filter = config["fragment_len_filter"],
        spikein_scale = config["spikein_scale"],
        regions = REGIONS,
        memG = getmemG("bam2bg")
    threads: getthreads("bam2bg")
    envmodules:
        TOOLS["bedtools"],
        TOOLS["samtools"],
        TOOLS["ucsc"]
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

        if [[ "{params.spikein}" == "NONE" ]];then
            echo "No spike-in scale was used"
            spikein_scale=1
        elif [[ "{params.spikein}" == "LIBRARY" ]];then
            spikein_scale=`cat {input.library_file} | grep {params.replicate} | grep {params.dupstatus} | cut -f2 -d" " | head -n1`
            echo "The spike-in is generated from the library size"
        else
            spikein_readcount=$(while read a b;do awk -v a=$a '{{if ($1==a) {{print $3}}}}' {input.bamidxstats};done < {input.spikein_len} | awk '{{sum=sum+$1}}END{{print sum}}')

            # if the spikein_readcount is below threshold, then there is not enough of the spikein control to use
            total_count=$(awk '{{sum+=$3; sum+=$4;}}END{{print sum;}}' {input.bamidxstats})
            spikein_percent=`echo "scale=6 ; $spikein_readcount / $total_count * 100" | bc`;\

            if [[ $spikein_percent < 0.001 ]]; then
                echo "The spikein percentage of {input.bam} was below the threshold (0.001%) at $spikein_percent%. The spikein_scale was set to 1."
                spikein_scale=1
            else
                spikein_scale=$(echo "{params.spikein_scale} / $spikein_readcount" | bc -l)
            fi
        fi

        # create fragments file
        samtools view -b -@{threads} {input.bam} {params.regions} | samtools sort -n -@{threads} -T $TMPDIR -o ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.bam -
        bedtools bamtobed -bedpe -i ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.bam > ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.bed

        awk -v fl={params.fragment_len_filter} '{{ if ($1==$4 && $6-$2 < fl) {{print $0}}}}' ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.bed > ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.clean.bed

        cut -f 1,2,6 ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.clean.bed | \\
            LC_ALL=C sort --buffer-size={params.memG} --parallel={threads} --temporary-directory=$TMPDIR -k1,1 -k2,2n -k3,3n > {output.fragments_bed}

        # run bedtools
        bedtools genomecov -bg -scale $spikein_scale -i {output.fragments_bed} -g {input.genome_len} > {output.bg}

        # create bigwig
        bedGraphToBigWig {output.bg} {input.genome_len} {output.bw}

        # add to YAML
        echo "spikein_scaling_factor=$spikein_scale" > {output.sf_yaml}
        """

rule deeptools_bw:
    input:
        bam = join(RESULTSDIR,"bam","{replicate}.{dupstatus}.bam"),
        genome_len = join(BOWTIE2_INDEX,"genome.len")
    output:
        clean_bam = join(RESULTSDIR,"deeptools","clean","{replicate}.{dupstatus}.clean.bam"),
        clean_bai = join(RESULTSDIR,"deeptools","clean","{replicate}.{dupstatus}.clean.bam.bai"),
        bw = join(RESULTSDIR,"deeptools","clean","{replicate}.{dupstatus}.clean.bigwig"),
    envmodules:
        TOOLS["samtools"],
        TOOLS["deeptools"],
    shell:
        """
        genome_size=`cat {input.genome_len} | awk '{{ sum+=$2 }} END{{ print sum }}'`
        samtools view -b -@4 {input.bam} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | samtools sort -@4 -o {output.clean_bam}
        samtools index {output.clean_bam}
        bamCoverage --bam {output.clean_bam} -o {output.bw} --binSize 25 --smoothLength 75 --numberOfProcessors 32 --normalizeUsing RPGC --effectiveGenomeSize $genome_size --centerReads
        """

rule deeptools_prep:
    input:
        bw = expand(join(RESULTSDIR,"deeptools","clean","{replicate}.{dupstatus}.clean.bigwig"), replicate=REPLICATES,dupstatus=DUPSTATUS),
    output:
        deeptools_prep = temp(expand(join(RESULTSDIR,"deeptools","clean", "{group}.{dupstatus}.deeptools_prep"),group=[a+b for a in TREATMENTS+["all_samples"] for b in ["", ".prot"]],dupstatus=DUPSTATUS)),
    run:
        for dupstatus in DUPSTATUS:
            for i in ["","prot."]:
                for item in TREATMENT_CONTROL_LIST:
                    labels = item.split('_vs_')
                    treatment, control = labels
                    bws = [join(RESULTSDIR,"deeptools","clean","{}.{}.clean.bigwig".format(item, dupstatus)) for item in labels]
                    f=open(join(RESULTSDIR,"deeptools","clean","{}.{}{}.deeptools_prep".format(treatment, i, dupstatus) ),'w')
                    f.write("{}\n".format(dupstatus))
                    f.write("{}\n".format(" ".join(bws)))
                    f.write("{}\n".format(" ".join(labels)))
                    f.close()

        for dupstatus in DUPSTATUS:
            labels = []
            for c in CONTROL_to_TREAT_DICT:
                labels += CONTROL_to_TREAT_DICT[c]
                labels.append(c)
            bws = [join(RESULTSDIR,"deeptools","clean","{}.{}.clean.bigwig".format(item, dupstatus)) for item in labels]

            for i in ["","prot."]:
                f_all=open(join(RESULTSDIR,"deeptools","clean","all_samples.{}{}.deeptools_prep".format(i, dupstatus)),'w')
                f_all.write("all_samples\n")
                f_all.write("{}\n".format(" ".join(bws)))
                f_all.write("{}\n".format(" ".join(labels)))
                f_all.close()

rule deeptools_mat:
    input:
        deeptools_prep = join(RESULTSDIR,"deeptools","clean", "{group}.{dupstatus}.deeptools_prep"),
        bw = expand(join(RESULTSDIR,"deeptools","clean","{replicate}.{dupstatus}.clean.bigwig"), replicate=REPLICATES,dupstatus=DUPSTATUS),
    output:
        metamat=temp(join(RESULTSDIR,"deeptools","clean", "{group}.{dupstatus}.metagene.mat.gz")),
        TSSmat=temp(join(RESULTSDIR,"deeptools","clean","{group}.{dupstatus}.TSS.mat.gz")),
        bed=temp(join(RESULTSDIR,"deeptools","clean","{group}.{dupstatus}.geneinfo.bed")),
    params:
        prebed="/data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_basic/geneinfo.bed",
        pythonver=TOOLS["python3"],
        deeptoolsver=TOOLS["deeptools"],
    threads:
        getthreads("deeptools_mat"),
    run:
        import re
        commoncmd="module load {params.deeptoolsver}; module load {params.pythonver};"
        listfile=list(map(lambda z:z.strip().split(),open(input.deeptools_prep,'r').readlines()))
        ext=listfile[0][0]
        bws=listfile[1]
        labels=listfile[2]
        if "prot" in wildcards.group:
            cmd1="grep --line-buffered 'protein_coding' "+ params.prebed  +" | awk -v OFS='\t' -F'\t' '{{print $1, $2, $3, $5, \".\", $4}}' > "+output.bed
        else:
            cmd1="awk -v OFS='\t' -F'\t' '{{print $1, $2, $3, $5, \".\", $4}}' "+params.prebed+" > "+output.bed
        cmd2="computeMatrix scale-regions -S "+" ".join(bws)+" -R "+output.bed+" -p 16 --upstream 1000 --regionBodyLength 2000 --downstream 1000 --skipZeros -o "+output.metamat+" --samplesLabel "+" ".join(labels)
        cmd3="computeMatrix reference-point -S "+" ".join(bws)+" -R "+output.bed+" -p 16 --referencePoint TSS --upstream 3000 --downstream 3000 --skipZeros -o "+output.TSSmat+" --samplesLabel "+" ".join(labels)
        shell(commoncmd+cmd1)
        shell(commoncmd+cmd2)
        shell(commoncmd+cmd3)

rule deeptools_heatmap:
    input:
        metamat=rules.deeptools_mat.output.metamat,
        TSSmat=rules.deeptools_mat.output.TSSmat,
    output:
        metaheat=join(RESULTSDIR,"deeptools","{group}.{dupstatus}.metagene_heatmap.pdf"),
        TSSheat=join(RESULTSDIR,"deeptools","{group}.{dupstatus}.TSS_heatmap.pdf"),
        metaline=join(RESULTSDIR,"deeptools","{group}.{dupstatus}.metagene_profile.pdf"),
        TSSline=join(RESULTSDIR,"deeptools","{group}.{dupstatus}.TSS_profile.pdf"),
    envmodules:
        TOOLS["deeptools"],
    shell:
        """
        plotHeatmap -m {input.metamat} -out {output.metaheat} --colorMap 'PuOr_r' --zMin auto --zMax auto --yAxisLabel 'average RPGC' --regionsLabel 'genes' --legendLocation 'none'
        plotHeatmap -m {input.TSSmat} -out {output.TSSheat} --colorMap 'RdBu_r' --zMin auto --zMax auto --yAxisLabel 'average RPGC' --regionsLabel 'genes' --legendLocation 'none'
        plotProfile -m {input.metamat} -out {output.metaline} --plotHeight 15 --plotWidth 15 --perGroup --yAxisLabel 'average RPGC' --plotType 'se' --legendLocation upper-right
        plotProfile -m {input.TSSmat} -out {output.TSSline} --plotHeight 15 --plotWidth 15 --perGroup --yAxisLabel 'average RPGC' --plotType 'se' --legendLocation upper-left
        """

rule cov_correlation:
    """
    Create replicate correlation plots from filtered BAM files
    """
    input:
        bams=expand(join(RESULTSDIR,"bam","{replicate}.{{dupstatus}}.bam"),replicate=REPLICATES),
        align_table=join(RESULTSDIR,"alignment_stats","alignment_stats.tsv")
    output:
        counts=join(RESULTSDIR,"deeptools","all.{dupstatus}.readCounts.npz"),
        pearson_corr=join(RESULTSDIR,"deeptools","all.{dupstatus}.PearsonCorr.tab"),
        pearson_plot=join(RESULTSDIR,"deeptools","all.{dupstatus}.PearsonCorr.png"),
        pca=join(RESULTSDIR,"deeptools","all.{dupstatus}.PCA.tab"),
        hc=join(RESULTSDIR,"deeptools","all.{dupstatus}.Pearson_heatmap.png"),
        pca_format=join(RESULTSDIR,"deeptools","all.{dupstatus}.PearsonPCA.png")
    params:
        rscript=join(SCRIPTSDIR,"_plot_correlation.R"),
        dupstatus="{dupstatus}"
    container: config['containers']['carlisle_r']
    threads: getthreads("cov_correlation")
    shell:
        """
        # Calculate genome-wide coverage
        multiBamSummary bins \
         --bamfiles {input.bams} \
         --smartLabels \
         -out {output.counts} \
         -p {threads}

        # Plot heatmap - Pearson
        plotCorrelation \
         -in {output.counts} \
         --corMethod pearson --skipZeros \
         --whatToPlot heatmap --plotNumbers \
         -o {output.pearson_plot} \
         --removeOutliers \
         --outFileCorMatrix {output.pearson_corr}

        # Plot PCA
        plotPCA -in {output.counts}  --outFileNameData {output.pca}

        # Plot heatmap and PCA (formatted)
        Rscript {params.rscript} {output.pearson_corr} {output.pca} {input.align_table} {params.dupstatus} {output.hc} {output.pca_format}
        """
