
def get_input_fastqs(wildcards):
    d = dict()
    d["R1"] = replicateName2R1[wildcards.replicate]
    d["R2"] = replicateName2R2[wildcards.replicate]
    return d

# check adapters
check_readaccess(config["adapters"])

rule trim:
# """
# Remove adapters using cutadapt:
# * min read length is 35
# * min avg. sliding window quality score of 10 per 10 bp window is required
# """
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
# """
# Align using bowtie:
# * use --dovetail option via "bowtie2_parameters" in config.yaml. This is recommended for 
# Cut and Run where overlapping R1 and R2 alignments are expected
# BAM is sorted and indexed, and stats are collected using flagstat and idxstats
# """
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
# """
# Raw alignment BAMs are filtered for:
# * duplicates only in spikein regions for no_dedup DUPSTATUS
# * all duplicates removed for all regions for dedup DUPSTATUS
# * alignments which are not proper pairs are removed.
# * alignments with fragment length larger than "fragment_len_filter" from config.yaml are removed.
# """
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
        TOOLS["python37"],
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
            mkdir -p ${{TMPDIR}}/{params.replicate}_{params.dupstatus}_picardtmp

            java -Xmx{params.memg} -jar $PICARDJARPATH/picard.jar MarkDuplicates \\
            --INPUT {input.bam} \\
            --OUTPUT ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.tmp1.bam \\
            --ASSUME_SORT_ORDER coordinate \\
            --TMP_DIR ${{TMPDIR}}/{params.replicate}_{params.dupstatus}_picardtmp \\
            --CREATE_INDEX true \\
            --METRICS_FILE {output.bam}.dupmetrics
            
            python {params.pyscript} --inputBAM ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.tmp1.bam \\
            --outputBAM ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.bam \\
            --fragmentlength {params.fragment_len_filter} \\
            --removemarkedduplicates

        else
            # deduplicate only the spikeins
            genome_regions=$(cut -f1 {input.genome_len} | tr '\\n' ' ')
            samtools view -@{threads} -b {input.bam} $genome_regions | \\
                samtools sort -@{threads} -T $TMPDIR -o ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.genome.bam -    

            spikein_regions=$(cut -f1 {input.spikein_len} | tr '\\n' ' ')
            samtools view -@{threads} -b {input.bam} $spikein_regions | \\
                samtools sort -@{threads} -T $TMPDIR -o ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.spikein.bam 
            
            mkdir -p ${{TMPDIR}}/{params.replicate}_{params.dupstatus}_picardtmp
            java -Xmx{params.memg} -jar $PICARDJARPATH/picard.jar MarkDuplicates \\
            --INPUT ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.spikein.bam \\
            --OUTPUT ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.spikein.tmp1.bam \\
            --ASSUME_SORT_ORDER coordinate \\
            --TMP_DIR ${{TMPDIR}}/{params.replicate}_{params.dupstatus}_picardtmp \\
            --CREATE_INDEX true \\
            --METRICS_FILE {output.bam}.spikeinonly.dupmetrics

            python {params.pyscript} --inputBAM ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.spikein.tmp1.bam \\
            --outputBAM ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.spikein.bam \\
            --fragmentlength 1000000 --removemarkedduplicates

            samtools merge -@{threads} -O BAM -o ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.bam \\
                ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.genome.bam \\
                ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.spikein.bam

        fi

        samtools sort -T ${{TMPDIR}} -@{threads} -o {output.bam} ${{TMPDIR}}/{params.replicate}.{params.dupstatus}.filtered.bam
        samtools index -@{threads} {output.bam}
        samtools flagstat {output.bam} > {output.bamflagstat}
        samtools idxstats {output.bam} > {output.bamidxstats}
        """

rule alignstats:
# """
# Number of reads stats collected in per_replicate YAML file:
# * nreads in fastq
# * nreads aligned (to genome and spikein)
# * nreads aligned (to genome and spikein) after fragment_len_filter filtering no_dedup DUPSTATUS
# * nreads aligned (to genome and spikein) after fragment_len_filter and duplicate filtering dedup DUPSTATUS
# """
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
        TOOLS["python37"],
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

localrules: gather_alignstats
rule gather_alignstats:
    input:
        stats=expand(join(RESULTSDIR,"alignment_stats","{replicate}.alignment_stats.yaml"),replicate=REPLICATES),
    output:
        join(RESULTSDIR,"alignment_stats","alignment_stats.tsv")
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
        --outTable {output}
        """

localrules: create_library_norm_scales
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

def get_library_input(wildcards):
    if (NORM_METHOD=="LIBRARY"):
        stats_file=join(RESULTSDIR,"alignment_stats","library_scale.tsv")
    else:
        stats_file=join(RESULTSDIR,"alignment_stats","alignment_stats.tsv")
    return(stats_file)

rule bam2bg:
# """
# Converted filtered BAM files to bedgraph and bigwig formats. SEACR needs bedgraph files as input.
# sf = Constant / [ Nreads aligning to spikein (deduped)] where Constant is defined as "spikein_scale" in config.yaml
# The above sf (scaling factor) is used to scale the bedgraph file. Scaled bedgraph is then converted to bigwig.
# """
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