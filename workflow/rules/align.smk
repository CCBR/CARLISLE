
def get_input_fastqs(wildcards):
    d = dict()
    d["R1"] = replicateName2R1[wildcards.replicate]
    d["R2"] = replicateName2R2[wildcards.replicate]
    return d

# check adapters
check_readaccess(config["adapters"])

rule trim:
    input:
        unpack(get_input_fastqs)
    output:
        R1 = join(RESULTSDIR,"tmp","trim","{replicate}.R1.trim.fastq.gz"),
        R2 = join(RESULTSDIR,"tmp","trim","{replicate}.R2.trim.fastq.gz"),
    params:
        adapters = config["adapters"],
    threads: getthreads("trim")
    envmodules:
        TOOLS["cutadapt"],
    shell:"""
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
    input:
        R1 = join(RESULTSDIR,"tmp","trim","{replicate}.R1.trim.fastq.gz"),
        R2 = join(RESULTSDIR,"tmp","trim","{replicate}.R2.trim.fastq.gz"),
        bt2 = join(BOWTIE2_INDEX,"ref.1.bt2")
    output:
        bam=join(RESULTSDIR,"tmp","bam","{replicate}.bam"),
        bai=join(RESULTSDIR,"tmp","bam","{replicate}.bam.bai"),
        bamflagstat=join(RESULTSDIR,"tmp","bam","{replicate}.bam.flagstat"),
        bamidxstats=join(RESULTSDIR,"tmp","bam","{replicate}.bam.idxstats"),
    params:
        replicate = "{replicate}",
        bowtie2_parameters = config["bowtie2_parameters"],
        bt2_base = join(BOWTIE2_INDEX,"ref"),
        pyscript = join(SCRIPTSDIR,"_filter_bam.py")
    threads: getthreads("align")
    envmodules:
        TOOLS["bowtie2"],
        TOOLS["samtools"],
    shell:"""
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
    samtools view -bS - |  \\
    samtools sort -T ${{TMPDIR}} -@{threads} -o {output.bam}
samtools index {output.bam}
samtools flagstat {output.bam} > {output.bamflagstat}
samtools idxstats {output.bam} > {output.bamidxstats}
"""

rule filter:
    input:
        bam=join(RESULTSDIR,"tmp","bam","{replicate}.bam"),
        bai=join(RESULTSDIR,"tmp","bam","{replicate}.bam.bai")
    output:
        bam=join(RESULTSDIR,"bam","{replicate}.bam"),
        bai=join(RESULTSDIR,"bam","{replicate}.bam.bai"),
        bamflagstat=join(RESULTSDIR,"bam","{replicate}.bam.flagstat"),
        bamidxstats=join(RESULTSDIR,"bam","{replicate}.bam.idxstats"),
    params:
        replicate = "{replicate}",
        bowtie2_parameters = config["bowtie2_parameters"],
        bt2_base = join(BOWTIE2_INDEX,"ref"),
        pyscript = join(SCRIPTSDIR,"_filter_bam.py")
    threads: getthreads("filter")
    envmodules:
        TOOLS["bowtie2"],
        TOOLS["samtools"],
        TOOLS["python37"]
    shell:"""
set -exo pipefail
if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then 
    TMPDIR="/lscratch/$SLURM_JOB_ID"
else
    dirname=$(basename $(mktemp))
    TMPDIR="/dev/shm/$dirname"
    mkdir -p $TMPDIR
fi
outbam_bn=$(basename {output.bam})
python {params.pyscript} --inputBAM {input.bam}\\
    --outputBAM ${{TMPDIR}}/${{outbam_bn%.*}}.filtered.bam
samtools sort -T ${{TMPDIR}} -@{threads} -o {output.bam} ${{TMPDIR}}/${{outbam_bn%.*}}.filtered.bam
samtools flagstat {output.bam} > {output.bamflagstat}
samtools idxstats {output.bam} > {output.bamidxstats}
"""
