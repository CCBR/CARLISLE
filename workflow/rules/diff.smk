localrules:contrast_init

rule contrast_init:
    input:
        expand(join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.narrowPeak"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        expand(join(RESULTSDIR,"peaks","macs2","{replicate}","{replicate}.{dupstatus}_peaks.broadPeak"),replicate=REPLICATES,dupstatus=DUPSTATUS),
        expand([join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.dedup.norm.stringent.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        expand([join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.dedup.norm.relaxed.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        expand([join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.no_dedup.norm.stringent.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        expand([join(RESULTSDIR,"peaks","seacr","{treatment}_vs_{control}","{treatment}_vs_{control}.no_dedup.norm.relaxed.bed")],zip,treatment=TREATMENTS,control=CONTROLS),
        join(RESULTSDIR,"replicate_sample.tsv")
    output:
        outtsv=join(RESULTSDIR,"peaks","contrasts","bed_bedgraph_paths.tsv"),
    params:
        resultsdir = RESULTSDIR,
    shell:"""
while read replicate sample;do
    for dupstatus in dedup no_dedup;do
        bedgraph=$(find {params.resultsdir} -iname "*${{replicate}}.${{dupstatus}}.bedgraph")
        for peaktype in narrowPeak broadPeak norm.relaxed.bed norm.stringent.bed;do
            if [[ $dupstatus == "dedup" ]];then
                bed=$(find {params.resultsdir} -iname "*${{replicate}}*${{dupstatus}}*${{peaktype}}" |grep -v no_dedup)
            else
                bed=$(find {params.resultsdir} -iname "*${{replicate}}*${{dupstatus}}*${{peaktype}}")
            fi
            echo -ne "$replicate\\t$sample\\t$dupstatus\\t$peaktype\\t$bed\\t$bedgraph\\n"
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
    shell:"""
for condition in {params.condition1} {params.condition2}
do
    while read s g d p bed bg;do
        if [[ "$g" == "$condition" ]];then
            if [[ "$d" == "{params.ds}" ]];then
                if [[ "$p" == "{params.pt}" ]];then
                    echo -ne "$g\t$s\t$bed\t$bg\n"
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
        si=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_sampleinfo.txt"),
    params:
        pyscript=join(SCRIPTSDIR,"_make_counts_matrix.py"),
        condition1 = "{c1}",
        condition2 = "{c2}",
        ds = "{ds}",
        pt = "{pt}",
    envmodules:
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
python {params.pyscript} --bedbedgraph {input.inputs} --tmpdir $TMPDIR --countsmatrix {output.cm} --sampleinfo {output.si}
"""

rule DESeq:
    input:
        cm=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_countsmatrix.txt"),
        si=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_sampleinfo.txt"),
    output:
        results=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_diffresults.txt"),
        html=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_diffanalysis.html"),
    params:
        rscript=join(SCRIPTSDIR,"_diff_markdown_wrapper.R"),
        rmd=join(SCRIPTSDIR,"_diff_markdown.Rmd"),
        condition1 = "{c1}",
        condition2 = "{c2}",
        ds = "{ds}",
        pt = "{pt}",
    envmodules:
        TOOLS["R"]
    shell:"""
set -exo pipefail
if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then 
    TMPDIR="/lscratch/$SLURM_JOB_ID"
else
    dirname=$(basename $(mktemp))
    TMPDIR="/dev/shm/$dirname"
    mkdir -p $TMPDIR
fi
Rscript {params.rscript} \\
    --rmd {params.rmd} \\
    --countsmatrix {input.cm} \\
    --sampleinfo {input.si} \\
    --dupstatus {params.ds} \\
    --condition1 {params.condition1} \\
    --condition2 {params.condition2} \\
    --results {output.results} \\
    --report {output.html}
"""

rule diffbb:
    input:
        results=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_diffresults.txt"),
        genome_len = join(BOWTIE2_INDEX,"genome.len"),
    output:
        bed=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_diffresults.bed"),
        bigbed=join(RESULTSDIR,"peaks","contrasts","{c1}_vs_{c2}__{ds}__{pt}","{c1}_vs_{c2}__{ds}__{pt}_diffresults.bigbed")
    params:
        fdr=FDRCUTOFF,
        lfc=LFCCUTOFF,
        randstr = RANDOMSTR,
    envmodules:
        TOOLS["ucsc"]
    shell:"""
set -exo pipefail
if [[ -d "/lscratch/$SLURM_JOB_ID" ]]; then 
    TMPDIR="/lscratch/$SLURM_JOB_ID"
else
    dirname=$(basename $(mktemp))
    TMPDIR="/dev/shm/$dirname"
    mkdir -p $TMPDIR
fi
awk -F"\\t" -v OFS="\\t" -v fdr={params.fdr} -v lfc={params.lfc} '{{if ($7<fdr && $3>lfc){{print $8,$9,$10,".","0","+","0","0","255,0,0"}}}}' {input.results} > ${{TMPDIR}}/{params.randstr}.tmp_up
awk -F"\\t" -v OFS="\\t" -v fdr={params.fdr} -v lfc={params.lfc} '{{if ($7<fdr && $3< -1*lfc){{print $8,$9,$10,".","0","-","0","0","0,255,0"}}}}' {input.results} > ${{TMPDIR}}/{params.randstr}.tmp_down
cat ${{TMPDIR}}/{params.randstr}.tmp_up ${{TMPDIR}}/{params.randstr}.tmp_down | sort -S 4G -T /dev/shm -k1,1 -k2,2n > {output.bed}
rm -f ${{TMPDIR}}/{params.randstr}.tmp_up ${{TMPDIR}}/{params.randstr}.tmp_down
bedToBigBed -type=bed9 {output.bed} {input.genome_len} {output.bigbed}
"""
