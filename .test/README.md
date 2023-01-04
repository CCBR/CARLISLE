# 12/16/22
New test data is added to pipeline

# sample location
/data/CCBR/projects/ccbr1155/CS031014/fastq

# sample ids mapping; original to renamed
VS586.R1.fastq.gz,HN6_IgG_rabbit_negative_control_1.R1.fastq.gz
VS586.R2.fastq.gz,HN6_IgG_rabbit_negative_control_1.R2.fastq.gz
VS588.R1.fastq.gz,HN6_H3K4me3_1.R1.fastq.gz
VS588.R2.fastq.gz,HN6_H3K4me3_1.R2.fastq.gz
VS589.R1.fastq.gz,HN6_H3K4me3_2.R1.fastq.gz
VS589.R2.fastq.gz,HN6_H3K4me3_2.R2.fastq.gz
VS591.R1.fastq.gz,53_H3K4me3_1.R1.fastq.gz
VS591.R2.fastq.gz,53_H3K4me3_1.R2.fastq.gz
VS592.R1.fastq.gz,53_H3K4me3_2.R1.fastq.gz
VS592.R2.fastq.gz,53_H3K4me3_2.R2.fastq.gz
VS594.R1.fastq.gz,HN6_H4K20me3_1.R1.fastq.gz
VS594.R2.fastq.gz,HN6_H4K20me3_1.R2.fastq.gz
VS595.R1.fastq.gz,HN6_H4K20me3_2.R1.fastq.gz
VS595.R2.fastq.gz,HN6_H4K20me3_2.R2.fastq.gz
VS597.R1.fastq.gz,53_H4K20m3_1.R1.fastq.gz
VS597.R2.fastq.gz,53_H4K20m3_1.R2.fastq.gz
VS598.R1.fastq.gz,53_H4K20m3_2.R1.fastq.gz

# subset each of the samples 
fq_original_loc="/data/CCBR/projects/ccbr1155/CS031014/fastq";\
test_loc="/data/CCBR_Pipeliner/Pipelines/CARLISLE_dev/.test";\
fq_list=("HN6_IgG_rabbit_negative_control_1.R1.fastq.gz" "HN6_IgG_rabbit_negative_control_1.R2.fastq.gz" "HN6_H3K4me3_1.R1.fastq.gz" "HN6_H3K4me3_1.R2.fastq.gz" "HN6_H3K4me3_2.R1.fastq.gz" "HN6_H3K4me3_2.R2.fastq.gz" "53_H3K4me3_1.R1.fastq.gz" "53_H3K4me3_1.R2.fastq.gz" "53_H3K4me3_2.R1.fastq.gz" "53_H3K4me3_2.R2.fastq.gz" "HN6_H4K20me3_1.R1.fastq.gz" "HN6_H4K20me3_1.R2.fastq.gz" "HN6_H4K20me3_2.R1.fastq.gz" "HN6_H4K20me3_2.R2.fastq.gz" "53_H4K20m3_1.R1.fastq.gz" "53_H4K20m3_1.R2.fastq.gz" "53_H4K20m3_2.R1.fastq.gz" "53_H4K20m3_2.R2.fastq.gz");\
for fq in ${fq_list[@]}; do echo $fq; remove_ext=`echo $fq | sed "s/.gz//g"`; zcat $fq_original_loc/$fq | head -1000000 > $test_loc/$remove_ext; gzip $test_loc/$remove_ext; done