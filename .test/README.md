################################################
# CHANGELOG
### samples retired 2/6/24
<!-- VS589.R1.fastq.gz,HN6_H3K4me3_2.R1.fastq.gz
VS589.R2.fastq.gz,HN6_H3K4me3_2.R2.fastq.gz
VS594.R1.fastq.gz,HN6_H4K20me3_1.R1.fastq.gz
VS594.R2.fastq.gz,HN6_H4K20me3_1.R2.fastq.gz
VS595.R1.fastq.gz,HN6_H4K20me3_2.R1.fastq.gz
VS595.R2.fastq.gz,HN6_H4K20me3_2.R2.fastq.gz
VS597.R1.fastq.gz,53_H4K20m3_1.R1.fastq.gz
VS597.R2.fastq.gz,53_H4K20m3_1.R2.fastq.gz
VS598.R1.fastq.gz,53_H4K20m3_2.R1.fastq.gz -->
### 12/16/22
New test data is added to pipeline

################################################
# HG38
# sample location
/data/CCBR/projects/ccbr1155/CS031014/fastq

# sample ids mapping; original to renamed
VS586.R1.fastq.gz,HN6_IgG_rabbit_negative_control_1.R1.fastq.gz
VS586.R2.fastq.gz,HN6_IgG_rabbit_negative_control_1.R2.fastq.gz
VS591.R1.fastq.gz,53_H3K4me3_1.R1.fastq.gz
VS591.R2.fastq.gz,53_H3K4me3_1.R2.fastq.gz
VS592.R1.fastq.gz,53_H3K4me3_2.R1.fastq.gz
VS592.R2.fastq.gz,53_H3K4me3_2.R2.fastq.gz
VS588.R1.fastq.gz,HN6_H3K4me3_1.R1.fastq.gz
VS588.R2.fastq.gz,HN6_H3K4me3_1.R2.fastq.gz

# subset each of the samples 
fq_original_loc="/data/CCBR/projects/ccbr1155/CS031014/fastq";\
test_loc="/data/CCBR_Pipeliner/Pipelines/CARLISLE_dev/.test";\
fq_list=("HN6_IgG_rabbit_negative_control_1.R1.fastq.gz" "HN6_IgG_rabbit_negative_control_1.R2.fastq.gz" "HN6_H3K4me3_1.R1.fastq.gz" "HN6_H3K4me3_1.R2.fastq.gz" "HN6_H3K4me3_2.R1.fastq.gz" "HN6_H3K4me3_2.R2.fastq.gz" "53_H3K4me3_1.R1.fastq.gz" "53_H3K4me3_1.R2.fastq.gz" "53_H3K4me3_2.R1.fastq.gz" "53_H3K4me3_2.R2.fastq.gz" "HN6_H4K20me3_1.R1.fastq.gz" "HN6_H4K20me3_1.R2.fastq.gz" "HN6_H4K20me3_2.R1.fastq.gz" "HN6_H4K20me3_2.R2.fastq.gz" "53_H4K20m3_1.R1.fastq.gz" "53_H4K20m3_1.R2.fastq.gz" "53_H4K20m3_2.R1.fastq.gz" "53_H4K20m3_2.R2.fastq.gz");\
for fq in ${fq_list[@]}; do echo $fq; remove_ext=`echo $fq | sed "s/.gz//g"`; zcat $fq_original_loc/$fq | head -1000000 > $test_loc/$remove_ext; gzip $test_loc/$remove_ext; done

################################################
# MM10

fq_original_loc="/data/sevillas2/carlisle/mm10";\
test_loc="/data/CCBR_Pipeliner/Pipelines/CARLISLE/v2.4.1/.test";\
fq_list=("RECQ1_MEF_WT_Unt1_S1_R1_001.fastq.gz" "RECQ1_MEF_WT_Unt1_S1_R2_001.fastq.gz" "RECQ1_MEF_WT_Unt2_S2_R1_001.fastq.gz" "RECQ1_MEF_WT_Unt2_S2_R2_001.fastq.gz" "RECQ1_MEF_KO_Unt5_S5_R1_001.fastq.gz" "RECQ1_MEF_KO_Unt5_S5_R2_001.fastq.gz" "IgG-1_S1_L001_R1_001.fastq.gz" "IgG-1_S1_L001_R2_001.fastq.gz");\
for fq in ${fq_list[@]}; do echo $fq; remove_ext=`echo $fq |sed "s/_S[0-9]_R/.R/g" | sed "s/_001.fastq.gz.gz/.gz/g" | sed "s/_/./g"`; zcat $fq_original_loc/$fq | head -1000000 > $test_loc/$remove_ext; gzip $test_loc/$remove_ext; done