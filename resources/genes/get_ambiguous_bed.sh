wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
gunzip gap.txt.gz
cut -f 2,3,4 gap.txt | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' > hg38.ambiguous.bed
rm -f gap.txt*

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
gunzip gap.txt.gz
cut -f2,3,4 gap.txt | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' > hg19.ambiguous.bed
rm -f gap.txt*

wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/gap.txt.gz
gunzip gap.txt.gz
cut -f2,3,4 gap.txt | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' > mm10.ambiguous.bed
rm -f gap.txt*

wget http://hgdownload.cse.ucsc.edu/goldenPath/mm39/database/gap.txt.gz
gunzip gap.txt.gz
cut -f2,3,4 gap.txt | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' > mm39.ambiguous.bed
rm -f gap.txt*
