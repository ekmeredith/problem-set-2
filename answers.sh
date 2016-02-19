#! /usr/bin/env bash

# Question 1 Use BEDtools intersect to identify the size of the largest overlap between
# CTCF and H3K4me3 locations.

tfbs="$HOME/Desktop/class/data-sets/bed/encode.tfbs.chr22.bed"
h3="$HOME/Desktop/class/data-sets/bed/encode.h3k4me3.hela.chr22.bed"
answer_1=$(cat $tfbs \
    |awk '($4 =="CTCF")' \
    |bedtools intersect -a - -b $h3 -wo \
    |awk '{print $NF}' \
    | sort -k1nr \
    | head -n1)
echo "answer-1: $answer_1"

# Question 2 Use BEDtools to calculate the GC content of nucleotides 19,000,000 to
#19,000,500 on chr22 of `hg19` genome build. Report the GC content
#as a fraction (e.g., 0.50).
genome="$HOME/Desktop/class/data-sets/fasta/hg19.chr22.fa"
echo -e "chr22\t19000000\t19000500" > region.bed
answer_2=$(bedtools nuc -fi $genome -bed region.bed | awk '(NR!=1)'\
|cut -f5)

echo "answer-2: $answer_2"

# Question 3 Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e.,
#interval) that has the largest mean signal in `ctcf.hela.chr22.bg.gz`.

signal="$HOME/Desktop/class/data-sets/bedtools/ctcf.hela.chr22.bg"
tfbs="$HOME/Desktop/class/data-sets/bed/encode.tfbs.chr22.bed"
cat $tfbs | awk '$4=="CTCF"'|cut -f1,2,3|sortBed -i - > intervals.bed
answer_3=$(bedtools map -a intervals.bed  -b $signal -c 4  -o mean \
    |awk '$4 != "." ' \
    |sort -k4nr \
    |head -n1 \
    |awk '{print $3-$2}')

echo "answer-3: $answer_3"

# Question 4
#Use BEDtools to identify the gene promoter (defined as 1000 bp upstream of
#a TSS) with the highest median signal in `ctcf.hela.chr22.bg.gz`.  Report
#the gene name (e.g., 'ABC123')

signal="$HOME/Desktop/class/data-sets/bedtools/ctcf.hela.chr22.bg"
tss="$HOME/Desktop/class/data-sets/bed/tss.hg19.chr22.bed"
genome="$HOME/Desktop/class/data-sets/genome/hg19.genome"
bedtools flank -i $tss -g $genome -l 1000 -r 0 -s |bedtools sort  -i - > promoters.bed

answer_4=$(bedtools map -a promoters.bed -b $signal -c 4 -o median \
    |sort -k7nr \
    |head -n1 \
    |cut -f4)

echo "answer-4: $answer_4"

# Question 5 Use BEDtools to identify the longest interval on `chr22` that is not
#covered by `genes.hg19.bed.gz`. Report the interval like `chr1:100-500`.

genes="$HOME/Desktop/class/data-sets/bed/genes.hg19.bed"
genome="$HOME/Desktop/class/data-sets/genome/hg19.genome"
bedtools complement -i $genes -g $genome > notgenes.bed
answer_5=$(grep -w "chr22" notgenes.bed \
    |awk '{print $1,$2,$3, $3-$2}' \
    |sort -k4nr \
    |head -n1 \
    |awk 'BEGIN {OFS=""} {print $1, ":",$2,"-",$3}')
echo "answer-5: $answer_5"

# Question 6 (extra credit) Use one or more BEDtools that we haven't covered in class. Be creative
#I want to know if CTCF and CTCFL bind the overlapping intervals. I will use the
#bedtools jaccard to figure out the overlap, and compare this to the
#overlap of CTCF with the transcription factor MYC

tfbs="$HOME/Desktop/class/data-sets/bed/encode.tfbs.chr22.bed"
grep -w "CTCF" $tfbs > ctcf.bed
grep -w "CTCFL" $tfbs > ctcfl.bed
grep -w "MYC" $tfbs >myc.bed
cftl=$(bedtools jaccard -a ctcf.bed -b ctcfl.bed \
    |awk 'NR != 1' \
    |cut -f3)
myc=$(bedtools jaccard -a ctcf.bed -b myc.bed \
    |awk 'NR != 1' \
    |cut -f3)
echo "answer-6: CTCF and CTCFL have a Jaccard similarity score of $cftl.
For comparison,  CTCF and myc have a similarity of $myc. Therefore, CTCFL doesn't seem to bind the
same sites as CTCF"
