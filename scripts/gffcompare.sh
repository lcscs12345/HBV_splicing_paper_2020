#!/bin/bash

# comparing the gtf files of rep 1 and 2 using rep 1 as reference
# and merging BAM file replicates for Gviz's sashimi plot
mkdir ~/virus/doc/hbv/rnaseq/gffcompare/
for i in A2 B2 C2 D3; do \
  cd ~/virus/doc/hbv/rnaseq/gffcompare/
  gffcompare \
  -o ${i} \
  -r ~/virus/doc/hbv/rnaseq/star_hbv_rep1/${i}/${i}.gtf \
  ~/virus/doc/hbv/rnaseq/star_hbv/${i}/${i}.gtf
  join -j1 -t$'\t' \
  <(grep "=" ${i}.tracking | cut -f5 | awk 'BEGIN {FS="|"} {print $2}' | sort) \
  <(gtfToGenePred -genePredExt ${i}.annotated.gtf stdout | awk '$8!=1' | sort -k1,1) \
  > ${i}.txt
  cat ${i}.txt | genePredToGtf file stdin ${i}.gtf
  
  samtools merge ${i}.bam \
  ~/virus/doc/hbv/rnaseq/star_hbv_rep1/${i}/${i}.bam \
  ~/virus/doc/hbv/rnaseq/star_hbv/${i}/${i}.bam
  samtools index ${i}.bam
done


# inserting gaps
for i in B2 C2; do \
  awk 'BEGIN {FS=OFS="\t"} $4>=540 {print $1,$2,$3,$4+6,$5+6,$6,$7,$8,$9}' ${i}.gtf > ${i}.gtf2
  awk 'BEGIN {FS=OFS="\t"} $4==1 && $5>=540 {print $1,$2,$3,$4,$5+6,$6,$7,$8,$9}' ${i}.gtf >> ${i}.gtf2
  awk 'BEGIN {FS=OFS="\t"} $4==458 {print $1,$2,$3,$4,$5+6,$6,$7,$8,$9}' ${i}.gtf >> ${i}.gtf2
  awk '$4<540 && $4!=1 && $4!=458' ${i}.gtf >> ${i}.gtf2
  awk '$4==1 && $5<540' ${i}.gtf >> ${i}.gtf2
done


# merging and comparing gtf
echo -e 'A2.gtf\nB2.gtf2\nC2.gtf2' > gtf.ls
echo -e 'A2.gtf\nB2.gtf2\nC2.gtf2\tD3.gtf2' > gtf.ls2

gtfmerge union gtf.ls union.gtf
gffcompare -i gtf.ls -o compare

join -1 1 -2 4 -t$'\t' \
<(sort -k1,1 D3.cor) \
<(sort -k4,4 D3.gtf) \
| cut -f2- | join -1 1 -2 5 -t$'\t' \
-o 2.2,2.3,2.4,2.1,1.2,2.6,2.7,2.8,2.9 \
<(sort -k1,1 D3.cor) \
<(sort -k5,5 -) > D3.gtf2

gffcompare -i gtf.ls2 -o gffcompare


# calculating spliced transcripts fraction
cd ~/virus/doc/hbv/rnaseq/gffcompare
for i in 1 2 3 4; do \
  paste \
  <(cut -f5- gffcompare.tracking | awk -v n=${i} '{print $n}' | cut -f2 -d"|") \
  <(cut -f1 gffcompare.tracking) \
  | grep -v TCONS_00000012
done | grep -v "\-" > gffcompare.map

for i in A2 B2 C2 D3; do \
  join -j1 -t$'\t' \
  <(sed 's/|/\t/g' ${i}.tracking | grep "=" | cut -f4,7 | sort -k1,1) \
  <(grep TPM ~/virus/doc/hbv/rnaseq/star_hbv_rep1/${i}/${i}.gtf | awk 'BEGIN {OFS="\t"} {print $12,$NF}' | sed 's/[";]//g' | sort -k1,1) \
  | join -1 2 -2 1 -t$'\t' \
  <(sort -k2,2 -) \
  <(grep TPM ~/virus/doc/hbv/rnaseq/star_hbv/${i}/${i}.gtf | awk 'BEGIN {OFS="\t"} {print $12,$NF}' | sed 's/[";]//g' | sort -k1,1) \
  | awk -v s=${i} 'BEGIN {OFS="\t"} {print $1, s, $3+$4}'
done > hbv.tpm

for i in A2 B2 C2 D3; do \
  grep TPM ~/virus/doc/hbv/rnaseq/star_hbv_rep1/${i}/${i}.gtf | awk -v s=${i} '{print s "\t" $NF}' | sed 's/[";]//g' | datamash -g1 sum 2
  grep TPM ~/virus/doc/hbv/rnaseq/star_hbv/${i}/${i}.gtf | awk -v s=${i} '{print s "\t" $NF}' | sed 's/[";]//g' | datamash -g1 sum 2
done | datamash -g1 sum 2 > total.tpm

join -j1 -t$'\t' \
<(sort -k1,1 gffcompare.map) \
<(sort -k1,1 hbv.tpm) \
| join -1 1 -2 3 -t$'\t' \
total.tpm \
<(sort -k3,3 -) \
| awk 'BEGIN {OFS="\t"} {print $1,$4,$3,$5,$2}' \
| sed '1i genotype\tcons_sp\tsp\ttpm\ttotal_tpm' > sp.tpm
mv sp.tpm hbv.tpm
# it's fine to represent these TPM values as percentages
# in fact these TPM values should be divided by 2 (average TPM of rep1 and rep2).
