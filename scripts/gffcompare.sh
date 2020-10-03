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
join -1 1 -2 4 -t$'\t' \
<(awk 'NR>1 {print $2 "\t" $3}' pgrna.pos.txt | sort -k1,1) \
<(sort -k4,4 D3.gtf) \
| cut -f2- | join -1 1 -2 5 -t$'\t' \
-o 2.2,2.3,2.4,2.1,1.2,2.6,2.7,2.8,2.9 \
<(awk 'NR>1 {print $2 "\t" $3}' pgrna.pos.txt | sort -k1,1) \
<(sort -k5,5 -) > D3.gtf2

join -1 1 -2 4 -t$'\t' \
<(awk 'NR>1 {print $2 "\t" $3}' pgrna.pos.txt | sort -k1,1) \
<(sort -k4,4 ~/virus/doc/hbv/ERP013934/X02496.1.pgrna.gtf) \
| cut -f2- | join -1 1 -2 5 -t$'\t' \
-o 2.2,2.3,2.4,2.1,1.2,2.6,2.7,2.8,2.9 \
<(awk 'NR>1 {print $2 "\t" $3}' pgrna.pos.txt | sort -k1,1) \
<(sort -k5,5 -) | sed 's/X02496\.1/HBV/' > X02496.gtf

echo -e 'A2.gtf\nB2.gtf2\nC2.gtf2\nD3.gtf2\nX02496.gtf' > gtf.ls
gffcompare -i gtf.ls -o gffcompare

join -t$'\t'  \
<(sort map.txt) \
<(gtfToGenePred gffcompare.combined.gtf stdout | sort) \
| cut -f2- | sort -k1,1g \
| genePredToGtf file stdin combined.gtf


# calculating spliced transcripts fraction
for i in 1 2 3 4 5; do \
  paste \
  <(cut -f5- gffcompare.tracking | awk -v n=${i} '{print $n}' | cut -f2 -d"|") \
  <(cut -f1 gffcompare.tracking)
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
| awk 'BEGIN {OFS="\t"} {print $1,$4,$3,$5,$2}' > sp.tpm

join -1 1 -2 2 -t$'\t' \
-o 2.1,1.2,2.2,2.3,2.4,2.5 \
<(sort map.txt) <(sed '1d' sp.tpm | sort -k2,2) \
| sed '1i genotype\tsplice_variant\tTCONS\toId\tTPM\ttotal_TPM' > splice_variants.txt
sed -i 's/\tSP/\tKnown\tSP/;s/\tpSP/\tCCS\tpSP/;s/\tsplice_variant/\tstatus\tsplice_variant/;s/CCS\tpSP13/Novel\tpSP13/;s/CCS\tpSP14/Novel\tpSP14/' splice_variants.txt

awk 'BEGIN{OFS="\t"} {if($2~/Known/) {print $1,$3,$6} else {print $1,$2,$6}}' splice_variants.txt > pie.txt

paste \
<(datamash -s -H -g 1 sum 6 < splice_variants.txt | sed '1d') <(cut -f2 total.tpm) \
| awk 'BEGIN{OFS="\t"; print "Genotype\tSpliced\tUnspliced"} {print $1,$2,$3-$2}' > bar.txt
# it's fine to represent these TPM values as percentages
# in fact these TPM values should be divided by 2 (average TPM of rep1 and rep2).
