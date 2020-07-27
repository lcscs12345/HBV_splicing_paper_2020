#!/bin/bash

# scoring HBV 5'ss and 3'ss
cd ~/virus/doc/hbv/rnaseq/
mkdir maxent
cd ~/virus/doc/hbv/rnaseq/maxent
for i in A2 B2 C2 D3; do \
  (
  pyfasta split -k9 -n1 -o8 ~/virus/ref/hbv/pgrna/${i}/${i}.pgrna.fa
  mv ~/virus/ref/hbv/pgrna/${i}/${i}.pgrna.split.9mer.8overlap.fa .
  sed -n 'n;p' ${i}.pgrna.split.9mer.8overlap.fa \
  | while read i; do \
    cat <(echo ">pgrna") <(echo ${i}) \
    | maxentscan_score5.pl -
  done > ${i}.score5.txt
  sed -i '$d' ${i}.score5.txt
  pyfasta split -k23 -n1 -o22 ~/virus/ref/hbv/pgrna/${i}/${i}.pgrna.fa
  mv ~/virus/ref/hbv/pgrna/${i}/${i}.pgrna.split.23mer.22overlap.fa .
  sed -n 'n;p' ${i}.pgrna.split.23mer.22overlap.fa \
  | while read i; do \
    cat <(echo ">pgrna") <(echo ${i}) \
    | maxentscan_score3.pl -
  done > ${i}.score3.txt
  sed -i '$d' ${i}.score3.txt
  ) &
done
wait

# adding ss position
for i in A2 B2 C2 D3; do \
  cat -n ${i}.score5.txt \
  | awk 'BEGIN {OFS="\t"} {print $1+2,$2,$3,"5ss"}' \
  | sed '$d' \
  | sed '1i pos\tintron.motif\tmaxent.score\tss' > ${i}.pos.score5.txt

  cat -n ${i}.score3.txt \
  | awk 'BEGIN {OFS="\t"} {print $1+20,$2,$3,"3ss"}' \
  | head -n -4 \
  | sed '1i pos\tintron.motif\tmaxent.score\tss' > ${i}.pos.score3.txt
done

# getting junction read counts
for i in A2 B2 C2 D3; do \
  cat \
  ~/virus/doc/hbv/rnaseq/star_hbv/${i}/SJ.out.tab \
  ~/virus/doc/hbv/rnaseq/star_hbv_rep1/${i}/SJ.out.tab \
  | grep ^${i} \
  | awk 'BEGIN {FS="\t"; OFS="|"} $4==1 && $7>9 && $9>=25 {print $1,$2-1,$3+1,".","+" "\t" $7}' \
  | datamash -s -g1 sum 2 | sed 's/|/\t/g' | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$6,$5}' \
  | sort -k2,2n -k3,3n > ${i}.SJ.out.bed
done

# transforming bed for merging with maxent
for i in A2 B2 C2 D3; do \
  cat \
  <(cut -f2,5 ${i}.SJ.out.bed | sed 's/$/\t5ss/') \
  <(cut -f3,5 ${i}.SJ.out.bed | sed 's/$/\t3ss/') \
  | datamash -s -g1,3 sum 2 \
  | sed "1i ${i}_pos\tss\t${i}_reads" > ${i}.SJ.out.txt
done
sed -i 's/2677/2674/' D3.SJ.out.txt

# filtering and adding gap to SJ.out.tab for ss comparison
cd ~/virus/doc/hbv/rnaseq/maxent/
cat A2.SJ.out.bed > A2.SJ.out.offset.bed
awk 'BEGIN {FS=OFS="\t"} $2<540 {print $1,$2,$3,$4,$5,"+"}' B2.SJ.out.bed > B2.SJ.out.offset.bed
awk 'BEGIN {FS=OFS="\t"} $2>=540 {print $1,$2+6,$3+6,$4,$5,"+"}' B2.SJ.out.bed >> B2.SJ.out.offset.bed
awk 'BEGIN {FS=OFS="\t"} $2<540 {print $1,$2,$3,$4,$5,"+"}' C2.SJ.out.bed > C2.SJ.out.offset.bed
awk 'BEGIN {FS=OFS="\t"} $2>=540 {print $1,$2+6,$3+6,$4,$5,"+"}' C2.SJ.out.bed >> C2.SJ.out.offset.bed
# D3
join -1 1 -2 2 -t$'\t' \
-o 2.1,1.2,2.3,2.4,2.5,2.6 \
<(sort -k1,1 ~/virus/doc/hbv/rnaseq/gffcompare/D3.cor) \
<(sort -k2,2 D3.SJ.out.bed) \
| join -1 1 -2 3 -t$'\t' \
-o 2.1,2.2,1.2,2.4,2.5,2.6 \
<(sort -k1,1 ~/virus/doc/hbv/rnaseq/gffcompare/D3.cor) \
<(sort -k3,3 -) | sort -k2,2n -k3,3n > D3.SJ.out.offset.bed
cat \
<(awk '{print $2 "\t5ss"}' *.SJ.out.offset.bed) \
<(awk '{print $3 "\t3ss"}' *.SJ.out.offset.bed) \
| sort -u -k1,1n \
| sed '1i A2_pos\tss' > SJ.out.offset.txt

# getting correspondence table of positions
for i in B2 C2 D3; do
  mafft \
  --addfull ~/virus/ref/hbv/pgrna/${i}/${i}.pgrna.fa \
  --mapout ~/virus/ref/hbv/pgrna/A2/A2.pgrna.fa
  sed '1,2d;s/, /\t/g' ~/virus/ref/hbv/pgrna/${i}/${i}.pgrna.fa.map \
  | cut -f2- | sed "1i ${i}_pos\tA2_pos" > ${i}.pgrna.map
done



# human 5'ss
cd ~/virus/doc/hbv/rnaseq
for i in A2 B2 C2; do \
  awk 'BEGIN {FS=OFS="\t"} $4==1 && $7>9 && $9>=25 {print $1,$2-4,$2+5,$9,$7,"+"}' ${i}.SJ.out.tab \
  | grep -v ${i} > ${i}.SJ.hg19.ss5.bed
  awk 'BEGIN {FS=OFS="\t"} $4==2 && $7>9 && $9>=25 {print $1,$3-6,$3+3,$9,$7,"-"}' ${i}.SJ.out.tab \
  | grep -v ${i} >> ${i}.SJ.hg19.ss5.bed
done
for y in A2 B2 C2; do \
  paste \
    <(cut -f5 ${y}/SJ.hg19.ss5.bed) \
    <(bedtools getfasta -fi ~/riboseq/ref/hg19/hg19.chromFa.fa -bed ${y}/SJ.hg19.ss5.bed -s -fo stdout -tab) \
  | cut -f1,3 \
  | while read i j; do \
    for x in `seq ${i}`; do \
      echo ${j}
    done
  done \
  | tr a-z A-Z \
  | sed 's/T/U/g' \
  | cat -n \
  | awk '{print ">hg19|"$1 "\n" $2}' > ${y}/hg19.ss5.fa
done

# human 3'ss
for i in A2 B2 C2; do \
  awk 'BEGIN {FS=OFS="\t"} $4==1 && $7>9 && $9>=25 {print $1,$3-20,$3+3,$9,$7,"+"}' ${i}.SJ.out.tab \
  | grep -v ${i} > ${i}.SJ.hg19.ss3.bed
  awk 'BEGIN {FS=OFS="\t"} $4==2 && $7>9 && $9>=25 {print $1,$2-4,$2+19,$9,$7,"-"}' ${i}.SJ.out.tab \
  | grep -v ${i} >> ${i}.SJ.hg19.ss3.bed
done
for y in A2 B2 C2; do \
  paste \
    <(cut -f5 ${y}/SJ.hg19.ss3.bed) \
    <(bedtools getfasta -fi ~/riboseq/ref/hg19/hg19.chromFa.fa -bed ${y}/SJ.hg19.ss3.bed -s -fo stdout -tab) \
  | cut -f1,3 \
  | while read i j; do \
    for x in `seq ${i}`; do \
      echo ${j}
    done
  done \
  | tr a-z A-Z \
  | sed 's/T/U/g' \
  | cat -n \
  | awk '{print ">hg19|"$1 "\n" $2}' > ${y}/hg19.ss3.fa
done



# plotting and comparing sequence logo
# HBV 5'ss
cd ~/virus/doc/hbv/rnaseq/maxent
for y in A2 B2 C2 D3; do \
  paste \
    <(cat ${y}.SJ.out.bed \
    | datamash -s -g 2 sum 5) \
    <(cat ${y}.SJ.out.bed \
    | datamash -s -g 2 sum 5 \
    | awk -v y="$y" 'BEGIN {OFS="\t"} {print y,$1-3,$1+6}' \
    | bedtools getfasta -fi ~/virus/ref/hbv/pgrna/${y}/${y}.pgrna.fa -bed stdin -fo stdout -tab) \
  | cut -f2,4 \
  | while read i j; do \
    for x in `seq ${i}`; do \
      echo ${j}
    done
  done \
  | tr a-z A-Z \
  | cat -n \
  | sed 's/T/U/g' \
  | awk -v y="$y" '{print ">" y "|" $1 "\n" $2}' > ${y}.ss5.fa
done

# HBV 3'ss
for y in A2 B2 C2 D3; do \
  paste \
    <(cat ${y}.SJ.out.bed \
    | datamash -s -g 3 sum 5) \
    <(cat ${y}.SJ.out.bed \
    | datamash -s -g 3 sum 5 \
    | awk -v y="$y" 'BEGIN {OFS="\t"} {print y,$1-21,$1+2}' \
    | bedtools getfasta -fi ~/virus/ref/hbv/pgrna/${y}/${y}.pgrna.fa -bed stdin -fo stdout -tab) \
  | cut -f2,4 \
  | while read i j; do \
    for x in `seq ${i}`; do \
      echo ${j}
    done
  done \
  | tr a-z A-Z \
  | cat -n \
  | sed 's/T/U/g' \
  | awk -v y="$y" '{print ">" y "|" $1 "\n" $2}' > ${y}.ss3.fa
done

# human 5'ss and 3'ss
for i in A2 B2 C2 D3; do \
  weblogo -f ${i}/${i}.ss5.fa -o ${i}/${i}.ss5.pdf -D fasta -A rna -F pdf --units probability
  weblogo -f ${i}/${i}.ss3.fa -o ${i}/${i}.ss3.pdf -D fasta -A rna -F pdf --units probability
done
cat *2/hg19.ss5.fa > hg19.ss5.fa
cat *2/hg19.ss3.fa > hg19.ss3.fa
weblogo -f hg19.ss5.fa -o hg19.ss5.pdf -D fasta -A rna -F pdf --units probability
weblogo -f hg19.ss3.fa -o hg19.ss3.pdf -D fasta -A rna -F pdf --units probability
