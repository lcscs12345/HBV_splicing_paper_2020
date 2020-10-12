#!/bin/bash

# completed Splicing Index (coSI)

# setting up directories
mkdir ~/virus/doc/hbv/rnaseq/ipsa_out ~/virus/doc/hbv/rnaseq/ipsa_out/hbv ~/virus/doc/hbv/rnaseq/ipsa_out/hg19
OUT=/Network/Servers/biocldap.otago.ac.nz/Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/ipsa_out
IPSA=/Network/Servers/biocldap.otago.ac.nz/Volumes/userdata/student_users/chunshenlim/virus/src/ipsa

# HBV
# building indexs for ipsa
cd ${IPSA}
for i in A2 B2 C2 D3; do \
  sed "s/^HBV/${i}/" ~/virus/doc/hbv/rnaseq/gffcompare/${i}.gtf > ${OUT}/hbv/${i}.gtf
  cp ~/virus/ref/hbv/pgrna/${i}/${i}.pgrna.fa ${OUT}/hbv/
done
cat ${OUT}/hbv/*.gtf > ${OUT}/hbv/hbv.gtf
perl Perl/transcript_elements.pl - < ${OUT}/hbv/hbv.gtf > ${OUT}/hbv/hbv.gfx
cd ${OUT}/hbv/
${IPSA}/maptools-3.2/bin/transf -dir ./ -idx hbv.idx -dbx hbv.dbx 

# merging the BAM files of HBV genotypes
for i in A2 B2 C2 D3; do \
  samtools view -h ~/virus/doc/hbv/rnaseq/star_hbv_rep1/${i}/${i}.bam | sed "s/SN:HBV/SN:${i}/;s/\tHBV/\t${i}/" \
  | samtools view -bS - > ${i}_rep1.bam
  samtools view -h ~/virus/doc/hbv/rnaseq/star_hbv/${i}/${i}.bam | sed "s/SN:HBV/SN:${i}/;s/\tHBV/\t${i}/" \
  | samtools view -bS - > ${i}_rep2.bam
done
samtools merge hbv.bam *.bam
rm *_rep*

#running ipsa
${IPSA}/sjcount-3.1/sjcount -bam hbv.bam \
-ssc ssc.1.tsv -log ssj.1.log -ssj ssj.1.tsv \
-nbins 51 -read1 1 -quiet
awk -f ${IPSA}/awk/aggregate.awk -v degree=0 -v readLength=51 -v margin=5 -v prefix= -v logfile=ssc.2.log ssc.1.tsv > ssc.2.tsv  
awk -f ${IPSA}/awk/aggregate.awk -v degree=1 -v readLength=51 -v margin=5 -v prefix= -v logfile=ssj.2.log ssj.1.tsv > ssj.2.tsv  
cd ${IPSA}
perl Perl/annotate.pl -deltaSS 10 \
-annot ${OUT}/hbv/hbv.gfx \
-dbx ${OUT}/hbv/hbv.dbx \
-idx ${OUT}/hbv/hbv.idx \
-in ${OUT}/hbv/ssj.2.tsv \
> ${OUT}/hbv/ssj.3.tsv

cd ${OUT}/hbv/
awk -f ${IPSA}/awk/choose_strand.awk ssj.3.tsv > ssj.4.tsv
awk -f ${IPSA}/awk/constrain_ssc.awk -v jncfile=ssj.4.tsv ssc.2.tsv > ssc.4.tsv
awk '$4>=1.5' ssc.4.tsv > ssc.6.tsv
awk '$4>=1.5 && $5>=0' ssj.4.tsv > ssj.6.tsv

cd ${IPSA}
perl Perl/zeta.pl \
-annot ${OUT}/hbv/hbv.gfx \
-ssc ${OUT}/hbv/ssc.6.tsv \
-ssj ${OUT}/hbv/ssj.6.tsv -mincount 10 \
> ${OUT}/hbv/zeta.gff

cd ${OUT}/hbv/
awk 'BEGIN{OFS="\t"} $3~/exon/ {print $1,$4-1,$5,$10,".",$7}' zeta.gff | sed 's/[";]//g' > exon.zeta.bed
awk 'BEGIN{OFS="\t"} $3~/intron/ {print $1,$4-1,$5,$10,$12,$7}' zeta.gff | sed 's/[";]//g' > intron.zeta.bed


# human
# building indexs for ipsa
cd ${IPSA}
perl Perl/transcript_elements.pl - < ~/riboseq/ref/hg19/gencode.v19.annotation.gtf > ~/riboseq/ref/hg19/gencode.v19.annotation.gfx
# chroms directory contains fasta files of indivdual human chromosome
cd ~/riboseq/ref/hg19/
${IPSA}/maptools-3.2/bin/transf -dir chroms/ -idx hg19.idx -dbx hg19.dbx

# merging pUC57 re1 and re2
cd ${IPSA}/hg19/
samtools merge pUC57.bam \
~/virus/doc/hbv/rnaseq/star_hbv_rep1/pUC57/Aligned.sortedByCoord.out.rmdup.uniq.bam \
~/virus/doc/hbv/rnaseq/star_hbv/pUC57/Aligned.sortedByCoord.out.rmdup.uniq.bam

#running ipsa
${IPSA}/sjcount-3.1/sjcount -bam pUC57.bam \
-ssc ssc.1.tsv -log ssj.1.log -ssj ssj.1.tsv \
-nbins 51 -read1 1 -quiet
awk -f ${IPSA}/awk/aggregate.awk -v degree=0 -v readLength=51 -v margin=5 -v prefix= -v logfile=ssc.2.log ssc.1.tsv > ssc.2.tsv  
awk -f ${IPSA}/awk/aggregate.awk -v degree=1 -v readLength=51 -v margin=5 -v prefix= -v logfile=ssj.2.log ssj.1.tsv > ssj.2.tsv  
cd ${IPSA}
perl Perl/annotate.pl -deltaSS 10 \
-annot ~/riboseq/ref/hg19/gencode.v19.annotation.gfx \
-dbx ~/riboseq/ref/hg19/hg19.dbx \
-idx ~/riboseq/ref/hg19/hg19.idx \
-in ${OUT}/hg19/ssj.2.tsv \
> ${OUT}/hg19/ssj.3.tsv

cd ${OUT}/hg19/
awk -f ${IPSA}/awk/choose_strand.awk ssj.3.tsv > ssj.4.tsv
awk -f ${IPSA}/awk/constrain_ssc.awk -v jncfile=ssj.4.tsv ssc.2.tsv > ssc.4.tsv
awk '$4>=1.5' ssc.4.tsv > ssc.6.tsv
awk '$4>=1.5 && $5>=0' ssj.4.tsv > ssj.6.tsv

cd ${IPSA}
perl Perl/zeta.pl \
-annot ~/riboseq/ref/hg19/gencode.v19.annotation.gfx \
-ssc ${OUT}/hg19/ssc.6.tsv \
-ssj ${OUT}/hg19/ssj.6.tsv -mincount 10 \
> ${OUT}/hg19/zeta.gff

cd ${OUT}/hg19/
awk 'BEGIN{OFS="\t"} $3~/exon/ {print $1,$4-1,$5,$10,".",$7}' zeta.gff | sed 's/[";]//g' > exon.zeta.pUC57.bed
awk 'BEGIN{OFS="\t"} $3~/intron/ {print $1,$4-1,$5,$10,$12,$7}' zeta.gff | sed 's/[";]//g' > intron.zeta.pUC57.bed
