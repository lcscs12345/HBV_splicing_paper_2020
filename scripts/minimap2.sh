#!/bin/bash

# generating a gtf annotation file for pgRNAs using PacBio CCS reads
# the fastq files are available at https://www.ebi.ac.uk/ena/browser/view/ERP013934

# aligning CCS reads to HBV genotype D pgRNA
cd ~/virus/doc/hbv/ERP013934/
minimap2 \
-a -xsplice -C5 -O6,24 -B4 -uf \
-L --cs=long \
X02496.1.fa \
m141122_094445_42243_c100695392550000001823149003241534_s1_p0.ccs.fastq \
T1.ccs.fastq \
T2.ccs.fastq \
T5.ccs.fastq \
| samtools view -bS - \
| samtools sort -o ccs.bam
samtools index ccs.bam

# generating a junction file
conda activate 2passpipeline
2passtools score \
-j 4 -d 10 -m 'GTAG' \
-w 128 -k 6 -lt 0.1 -ht 0.6 --stranded \
-s 12345 \
-f X02496.1.fa -o score.bed ccs.bam
awk '$4~/GTAG/ && $6~/+/ && $5>10' score.bed | cut -f-6 > score.plus.bed
samtools view -hS ccs.bam > ccs.sam
conda deactivate

# generating a 'reference' gtf file
conda activate flair_env
python ~/virus/src/flair/bin/sam_to_psl.py ccs.sam ccs.psl
python ~/virus/src/flair/bin/psl_to_bed.py ccs.psl ccs.bed
python ~/virus/src/flair/flair.py correct \
-q ccs.bed \
-g X02496.1.fa \
-j score.plus.bed
python ~/virus/src/flair/flair.py collapse \
-q flair_all_corrected.bed \
--filter comprehensive \
--keep_intermediate \
-g X02496.1.fa \
-r m141122_094445_42243_c100695392550000001823149003241534_s1_p0.ccs.fastq \
T1.ccs.fastq T2.ccs.fastq T5.ccs.fastq
conda deactivate

# patching the start and end positions
bedToGenePred flair.collapse.isoforms.bed stdout \
| sed 's|/|\t|;s|/|_|' | cut -f2- \
| awk 'BEGIN{OFS="\t"} NR>3 && NR<36 {print $1,$2,$3,"0","3207","3207","3207",$8,$9,$10}' \
| sed 's/\t6,/\t0,/;s/\t207,/\t0,/;s/\t623,/\t0,/;s/,[0-9]\+,$/,3207,/' \
| awk '!seen[$9,$10]++' | sort -k8,8g \
| genePredToGtf file stdin stdout > X02496.1.pgrna.gtf

# mapping D3 pgRNA coordinates to standard coordinates (centered at the EcoRI site)
mafft \
--addfull ~/virus/ref/hbv/pgrna/D3/D3.pgrna.fa \
--mapout ~/virus/ref/hbv/pgrna/C2/C2.pgrna.fa
sed '1,2d;s/, /\t/g' ~/virus/ref/hbv/pgrna/D3/D3.pgrna.fa.map \
| cut -f2- | sed "1i D3_pos\tC2_pos" > ~/virus/doc/hbv/rnaseq/maxent/pgrna.map
awk '$3~/exon/ {print $4 "\n" $5}' X02496.1.pgrna.gtf | sort -u \
| join -t$'\t' \
- <(sed '1d' ~/virus/doc/hbv/rnaseq/maxent/pgrna.map | sort -k1,1) \
| join -1 2 -2 1 -t$'\t' \
- <(sed '1d' ~/virus/doc/hbv/rnaseq/maxent/C2.pgrna.map | sort -k1,1) \
| awk 'BEGIN{OFS="\t"} {if($1+1813<3215) {print $0,$1+1813} else if ($1+1813>3215) {print $0,$1+1813-3215}}' \
| sed '1i C2\tD3\tA2\tlabel' > ~/virus/doc/hbv/rnaseq/gffcompare/pgrna.pos.txt

# creating genotype specific gtf files
awk 'NR>2 {print "-e \x027s|\\t" $2 "\\t|\\t" $3 "\\t|\x027 \\"}' ~/virus/doc/hbv/rnaseq/gffcompare/pgrna.pos.txt \
| sed '1i sed \\' | sed '$ s|.$|X02496\.1\.pgrna\.gtf > ~/virus/ref/hbv/pgrna/A2/A2.pgrna.gtf|' | sh
sed -i s'/X02496\.1/A2/' ~/virus/ref/hbv/pgrna/A2/A2.pgrna.gtf
awk 'NR>2 {print "-e \x027s|\\t" $2 "\\t|\\t" $1 "\\t|\x027 \\"}' ~/virus/doc/hbv/rnaseq/gffcompare/pgrna.pos.txt \
| sed '1i sed \\' | sed '$ s|.$|X02496\.1\.pgrna\.gtf > ~/virus/ref/hbv/pgrna/B2/B2.pgrna.gtf|' | sh
sed -i s'/X02496\.1/B2/' ~/virus/ref/hbv/pgrna/B2/B2.pgrna.gtf
awk 'NR>2 {print "-e \x027s|\\t" $2 "\\t|\\t" $1 "\\t|\x027 \\"}' ~/virus/doc/hbv/rnaseq/gffcompare/pgrna.pos.txt \
| sed '1i sed \\' | sed '$ s|.$|X02496\.1\.pgrna\.gtf > ~/virus/ref/hbv/pgrna/C2/C2.pgrna.gtf|' | sh
sed -i s'/X02496\.1/C2/' ~/virus/ref/hbv/pgrna/C2/C2.pgrna.gtf
sed s'/X02496\.1/D3/' X02496.1.pgrna.gtf > ~/virus/ref/hbv/pgrna/D3/D3.pgrna.gtf
