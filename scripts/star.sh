#!/bin/bash

# generating star index
for i in A2 B2 C2 D3; do \
 cd ~/virus/ref/hbv/pgrna/${i}/;
 mkdir star_index;
 STAR --runMode genomeGenerate \
 --runThreadN 20 \
 --genomeDir star_index \
 --genomeFastaFiles ~/riboseq/ref/hg19/hg19.chromFa.fa ${i}.pgrna.fa \
 --genomeSAindexNbases 5 \
 --genomeChrBinNbits 11 \
 --sjdbGTFfile ~/riboseq/ref/hg19/gencode.v19.annotation.gtf \
 --sjdbOverhang 100
done


# aligning reads
cd ~/virus/doc/hbv/rnaseq/

# sample
# rep 1
mkdir -p star_hbv_rep1
for i in A2 B2 C2 D3; do \
 cd ~/virus/doc/hbv/rnaseq/star_hbv_rep1
 mkdir -p ${i}
 cd ${i}
 STAR --runThreadN 20 \
 --runMode alignReads \
 --twopassMode Basic \
 --twopass1readsN -1 \
 --genomeDir ~/virus/ref/hbv/pgrna/${i}/star_index \
 --readFilesIn /Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_rep1_out/${i}_rep1_1.fastq.gz,/Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_rep1_out/${i}_rep2_1.fastq.gz /Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_rep1_out/${i}_rep1_2.fastq.gz,/Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_rep1_out/${i}_rep2_2.fastq.gz \
 --readFilesCommand zcat \
 --outFilterType BySJout \
 --alignIntronMin 20 \
 --alignIntronMax 1000000 \
 --alignMatesGapMax 1000000 \
 --alignSJoverhangMin 8 \
 --alignSJDBoverhangMin 1 \
 --chimSegmentMin 15 \
 --chimJunctionOverhangMin 15 \
 --sjdbOverhang 100 \
 --sjdbScore 1 \
 --outSAMstrandField intronMotif \
 --outSAMtype BAM SortedByCoordinate \
 --limitBAMsortRAM 2815365604;

 cd ~/virus/doc/hbv/rnaseq/star_hbv_rep1/${i}
 samtools rmdup Aligned.sortedByCoord.out.bam Aligned.sortedByCoord.out.rmdup.bam
 samtools index Aligned.sortedByCoord.out.rmdup.bam
 samtools view -q 255 -bh Aligned.sortedByCoord.out.rmdup.bam > Aligned.sortedByCoord.out.rmdup.uniq.bam
 samtools index Aligned.sortedByCoord.out.rmdup.uniq.bam
 samtools view -h Aligned.sortedByCoord.out.rmdup.uniq.bam ${i} | sed "s/${i}/HBV/" | samtools view -S -bh - > ${i}.bam
 samtools index ${i}.bam
 stringtie \
 Aligned.sortedByCoord.out.rmdup.uniq.bam \
 --rf -G ~/virus/ref/hbv/pgrna/${i}/${i}.pgrna.gtf \
 -p 20 -f 0 -l ${i} -j 2 \
 -A gene_abund.tab > stringtie.gtf
 grep ^${i} stringtie.gtf | sed "s/${i}/HBV/" | awk '$7~/+/' > ${i}.gtf
done

# rep 2
mkdir -p star_hbv
for i in A2 B2 C2 D3; do \
 cd ~/virus/doc/hbv/rnaseq/star_hbv
 mkdir -p ${i}
 cd ${i}
 STAR --runThreadN 20 \
 --runMode alignReads \
 --twopassMode Basic \
 --twopass1readsN -1 \
 --genomeDir ~/virus/ref/hbv/pgrna/${i}/star_index \
 --readFilesIn /Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_out/${i}_rep2_1.fastq.gz /Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_out/${i}_rep2_2.fastq.gz \
 --readFilesCommand zcat \
 --outFilterType BySJout \
 --alignIntronMin 20 \
 --alignIntronMax 1000000 \
 --alignMatesGapMax 1000000 \
 --alignSJoverhangMin 8 \
 --alignSJDBoverhangMin 1 \
 --chimSegmentMin 15 \
 --chimJunctionOverhangMin 15 \
 --sjdbOverhang 100 \
 --sjdbScore 1 \
 --outSAMstrandField intronMotif \
 --outSAMtype BAM SortedByCoordinate \
 --limitBAMsortRAM 2815365604;

 cd ~/virus/doc/hbv/rnaseq/star_hbv/${i}
 samtools rmdup Aligned.sortedByCoord.out.bam Aligned.sortedByCoord.out.rmdup.bam
 samtools index Aligned.sortedByCoord.out.rmdup.bam
 samtools view -q 255 -bh Aligned.sortedByCoord.out.rmdup.bam > Aligned.sortedByCoord.out.rmdup.uniq.bam
 samtools index Aligned.sortedByCoord.out.rmdup.uniq.bam
 samtools view -h Aligned.sortedByCoord.out.rmdup.uniq.bam ${i} | sed "s/${i}/HBV/" | samtools view -S -bh - > ${i}.bam
 samtools index ${i}.bam
 stringtie \
 Aligned.sortedByCoord.out.rmdup.uniq.bam \
 --rf -G ~/virus/ref/hbv/pgrna/${i}/${i}.pgrna.gtf \
 -p 20 -f 0 -l ${i} -j 2 \
 -A gene_abund.tab > stringtie.gtf
 grep ^${i} stringtie.gtf | sed "s/${i}/HBV/" | awk '$7~/+/' > ${i}.gtf
done


# control
# pUC57 rep 1
cd ~/virus/doc/hbv/rnaseq/star_hbv_rep1
mkdir -p pUC57
cd pUC57
STAR --runThreadN 20 \
--runMode alignReads \
--twopassMode Basic \
--twopass1readsN -1 \
--genomeDir ~/riboseq/ref/hg19/star_rnaseq_index \
--readFilesIn /Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_rep1_out/pUC57_rep1_1.fastq.gz,/Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_rep1_out/pUC57_rep2_1.fastq.gz /Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_rep1_out/pUC57_rep1_2.fastq.gz,/Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_rep1_out/pUC57_rep2_2.fastq.gz \
--readFilesCommand zcat \
--outFilterType BySJout \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--sjdbOverhang 100 \
--sjdbScore 1 \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 2815365604;

samtools rmdup Aligned.sortedByCoord.out.bam Aligned.sortedByCoord.out.rmdup.bam
samtools index Aligned.sortedByCoord.out.rmdup.bam
samtools view -q 255 -bh Aligned.sortedByCoord.out.rmdup.bam > Aligned.sortedByCoord.out.rmdup.uniq.bam
samtools view -c -f 0x2 Aligned.sortedByCoord.out.rmdup.uniq.bam


# pUC57 rep 2
cd ~/virus/doc/hbv/rnaseq/star_hbv
mkdir -p pUC57
cd pUC57
STAR --runThreadN 20 \
--runMode alignReads \
--twopassMode Basic \
--twopass1readsN -1 \
--genomeDir ~/riboseq/ref/hg19/star_rnaseq_index \
--readFilesIn /Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_out/pUC57_rep2_1.fastq.gz /Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/skewer_out/pUC57_rep2_2.fastq.gz \
--readFilesCommand zcat \
--outFilterType BySJout \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--sjdbOverhang 100 \
--sjdbScore 1 \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 2815365604;

samtools rmdup Aligned.sortedByCoord.out.bam Aligned.sortedByCoord.out.rmdup.bam
samtools index Aligned.sortedByCoord.out.rmdup.bam
samtools view -q 255 -bh Aligned.sortedByCoord.out.rmdup.bam > Aligned.sortedByCoord.out.rmdup.uniq.bam
samtools view -c -f 0x2 Aligned.sortedByCoord.out.rmdup.uniq.bam
