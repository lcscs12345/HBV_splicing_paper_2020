#!/bin/bash

# rep 1
for i in A2 B2 C2 D3 pUC57; do \
  cd ~/HBV_splicing_paper_2020/doc/hbv/rnaseq/star_hbv_rep1/${i}
  cat \
  <(samtools view -H Aligned.sortedByCoord.out.rmdup.uniq.bam) \
  ~/HBV_splicing_paper_2020/ref/hg19/gencode.v19.rRNA.interval_list \
  > gencode.v19.rRNA.interval_list
  java -jar /usr/local/bin/picard.jar CollectRnaSeqMetrics \
        I=Aligned.sortedByCoord.out.rmdup.uniq.bam \
        O=picard.RNA_Metrics \
        REF_FLAT=~/HBV_splicing_paper_2020/ref/hg19/gencode.v19.annotation.txt \
        RIBOSOMAL_INTERVALS=gencode.v19.rRNA.interval_list \
        CHART_OUTPUT=picard.RNA_Metrics.pdf \
        STRAND_SPECIFICITY=NONE
done

cd ..
cat \
<(awk 'NR>6 && NR<9' A2/picard.RNA_Metrics | cut -f-7) \
<(awk 'NR==8' B2/picard.RNA_Metrics | cut -f-7) \
<(awk 'NR==8' C2/picard.RNA_Metrics | cut -f-7) \
<(awk 'NR==8' pUC57/picard.RNA_Metrics | cut -f-7) > picard.RNA_Metrics_rep1.txt


# rep 2
for i in A2 B2 C2 D3 pUC57; do \
  cd ~/HBV_splicing_paper_2020/doc/hbv/rnaseq/star_hbv/${i}
  cat \
  <(samtools view -H Aligned.sortedByCoord.out.rmdup.uniq.bam) \
  ~/HBV_splicing_paper_2020/ref/hg19/gencode.v19.rRNA.interval_list \
  > gencode.v19.rRNA.interval_list
  java -jar /usr/local/bin/picard.jar CollectRnaSeqMetrics \
        I=Aligned.sortedByCoord.out.rmdup.uniq.bam \
        O=picard.RNA_Metrics \
        REF_FLAT=~/HBV_splicing_paper_2020/ref/hg19/gencode.v19.annotation.txt \
        RIBOSOMAL_INTERVALS=gencode.v19.rRNA.interval_list \
        CHART_OUTPUT=picard.RNA_Metrics.pdf \
        STRAND_SPECIFICITY=NONE
done

cd ..
cat \
<(awk 'NR>6 && NR<9' A2/picard.RNA_Metrics | cut -f-7) \
<(awk 'NR==8' B2/picard.RNA_Metrics | cut -f-7) \
<(awk 'NR==8' C2/picard.RNA_Metrics | cut -f-7) \
<(awk 'NR==8' pUC57/picard.RNA_Metrics | cut -f-7) > picard.RNA_Metrics_rep2.txt