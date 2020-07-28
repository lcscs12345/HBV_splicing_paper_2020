#!/bin/bash

# gene counting
for i in A2 B2 C2 D3 pUC57; do \
  cd ~/virus/doc/hbv/rnaseq/star_hbv/${i}
  ~/riboseq/src/mmquant/mmquant_static \
  -a ~/riboseq/ref/hg19/gencode.v19.annotation.gtf \
  -r Aligned.sortedByCoord.out.bam \
  -s RF -e Y \
  -O ~/virus/doc/hbv/rnaseq/mmquant/${i}_rep2.mmquant.stat \
  -t 20 \
  -o ~/virus/doc/hbv/rnaseq/mmquant/${i}_rep2.mmquant.tsv
  sed '1d' ~/virus/doc/hbv/rnaseq/mmquant/${i}_rep1.mmquant.tsv \
  | sed "1i gene\t${i}_rep1" > ~/virus/doc/hbv/rnaseq/mmquant/${i}_rep1.txt 

  cd ~/virus/doc/hbv/rnaseq/star_hbv_rep1/${i} 
  ~/riboseq/src/mmquant/mmquant_static \
  -a ~/riboseq/ref/hg19/gencode.v19.annotation.gtf \
  -r Aligned.sortedByCoord.out.bam \
  -s RF -e Y \
  -O ~/virus/doc/hbv/rnaseq/mmquant/${i}_rep1.mmquant.stat \
  -t 20 \
  -o ~/virus/doc/hbv/rnaseq/mmquant/${i}_rep1.mmquant.tsv
  sed '1d' ${i}_rep2.mmquant.tsv \
  | sed "1i gene\t${i}_rep2" > ~/virus/doc/hbv/rnaseq/mmquant/${i}_rep2.txt
done


# merging counts data using the following R commands
setwd("~/virus/doc/hbv/rnaseq/mmquant/")
filenames <- gsub("\\.txt$","", list.files(pattern="\\.txt$"))
for(i in filenames){
  assign(i, read.table(paste(i, ".txt", sep=""), header=T))
}
d <- Reduce(function(...) merge(..., all=TRUE), 
            list(A2_rep1,A2_rep2,B2_rep1,B2_rep2,C2_rep1,C2_rep2,D3_rep1,D3_rep2,pUC57_rep1,pUC57_rep2))
write.table(d, file = "counts.mmquant.txt", sep = "\t", col.names = T, row.names = F, quote=FALSE)


# preparing input file for DESeq2
grep -v "\-\-" counts.mmquant.txt | grep -v NA | awk '$10>99 && $11>99' > counts.txt
join -j1 -t$'\t' \
gene.txt \
<(sed '1d' counts.txt | sed 's/\.[0-9]\+//' | sort -k1,1) \
| cat <(head -n1 counts.txt | sed 's/gene/id\tsymbol\tdescription/') - \
> counts.gene.txt