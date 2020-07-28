library(tidyverse)
library(DESeq2)
setwd("~/virus/doc/hbv/rnaseq/mmquant/")
countData <- read.table("counts.gene.txt", header=TRUE)
countData <- column_to_rownames(countData, var = "id")
countData <- subset(countData, select=-c(symbol,description))

colData <- data.frame(
  row.names=colnames(countData),
  condition=c("A2","A2","B2","B2","C2","C2","D3","D3","0pUC57","0pUC57"),
  replicate=c("rep1","rep2","rep1","rep2","rep1","rep2","rep1","rep2","rep1","rep2"))
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ replicate + condition)

vsd <- vst(dds)
pcaData <- plotPCA(vsd, intgroup=c("condition","replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("pca.pdf", width=4, height=3)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$replicate)
pcaData <- plotPCA(vsd, intgroup=c("condition","replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("pca.limma.pdf", width=4, height=3)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

dds <- DESeq(dds)

# A2
A2_pUC57 <- results(dds, contrast=c("condition","A2","0pUC57"), pAdjustMethod="fdr")
A2_pUC57 <- A2_pUC57[order(A2_pUC57$padj),]
write.table(as.data.frame(A2_pUC57), file="A2_pUC57.tsv", sep = "\t", col.names = T, row.names = T, quote=FALSE)
pdf("A2_pUC57.maplot.pdf", width = 3.5, height = 4)
plotMA(A2_pUC57, main="A2 vs pUC57", alpha=0.05, ylim=c(-1,1))
dev.off()
# B2
B2_pUC57 <- results(dds, contrast=c("condition","B2","0pUC57"), pAdjustMethod="fdr")
B2_pUC57 <- B2_pUC57[order(B2_pUC57$padj),]
write.table(as.data.frame(B2_pUC57), file="B2_pUC57.tsv", sep = "\t", col.names = T, row.names = T, quote=FALSE)
pdf("B2_pUC57.maplot.pdf", width = 3.5, height = 4)
plotMA(B2_pUC57, main="B2 vs pUC57", alpha=0.05, ylim=c(-1,1))
dev.off()
# C2
C2_pUC57 <- results(dds, contrast=c("condition","C2","0pUC57"), pAdjustMethod="fdr")
C2_pUC57 <- C2_pUC57[order(C2_pUC57$padj),]
write.table(as.data.frame(C2_pUC57), file="C2_pUC57.tsv", sep = "\t", col.names = T, row.names = T, quote=FALSE)
pdf("C2_pUC57.maplot.pdf", width = 3.5, height = 4)
plotMA(C2_pUC57, main="C2 vs pUC57", alpha=0.05, ylim=c(-1,1))
dev.off()
# D3
D3_pUC57 <- results(dds, contrast=c("condition","D3","0pUC57"), pAdjustMethod="fdr")
D3_pUC57 <- D3_pUC57[order(D3_pUC57$padj),]
write.table(as.data.frame(D3_pUC57), file="D3_pUC57.tsv", sep = "\t", col.names = T, row.names = T, quote=FALSE)
pdf("D3_pUC57.maplot.pdf", width = 3.5, height = 4)
plotMA(D3_pUC57, main="D3 vs pUC57", alpha=0.05, ylim=c(-1,1))
dev.off()

# plotting differentially expressed genes
b2sig <- c("SH3BP2","INHBE","AKNA","ADM2","FAM129A","FLNC",
           "SESN2","FABP1","ASNS","CHAC1","UNC5B","LARP6")
d <- filter(counts, symbol %in% b2sig) %>% 
  subset(select=-description) %>% 
  melt() %>%
  separate(variable, into=c("condition","replicate"), sep="_")
pdf("b2sig.pdf", width=3.5, height=12)
ggplot(d, aes(condition, value, color=replicate)) + 
  geom_point(alpha=0.5) + 
  expand_limits(x = 0, y = 0) +
  ylab("Normalized read counts") +
  facet_grid(symbol~., scales="free")
dev.off()

# writing output file
dds_counts <- estimateSizeFactors(dds)
counts <- counts(dds_counts, normalized=TRUE)
counts <- as.data.frame(counts)
counts <- rownames_to_column(counts, var="id")
write.table(counts, file="counts.deseq2.tsv", sep = "\t", col.names = T, row.names = F, quote=FALSE)