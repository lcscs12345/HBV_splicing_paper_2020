library(Gviz)
library(GenomicFeatures)
options(ucscChromosomeNames=FALSE)
HBV <- makeTxDbFromGFF("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/star_hbv_rep1/hbv.orf.gtf")
orf <- GeneRegionTrack(HBV,geneSymbol=TRUE,showId=TRUE)
spHBV <- makeTxDbFromGFF("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/gffcompare/combined.gtf")
sptx <- GeneRegionTrack(spHBV,transcriptAnnotation="symbol",showId=TRUE,fill="red")
A2 <- AlignmentsTrack("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/gffcompare/A2.bam")
B2 <- AlignmentsTrack("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/gffcompare/B2.bam")
C2 <- AlignmentsTrack("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/gffcompare/C2.bam")
D3 <- AlignmentsTrack("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/gffcompare/D3.bam")
pdf("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/gffcompare/hbv_sashimi.pdf",width = 6, height = 14)
plotTracks(c(sptx,orf,A2,B2,C2,D3), chromosome = "HBV",
           type = c("coverage","sashimi"), sashimiNumbers=TRUE, sashimiScore=100,
           sizes=c(10,4,20,20,20,20))
dev.off()


#maxent
sjo <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/SJ.out.offset.txt", header=TRUE, sep="\t")
b2 <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/B2.pgrna.map", header=TRUE, sep="\t")
c2 <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/C2.pgrna.map", header=TRUE, sep="\t")
d3 <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/D3.pgrna.map", header=TRUE, sep="\t")
df <- Reduce(function(...) merge(..., all=TRUE), list(sjo,b2,c2,d3))
df <- subset(df, ss!="NA")

#A2
sj <- subset(df, select=c(A2_pos,ss))
names(sj) <- c("pos","ss")
A2score5 <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/A2.pos.score5.txt", header=T, sep="\t")
A2score3 <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/A2.pos.score3.txt", header=T, sep="\t")
A2ss <- rbind(A2score5,A2score3)
d <- merge(sj,A2ss,by=c("pos","ss"))
library(ggplot2)
pdf("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/A2.score.pdf",width = 10, height = 3)
qplot(data=d, x=pos,xend=pos,y=0,yend=maxent.score) + 
  geom_segment(stat="identity",aes(colour=factor(ss))) +
  geom_label(aes(y=maxent.score-0.005*(max(maxent.score)-min(0)), x=pos, label=pos, colour=factor(ss)), hjust=0.5) + 
  geom_vline(xintercept=1) + geom_vline(xintercept=3246) +
  scale_x_continuous(limits = c(1,3246)) + scale_y_continuous(position = "right") +
  theme_classic()
dev.off()

#B2
sj <- subset(df, select=c(B2_pos,ss))
names(sj) <- c("pos","ss")
B2score5 <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/B2.pos.score5.txt", header=T, sep="\t")
B2score3 <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/B2.pos.score3.txt", header=T, sep="\t")
B2ss <- rbind(B2score5,B2score3)
d <- merge(sj,B2ss,by=c("pos","ss"))
library(ggplot2)
pdf("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/B2.score.pdf",width = 10, height = 3)
qplot(data=d, x=pos,xend=pos,y=0,yend=maxent.score) + 
  geom_segment(stat="identity",aes(colour=factor(ss))) +
  geom_label(aes(y=maxent.score-0.005*(max(maxent.score)-min(0)), x=pos, label=pos, colour=factor(ss)), hjust=0.5) + 
  geom_vline(xintercept=1) + geom_vline(xintercept=3246) +
  scale_x_continuous(limits = c(1,3246)) + scale_y_continuous(position = "right") +
  theme_classic()
dev.off()

#C2
sj <- subset(df, select=c(C2_pos,ss))
names(sj) <- c("pos","ss")
C2score5 <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/C2.pos.score5.txt", header=T, sep="\t")
C2score3 <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/C2.pos.score3.txt", header=T, sep="\t")
C2ss <- rbind(C2score5,C2score3)
d <- merge(sj,C2ss,by=c("pos","ss"))
library(ggplot2)
pdf("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/C2.score.pdf",width = 10, height = 3)
qplot(data=d, x=pos,xend=pos,y=0,yend=maxent.score) + 
  geom_segment(stat="identity",aes(colour=factor(ss))) +
  geom_label(aes(y=maxent.score-0.005*(max(maxent.score)-min(0)), x=pos, label=pos, colour=factor(ss)), hjust=0.5) + 
  geom_vline(xintercept=1) + geom_vline(xintercept=3246) +
  scale_x_continuous(limits = c(1,3246)) + scale_y_continuous(position = "right") +
  theme_classic()
dev.off()

#D3
sj <- subset(df, select=c(D3_pos,ss))
names(sj) <- c("pos","ss")
D3score5 <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/D3.pos.score5.txt", header=T, sep="\t")
D3score3 <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/D3.pos.score3.txt", header=T, sep="\t")
D3ss <- rbind(D3score5,D3score3)
d <- merge(sj,D3ss,by=c("pos","ss"))
library(ggplot2)
pdf("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/maxent/D3.score.pdf",width = 10, height = 3)
qplot(data=d, x=pos,xend=pos,y=0,yend=maxent.score) + 
  geom_segment(stat="identity",aes(colour=factor(ss))) +
  geom_label(aes(y=maxent.score-0.005*(max(maxent.score)-min(0)), x=pos, label=pos, colour=factor(ss)), hjust=0.5) + 
  geom_vline(xintercept=1) + geom_vline(xintercept=3246) +
  scale_x_continuous(limits = c(1,3246)) + scale_y_continuous(position = "right") +
  theme_classic()
dev.off()
