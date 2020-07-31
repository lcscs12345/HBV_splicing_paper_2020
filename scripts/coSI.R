# compare splicing frequencies
.libPaths(c("/Volumes/scratch/gardnerlab/cslim/R/" , .libPaths()))
library(tidyverse)
library(reshape2)
library(exactRankTests)

a2 <- read.table("~/virus/doc/hbv/rnaseq/maxent/A2.SJ.out.txt", header=TRUE, sep="\t")
b2 <- read.table("~/virus/doc/hbv/rnaseq/maxent/B2.SJ.out.txt", header=TRUE, sep="\t")
c2 <- read.table("~/virus/doc/hbv/rnaseq/maxent/C2.SJ.out.txt", header=TRUE, sep="\t")
d3 <- read.table("~/virus/doc/hbv/rnaseq/maxent/D3.SJ.out.txt", header=TRUE, sep="\t")
a2 <- subset(a2, select=-A2_pos) %>% melt()
b2 <- subset(b2, select=-B2_pos) %>% melt()
c2 <- subset(c2, select=-C2_pos) %>% melt()
d3 <- subset(d3, select=-D3_pos) %>% melt()
hbv <- rbind(a2,b2,c2,d3)
hbv$value <- as.numeric(hbv$value)
hbv$ss <- factor(hbv$ss, levels = c("5ss","3ss"))

ha2 <- read.table("~/virus/doc/hbv/rnaseq/maxent_hg19/A2.SJ.out.txt", header=TRUE, sep="\t")
hb2 <- read.table("~/virus/doc/hbv/rnaseq/maxent_hg19/B2.SJ.out.txt", header=TRUE, sep="\t")
hc2 <- read.table("~/virus/doc/hbv/rnaseq/maxent_hg19/C2.SJ.out.txt", header=TRUE, sep="\t")
hd3 <- read.table("~/virus/doc/hbv/rnaseq/maxent_hg19/D3.SJ.out.txt", header=TRUE, sep="\t")
h57 <- read.table("~/virus/doc/hbv/rnaseq/maxent_hg19/pUC57.SJ.out.txt", header=TRUE, sep="\t")
ha2 <- subset(ha2, select=-A2_pos) %>% melt() %>% mutate_all(str_replace_all, "_reads", "")
hb2 <- subset(hb2, select=-B2_pos) %>% melt() %>% mutate_all(str_replace_all, "_reads", "")
hc2 <- subset(hc2, select=-C2_pos) %>% melt() %>% mutate_all(str_replace_all, "_reads", "")
hd3 <- subset(hd3, select=-D3_pos) %>% melt() %>% mutate_all(str_replace_all, "_reads", "")
h57 <- subset(h57, select=-pUC57_pos) %>% melt() %>% mutate_all(str_replace_all, "_reads", "")
hg <- rbind(ha2,hb2,hc2,hd3,h57)
hg$value <- as.numeric(hg$value)
hg$ss <- factor(hg$ss, levels = c("5ss","3ss"))

pdf("~/virus/doc/hbv/rnaseq/maxent/sj.pdf", width = 3.25, height = 2)
ggplot(hbv, aes(ss,value)) + 
  geom_dotplot(binaxis='y', stackdir='center', fill='white') + 
  stat_summary(fun = "median",  geom = "crossbar", color = "black") + 
  facet_grid(cols = vars(variable)) + scale_y_log10()
dev.off()

pdf("~/virus/doc/hbv/rnaseq/maxent/sj.hg.pdf", width = 4, height = 2)
ggplot(hg, aes(ss,value)) + 
  geom_violin(trim=FALSE) + 
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
  facet_grid(cols = vars(variable)) + scale_y_log10()
dev.off()

perm.test(log(a2[which(a2$ss=='5ss'),]$value), log(a2[which(a2$ss=='3ss'),]$value), alternative='greater')
perm.test(log(b2[which(b2$ss=='5ss'),]$value), log(b2[which(b2$ss=='3ss'),]$value), alternative='greater')
perm.test(log(c2[which(c2$ss=='5ss'),]$value), log(c2[which(c2$ss=='3ss'),]$value), alternative='greater')
perm.test(log(d3[which(d3$ss=='5ss'),]$value), log(d3[which(d3$ss=='3ss'),]$value), alternative='greater')

t.test(log(a2[which(a2$ss=='5ss'),]$value), log(a2[which(a2$ss=='3ss'),]$value), alternative='greater')
t.test(log(b2[which(b2$ss=='5ss'),]$value), log(b2[which(b2$ss=='3ss'),]$value), alternative='greater')
t.test(log(c2[which(c2$ss=='5ss'),]$value), log(c2[which(c2$ss=='3ss'),]$value), alternative='greater')
t.test(log(d3[which(d3$ss=='5ss'),]$value), log(d3[which(d3$ss=='3ss'),]$value), alternative='greater')
t.test(log(hg[which(hg$ss=='5ss'),]$value), log(hg[which(hg$ss=='3ss'),]$value), alternative='greater')


# completed Splicing Index (coSI)
.libPaths(c("/Volumes/scratch/gardnerlab/cslim/R/" , .libPaths()))
library(tidyverse)
library(reshape2)
library(IDPmisc)
setwd("~/virus/doc/hbv/rnaseq/ipsa_out/")
i <- read.table("hbv/intron.zeta.bed", header=F)
names(i) <- c('chr','start','end','cosi3','cosi5','strand')
i <- NaRV.omit(i)
e <-  read.table("hbv/exon.zeta.bed", header=F)
names(e) <- c('chr','start','end','cosi','blank','strand')
e <- NaRV.omit(e)
hi <- read.table("hg19/intron.zeta.bed", header=F)
names(hi) <- c('chr','start','end','cosi3','cosi5','strand')
hi <- NaRV.omit(hi)
he <-  read.table("hg19/exon.zeta.bed", header=F)
names(he) <- c('chr','start','end','cosi','blank','strand')
he <- NaRV.omit(he)

# https://stackoverflow.com/questions/57128090/remove-baseline-color-for-geom-histogram
# Z.Lin's solution
# modified version of StatBin2 inherits from StatBin, except for an
# additional 2nd last line in compute_group() function
StatBin2 <- ggproto(
  "StatBin2", 
  StatBin,
  compute_group = function (data, scales, binwidth = NULL, bins = NULL, 
                            center = NULL, boundary = NULL, 
                            closed = c("right", "left"), pad = FALSE, 
                            breaks = NULL, origin = NULL, right = NULL, 
                            drop = NULL, width = NULL) {
    if (!is.null(breaks)) {
      if (!scales$x$is_discrete()) {
        breaks <- scales$x$transform(breaks)
      }
      bins <- ggplot2:::bin_breaks(breaks, closed)
    }
    else if (!is.null(binwidth)) {
      if (is.function(binwidth)) {
        binwidth <- binwidth(data$x)
      }
      bins <- ggplot2:::bin_breaks_width(scales$x$dimension(), binwidth, 
                                         center = center, boundary = boundary, 
                                         closed = closed)
    }
    else {
      bins <- ggplot2:::bin_breaks_bins(scales$x$dimension(), bins, 
                                        center = center, boundary = boundary, 
                                        closed = closed)
    }
    res <- ggplot2:::bin_vector(data$x, bins, weight = data$weight, pad = pad)
    
    # drop 0-count bins completely before returning the dataframe
    res <- res[res$count > 0, ] 
    
    res
  })

ggplot(e) + 
  geom_histogram(aes(cosi), binwidth = 0.05, color='black', fill=NA, stat = StatBin2) +
  xlim(c(-0.05,1.05)) +
  facet_grid(cols = vars(chr))
ggplot(he) + 
  geom_histogram(aes(cosi), binwidth = 0.05, color='black', fill=NA, stat = StatBin2) +
  xlim(c(-0.05,1.05))

pdf("~/virus/doc/hbv/rnaseq/ipsa_out/cosi.pdf", width = 3.12, height = 2)
ggplot(i) + 
  geom_histogram(aes(cosi3), binwidth = 0.05, boundary = 0, color='#f8766d', alpha=0.2, fill='#f8766d', stat = StatBin2) + 
  geom_histogram(aes(cosi5), binwidth = 0.05, boundary = 0, color='#00bfc4', alpha=0.2, fill='#00bfc4', stat = StatBin2) + 
  xlim(c(-0.05,1.05)) +
  ylim(c(0,13)) +
  facet_grid(cols = vars(chr)) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf("~/virus/doc/hbv/rnaseq/ipsa_out/cosi.hg.pdf", width = 1.41, height = 2)
ggplot(hi) + 
  geom_histogram(aes(cosi3), binwidth = 0.05, boundary = 0, color='#f8766d', alpha=0.2, fill='#f8766d', stat = StatBin2) + 
  geom_histogram(aes(cosi5), binwidth = 0.05, boundary = 0, color='#00bfc4', alpha=0.2, fill='#00bfc4', stat = StatBin2) + 
  xlim(c(-0.05,1.05)) +
  ylim(c(0,100000)) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()