library(tidyverse)
library(scales)
library(ggforce)

d <- read.table("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/gffcompare/splice_variants.txt", sep="\t", header=T)
d$spliced_variant <- factor(d$spliced_variant, 
                            levels = c("pSP7","pSP6","pSP5",
                                       "pSP4","SP7","pSP3",
                                       "pSP2","pSP1","SP6",
                                       "SP9","SP1","SP5"))

pdf("~/HBV_splicing_paper_2020/doc/hbv/rnaseq/gffcompare/pie.pdf", width = 8, height = 6)
# https://stackoverflow.com/questions/16184188/ggplot-facet-piechart-placing-text-in-the-middle-of-pie-chart-slices
# adapted from Claus Wilke's solution
dat_pies <- with(d, d[order(genotype, spliced_variant,decreasing=T),])
dat_pies <- left_join(dat_pies,
                      dat_pies %>% 
                        group_by(genotype) %>%
                        summarize(tpm_total = sum(TPM))) %>%
  group_by(genotype) %>%
  mutate(end_angle = 2*pi*cumsum(TPM)/tpm_total,      # ending angle for each pie slice
         start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
         mid_angle = 0.5*(start_angle + end_angle))   # middle of each pie slice, for the text label

rpie = 1 # pie radius

dat_pies <- mutate(dat_pies,
                   hjust = ifelse(mid_angle>pi, 1, 0),
                   vjust = ifelse(mid_angle<pi/2 | mid_angle>3*pi/2, 0, 1))
dat_pies$percent <- percent(dat_pies$TPM/dat_pies$tpm_total,accuracy = 1L)
rlabel = 1.05 * rpie # now we place labels outside of the pies

ggplot(dat_pies) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie,
                   start = start_angle, end = end_angle, fill = spliced_variant)) +
  geom_text(aes(x = rlabel*sin(mid_angle), y = rlabel*cos(mid_angle), label = percent,
                hjust = hjust, vjust = vjust)) +
  coord_fixed() +
  scale_x_continuous(limits = c(-1.5, 1.4), name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1, 1), name = "", breaks = NULL, labels = NULL) +
  facet_wrap(~genotype) + theme_void() + 
  scale_fill_brewer(palette="Set3") +  
  guides(fill = guide_legend(reverse = TRUE))

dev.off()