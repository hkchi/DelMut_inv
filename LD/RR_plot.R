## Plot recombination rate and centromeres
## Kaichi Huang 2021 Nov

library(tidyverse)

# Recombiantion rates
my.rr <- read_tsv("recombination_rate.tsv", col_names=F) %>% rename(chr=X1, bp=X2, recombination_rate=X3) 

# Centromeres
my.centro <- read_delim("../centromeres/LC075745_HA412.bed", delim =" ", col_names=F) %>%
  rename(chr=X1, start=X2, end=X3) %>%
  filter(grepl("Ha412HOChr[0-1][0-9]$", chr))

# Inversions
inv_info <- read_tsv("inversion_info.tsv", col_names=F)
names(inv_info) <- c("inv","chr","start","end")
inv_info$y <- c(5,5.1,5,5,5,4.9,5,5,5)
inv_info$y2 <- c(5.1,5.15,5.1,5.1,5.1,4.85,5.1,5.1,5.1)

# Plot
pdf("RR_plot.pdf", 8, 5)
ggplot() + theme_bw() +
  geom_tile(data=my.centro,aes(x=(start+end)/2e6, width=2, y=1, height=Inf), fill="lightsteelblue", alpha=0.4) +
  geom_point(data=my.rr, aes(bp/1e6,recombination_rate*1e6), col="plum1", alpha=0.5, size=.75) +
  geom_segment(data=inv_info, aes(x=start/1e6, xend=end/1e6, y=y, yend=y, col=species), size=1.5, alpha=1) +
  facet_wrap( ~ chr, nrow=3, scale = "free_x") +
  scale_color_manual(values=c("#FFB31C","#338732","deepskyblue1"), guide="none") +
  ylab("Recombination rate (cM/Mbp)") + xlab("Mbp") +
  theme(panel.grid=element_blank())
dev.off()
