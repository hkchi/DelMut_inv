## Plot window-based LD
## Kaichi Huang 2020 Dec

library(tidyverse)

inv_info <- read_tsv("inversion_info.tsv", col_names=F)
names(inv_info) <- c("inv","chr","start","end")

my.centro <- read_delim("../centromeres/LC075745_HA412.bed", delim =" ", col_names=F) %>%
    rename(chr=X1, start=X2, end=X3) %>%
    filter(grepl("Ha412HOChr[0-1][0-9]$", chr))

pdf("LD_plot.pdf", width=7, height=7)
for (i in 1:nrow(inv_info)) {
  inv <- inv_info$inv[i]
  inv_chr <- inv_info$chr[i]
  inv_start <- inv_info$start[i]
  inv_end <- inv_info$end[i]
  my.centro.chr <- my.centro %>% filter(chr==inv_chr)
  my.ld2.poly <- read_tsv(paste0("LD_",inv,"_poly",".window.ld"))
  my.ld2.homo <- read_tsv(paste0("LD_",inv,"_homo",".window.ld"))
  print(
    ggplot(my.ld2.poly, aes(x=win1/1e6, y=win2/1e6)) + theme_classic() +
      geom_tile(aes(fill=mean_r2)) +
      geom_tile(data=my.ld2.homo, aes(x=win2/1e6,y=win1/1e6,fill=mean_r2)) +
      scale_fill_gradientn(colours=c("grey95","blue","red"), values=c(0,0.5,1), name="LD") +
      geom_segment(mapping=aes(x=inv_start/1e6,xend=inv_end/1e6,y=-4,yend=-4), col="purple",size = 2.5) +
      geom_segment(mapping=aes(x=-4,xend=-4,y=inv_start/1e6,yend=inv_end/1e6), col="purple",size = 2.5) +
      geom_point(data=my.centro.chr, mapping=aes(x=(start+end)/2e6, y=-2), col="red", shape=18, size=3, alpha=0.6) +
      geom_point(data=my.centro.chr, mapping=aes(x=-2, y=(start+end)/2e6), col="red", shape=18, size=3, alpha=0.6) +
      scale_x_continuous(expand=c(0.02,0)) +
      scale_y_continuous(expand=c(0.02,0)) +
      coord_fixed(ratio = 1) +
      xlab("Mbp") + ylab("Mbp") +
      ggtitle(inv)
  )
}
dev.off()
