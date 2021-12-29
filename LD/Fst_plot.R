## Plots Fst along the inversions
## Kaichi Huang 2020 Dec

library(tidyverse)
library(gridExtra)

inv_info <- read_tsv("inversion_info.tsv", col_names=F)
names(inv_info) <- c("inv","chr","start","end")

plot.list <- list()
for (i in 1:nrow(inv_info)) {
  inv <- inv_info$inv[i]
  inv_chr <- inv_info$chr[i]
  inv_start <- inv_info$start[i]
  inv_end <- inv_info$end[i]
  my.fst.win1 <- read_tsv(paste(inv,"w1000000s200000.windowed.weir.fst",sep="."), col_names=T, na=c("-nan","")) %>%
    mutate(POS=(BIN_START+BIN_END)/2) %>% filter(WEIGHTED_FST>=0)
  my.fst.win2 <- read_tsv(paste(inv,"w100000s20000.windowed.weir.fst",sep="."), col_names=T, na=c("-nan","")) %>%
    mutate(POS=(BIN_START+BIN_END)/2) %>% filter(WEIGHTED_FST>=0)
  p <- ggplot() + theme_classic() +
      geom_line(data=my.fst.win2,mapping=aes(x=POS/1e6, y=WEIGHTED_FST), col="cyan",alpha=.5, na.rm=T) +
      geom_line(data=my.fst.win1,mapping=aes(x=POS/1e6, y=WEIGHTED_FST), col="brown3", alpha=.7, na.rm=T) +
      geom_vline(xintercept=start/1e6, linetype=2) +
      geom_vline(xintercept=end/1e6, linetype=2) +
      xlab(paste(chr,"(Mb)")) +
      {if(i<=5){ylab(expression(italic(F)[ST]))}else{ylab("")}} +
      ggtitle(inv)
  plot.list <- c(plot.list, list(p))
}

pdf("Fst_plot.pdf", width=10, height=10)
print(
  grid.arrange(
    plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],plot.list[[7]],plot.list[[8]],plot.list[[9]],
    layout_matrix = cbind(c(1,2,3,4,5), c(6,7,8,9,NA))
  )
)
dev.off()
