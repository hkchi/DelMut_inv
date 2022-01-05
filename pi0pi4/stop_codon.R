## Calculate and analyze Pnonsense/P0
## Kaichi Huang 2020 Dec

library(tidyverse)

win <- 500000

inv_info <- read_tsv("inversion.list", col_names=F)
names(inv_info) <- c("inv","species","chr","start","end")

my.rr <- read_tsv("recombination_rate.tsv", col_names=F) %>% rename(chr=X1, bp=X2, recombination_rate=X3) 

for (my.species in c("ann", "arg", "pet")) {
  # Sample list
  samples <- read_tsv( paste0(my.species,".samplelist.txt"), col_names=F)
  names(samples) <- c("sample")
  # Genotype table
  gt_table.0 <- read_delim(paste0(my.species,".0fold.table"), delim=" ", col_names=T)
  gt_table.stop <- read_delim(paste0(my.species,".stop.table"),delim=" ",col_names=T)
  # Calculated pi0/pi4 in sliding windows (samples from all populations)
  gt_table.0.s <- gt_table.0 %>% select(c(CHROM,POS,REF,ALT,samples$sample)) %>%
    unite(all_gt, samples$sample, sep="") %>%
    mutate(c_ref=str_count(all_gt, REF), c_alt=str_count(all_gt, ALT)) %>%
    mutate(poly=c_ref*c_alt) %>%
    select(CHROM, POS, poly) %>%
    filter(poly!=0)
  gt_table.stop.s <- gt_table.stop %>% select(c(CHROM,POS,REF,ALT,samples$sample)) %>%
    unite(all_gt, samples$sample, sep="") %>%
    mutate(c_ref=str_count(all_gt, REF), c_alt=str_count(all_gt, ALT)) %>%
    mutate(poly=c_ref*c_alt) %>%
    select(CHROM, POS, poly) %>%
    filter(poly!=0)
  gt_table.0.s.bin <- gt_table.0.s %>% mutate(BIN_START=(POS%/%win)*win) %>%
    group_by(CHROM, BIN_START) %>%
    summarize(count0=n())
  gt_table.stop.s.bin <- gt_table.stop.s %>% mutate(BIN_START=(POS%/%win)*win) %>%
    group_by(CHROM, BIN_START) %>%
    summarize(counts=n())
  my.pp.bin <- gt_table.stop.s.bin %>% inner_join(.,gt_table.0.s.bin) %>% 
    mutate(pp=counts/count0)
  write_tsv(my.pp.bin, paste0(my.species,".bin_pp.txt"))
  # Pnonsense/P0 ~ recombination rate
  my.pp.bin <- my.pp.bin %>% inner_join(.,my.rate)
  cor.test(log(my.pp.bin$pp), my.pp.bin$recombination_rate)
  # Compare inversions with genome-wide 
  for (i in inv_info) {
    species <- inv_info$species[i]
    inv_chr <- inv_info$chr[i]
    inv_start <- inv_info$start[i]
    inv_end <- inv_info$end[i]
    if (species == my.species) {
      my.pp <- my.pp.bin %>%
        mutate(category=case_when(recombination_rate*1e6 < 0.01 ~ "null",
                                  (recombination_rate*1e6 >= 0.01) & (recombination_rate*1e6 < 1) ~ "reduced",
                                  TRUE ~ "high")
        ) %>%
        mutate(category=factor(category, levels=c("null","reduced","high"))
        )
      my.pp.inv <- my.pp %>%
        filter(CHROM==inv_chr & BIN_START %in% seq((inv_start%/%win)*win,(inv_end%/%win)*win,win))
      if (length(which(my.pp.inv$category=="null"))>=10) {
        t.test(log(my.pp$pp[which(my.pp$category=="null")]), log(my.pp.inv$pp[which(my.pp.inv$category=="null")]))
      }
      if (length(which(my.pp.inv$category=="reduced"))>=10) {
        t.test(log(my.pp$pp[which(my.pp$category=="reduced")]), log(my.pp.inv$pp[which(my.pp.inv$category=="reduced")]))
      }
      if (length(which(my.pp.inv$category=="high"))>=10) {
        t.test(log(my.pp$pp[which(my.pp$category=="high")]), log(my.pp.inv$pp[which(my.pp.inv$category=="high")]))
      }
    }
  }
  # Calculate pp in sliding windows (samples homozygous for the minor arrangement)
  for (i in inv_info) {
    inv <- inv_info$inv[i]
    species <- inv_info$species[i]
    inv_chr <- inv_info$chr[i]
    inv_start <- inv_info$start[i]
    inv_end <- inv_info$end[i]
    if (species == my.species) {
      samples.m <- read_tsv(paste(inv,"minorList.txt",sep="."))
      gt_table.0.s <- gt_table.0 %>% select(c(CHROM,POS,REF,ALT,samples$sample)) %>%
        unite(all_gt, samples$sample, sep="") %>%
        mutate(c_ref=str_count(all_gt, REF), c_alt=str_count(all_gt, ALT)) %>%
        mutate(poly=c_ref*c_alt) %>%
        select(CHROM, POS, poly) %>%
        filter(poly!=0)
      gt_table.stop.s <- gt_table.stop %>% select(c(CHROM,POS,REF,ALT,samples$sample)) %>%
        unite(all_gt, samples$sample, sep="") %>%
        mutate(c_ref=str_count(all_gt, REF), c_alt=str_count(all_gt, ALT)) %>%
        mutate(poly=c_ref*c_alt) %>%
        select(CHROM, POS, poly) %>%
        filter(poly!=0)
      gt_table.0.s.bin <- gt_table.0.s %>% mutate(BIN_START=(POS%/%win)*win) %>%
        group_by(CHROM, BIN_START) %>%
        summarize(count0=n())
      gt_table.stop.s.bin <- gt_table.stop.s %>% mutate(BIN_START=(POS%/%win)*win) %>%
        group_by(CHROM, BIN_START) %>%
        summarize(counts=n())
      my.pp.bin <- gt_table.stop.s.bin %>% inner_join(.,gt_table.0.s.bin) %>% 
        mutate(pp=counts/count0)
      write_tsv(my.pp.bin, paste0(inv,"bin_pp.txt"))
      # Compare inversions with genome-wide
      my.pp <- my.pp.bin %>% filter(pi0!=0 & pi4!=0) %>% inner_join(.,my.rate) %>% 
        mutate(category=case_when(recombination_rate*1e6 < 0.01 ~ "null",
                                  (recombination_rate*1e6 >= 0.01) & (recombination_rate*1e6 < 1) ~ "reduced",
                                  TRUE ~ "high")
        ) %>%
        mutate(category=factor(category, levels=c("null","reduced","high"))
        )
      my.pp.inv <- my.pp %>%
        filter(CHROM==inv_chr & BIN_START %in% seq((inv_start%/%win)*win,(inv_end%/%win)*win,win))
      if (length(which(my.pp.inv$category=="null"))>=10) {
        t.test(log(my.pp$pp[which(my.pp$category=="null")]), log(my.pp.inv$pp[which(my.pp.inv$category=="null")]))
      }
      if (length(which(my.pp.inv$category=="reduced"))>=10) {
        t.test(log(my.pp$pp[which(my.pp$category=="reduced")]), log(my.pp.inv$pp[which(my.pp.inv$category=="reduced")]))
      }
      if (length(which(my.pp.inv$category=="high"))>=10) {
        t.test(log(my.pp$pp[which(my.pp$category=="high")]), log(my.pp.inv$pp[which(my.pp.inv$category=="high")]))
      }
    }
  }
}
