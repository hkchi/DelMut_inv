## Calculate and analyze pi0/pi4
## Kaichi Huang 2020 Nov

library(tidyverse)

win <- 500000

inv_info <- read_tsv("inversion.list", col_names=F)
names(inv_info) <- c("inv","species","chr","start","end")

my.rr <- read_tsv("recombination_rate.tsv", col_names=F) %>% rename(chr=X1, bp=X2, recombination_rate=X3) 

# Read in 0-fold and 4-fold site lists
zero_sites <- read_tsv("Han412.0fold.2.txt", col_names=F)
names(zero_sites) <- c("chr","pos")
four_sites <- read_tsv("Han412.4fold.2.txt", col_names=F)
names(four_sites) <- c("chr","pos")

for (my.species in c("ann", "arg", "pet")) {
  # Sample list
  samples <- read_tsv( paste0(my.species,".samplelist.txt"), col_names=F)
  names(samples) <- c("sample")
  # Genotype table
  gt_table.0 <- read_delim(paste0(my.species,".0fold.table"), delim=" ", col_names=T)
  gt_table.4 <- read_delim(paste0(my.species,".4fold.table"), delim=" ", col_names=T)
  # Calculate pi0/pi4 in sliding windows (samples from all populations)
  gt_table.0.s <- gt_table.0 %>% select(c(CHROM,POS,REF,ALT,samples$sample)) %>%
    unite(all_gt, samples$sample, sep="") %>%
    mutate(c_ref=str_count(all_gt, REF), c_alt=str_count(all_gt, ALT)) %>%
    mutate(n=c_ref+c_alt) %>%
    mutate(heterozygosity=2*(c_ref/n)*(c_alt/n)*n/(n-1)) %>%
    select(CHROM, POS, heterozygosity)
  gt_table.4.s <- gt_table.4 %>% select(c(CHROM,POS,REF,ALT,samples$sample)) %>%
    unite(all_gt, samples$sample, sep="") %>%
    mutate(c_ref=str_count(all_gt, REF), c_alt=str_count(all_gt, ALT)) %>%
    mutate(n=c_ref+c_alt) %>%
    mutate(heterozygosity=2*(c_ref/n)*(c_alt/n)*n/(n-1)) %>%
    select(CHROM, POS, heterozygosity)
  zero_sites.bin <- zero_sites %>% mutate(BIN_START=(pos%/%win)*win) %>%
    group_by(chr, BIN_START) %>%
    summarize(count=n())
  four_sites.bin <- four_sites %>% mutate(BIN_START=(pos%/%win)*win) %>%
    group_by(chr, BIN_START) %>%
    summarize(count=n())
  gt_table.0.s.bin <- gt_table.0.s %>% mutate(BIN_START=(POS%/%win)*win) %>%
    group_by(CHROM, BIN_START) %>%
    summarize(h0=sum(heterozygosity))
  gt_table.4.s.bin <- gt_table.4.s %>% mutate(BIN_START=(POS%/%win)*win) %>%
    group_by(CHROM, BIN_START) %>%
    summarize(h4=sum(heterozygosity))
  pi0.bin <- zero_sites.bin %>% rename(CHROM=chr) %>% inner_join(.,gt_table.0.s.bin) %>% mutate(pi0=h0/count) %>% select(CHROM,BIN_START,pi0)
  pi4.bin <- four_sites.bin %>% rename(CHROM=chr) %>% inner_join(.,gt_table.4.s.bin) %>% mutate(pi4=h4/count) %>% select(CHROM,BIN_START,pi4)
  my.pi0pi4.bin <- inner_join(pi0.bin,pi4.bin) %>% mutate(pi0pi4=pi0/pi4)
  write_tsv(my.pi0pi4.bin, paste0(my.species,".bin_pi0pi4.txt"))
  # pi0/pi4 ~ recombination rate
  my.pi0pi4.bin <- my.pi0pi4.bin %>% inner_join(.,my.rate)
  cor.test(log(my.pi0pi4.bin$pi0pi4), my.pi0pi4.bin$recombination_rate)
  # Compare inversions with genome-wide 
  for (i in inv_info) {
    species <- inv_info$species[i]
    inv_chr <- inv_info$chr[i]
    inv_start <- inv_info$start[i]
    inv_end <- inv_info$end[i]
    if (species == my.species) {
      my.pi0pi4 <- my.pi0pi4.bin %>%
        mutate(category=case_when(recombination_rate*1e6 < 0.01 ~ "null",
                                  (recombination_rate*1e6 >= 0.01) & (recombination_rate*1e6 < 1) ~ "reduced",
                                  TRUE ~ "high")
        ) %>%
        mutate(category=factor(category, levels=c("null","reduced","high"))
        )
      my.pi0pi4.inv <- my.pi0pi4 %>%
        filter(CHROM==inv_chr & BIN_START %in% seq((inv_start%/%win)*win,(inv_end%/%win)*win,win))
      if (length(which(my.pi0pi4.inv$category=="null"))>=10) {
        t.test(log(my.pi0pi4$pi0pi4[which(my.pi0pi4$category=="null")]), log(my.pi0pi4.inv$pi0pi4[which(my.pi0pi4.inv$category=="null")]))
      }
      if (length(which(my.pi0pi4.inv$category=="reduced"))>=10) {
        t.test(log(my.pi0pi4$pi0pi4[which(my.pi0pi4$category=="reduced")]), log(my.pi0pi4.inv$pi0pi4[which(my.pi0pi4.inv$category=="reduced")]))
      }
      if (length(which(my.pi0pi4.inv$category=="high"))>=10) {
       t.test(log(my.pi0pi4$pi0pi4[which(my.pi0pi4$category=="high")]), log(my.pi0pi4.inv$pi0pi4[which(my.pi0pi4.inv$category=="high")]))
      }
    }
  }
  # Calculate pi0pi4 in sliding windows (samples homozygous for the minor arrangement)
  for (i in inv_info) {
    inv <- inv_info$inv[i]
    species <- inv_info$species[i]
    inv_chr <- inv_info$chr[i]
    inv_start <- inv_info$start[i]
    inv_end <- inv_info$end[i]
    if (species == my.species) {
      samples.m <- read_tsv(paste(inv,"minorList.txt",sep="."))
      gt_table.0.s <- gt_table.0 %>% select(c(CHROM,POS,REF,ALT,samples.m$sample)) %>%
        unite(all_gt, samples.m$sample, sep="") %>%
        mutate(c_ref=str_count(all_gt, REF), c_alt=str_count(all_gt, ALT)) %>%
        mutate(n=c_ref+c_alt) %>%
        mutate(heterozygosity=2*(c_ref/n)*(c_alt/n)*n/(n-1)) %>%
        select(CHROM, POS, heterozygosity)
      gt_table.4.s <- gt_table.4 %>% select(c(CHROM,POS,REF,ALT,samples.m$sample)) %>%
        unite(all_gt, samples.m$sample, sep="") %>%
        mutate(c_ref=str_count(all_gt, REF), c_alt=str_count(all_gt, ALT)) %>%
        mutate(n=c_ref+c_alt) %>%
        mutate(heterozygosity=2*(c_ref/n)*(c_alt/n)*n/(n-1)) %>%
        select(CHROM, POS, heterozygosity)
      zero_sites.bin <- zero_sites %>% mutate(BIN_START=(pos%/%win)*win) %>%
        group_by(chr, BIN_START) %>%
        summarize(count=n())
      four_sites.bin <- four_sites %>% mutate(BIN_START=(pos%/%win)*win) %>%
        group_by(chr, BIN_START) %>%
        summarize(count=n())
      gt_table.0.s.bin <- gt_table.0.s %>% mutate(BIN_START=(POS%/%win)*win) %>%
        group_by(CHROM, BIN_START) %>%
        summarize(h0=sum(heterozygosity))
      gt_table.4.s.bin <- gt_table.4.s %>% mutate(BIN_START=(POS%/%win)*win) %>%
        group_by(CHROM, BIN_START) %>%
        summarize(h4=sum(heterozygosity))
      pi0.bin <- zero_sites.bin %>% rename(CHROM=chr) %>% inner_join(.,gt_table.0.s.bin) %>% mutate(pi0=h0/count) %>% select(CHROM,BIN_START,pi0)
      pi4.bin <- four_sites.bin %>% rename(CHROM=chr) %>% inner_join(.,gt_table.4.s.bin) %>% mutate(pi4=h4/count) %>% select(CHROM,BIN_START,pi4)
      my.pi0pi4.bin <- inner_join(pi0.bin,pi4.bin) %>% mutate(pi0pi4=pi0/pi4)
      write_tsv(my.pi0pi4.bin, paste0(inv,".bin_pi0pi4.txt"))
      # Compare inversions with genome-wide
      my.pi0pi4 <- my.pi0pi4.bin %>% filter(pi0!=0 & pi4!=0) %>% inner_join(.,my.rate) %>% 
        mutate(category=case_when(recombination_rate*1e6 < 0.01 ~ "null",
                                  (recombination_rate*1e6 >= 0.01) & (recombination_rate*1e6 < 1) ~ "reduced",
                                  TRUE ~ "high")
        ) %>%
        mutate(category=factor(category, levels=c("null","reduced","high"))
        )
      my.pi0pi4.inv <- my.pi0pi4 %>%
        filter(CHROM==inv_chr & BIN_START %in% seq((inv_start%/%win)*win,(inv_end%/%win)*win,win))
      if (length(which(my.pi0pi4.inv$category=="null"))>=10) {
        t.test(log(my.pi0pi4$pi0pi4[which(my.pi0pi4$category=="null")]), log(my.pi0pi4.inv$pi0pi4[which(my.pi0pi4.inv$category=="null")]))
      }
      if (length(which(my.pi0pi4.inv$category=="reduced"))>=10) {
        t.test(log(my.pi0pi4$pi0pi4[which(my.pi0pi4$category=="reduced")]), log(my.pi0pi4.inv$pi0pi4[which(my.pi0pi4.inv$category=="reduced")]))
      }
      if (length(which(my.pi0pi4.inv$category=="high"))>=10) {
        t.test(log(my.pi0pi4$pi0pi4[which(my.pi0pi4$category=="high")]), log(my.pi0pi4.inv$pi0pi4[which(my.pi0pi4.inv$category=="high")]))
      }
    }
  }
  # Deleterious load in 20 random samples
  my.pi0pi4.spe <- data.frame(name=character(), pi0=numeric(), pi4=numeric(), pi0pi4=numeric())
  total.4 <- nrow(four_sites)
  total.0 <- nrow(zero_sites)
  for (sam in sample(samples$sample,20)) {
    # crazy variable and function names :-)
    pi4 <- gt_table.4 %>%
      unite(gt, sam, sep="") %>%
      mutate(c_ref=str_count(gt, REF), c_alt=str_count(gt, ALT)) %>%
      mutate(heterozygosity=c_ref*c_alt) %>%
      summarize(pi=sum(heterozygosity)/total.4) %>%
      pull(pi)
    pi0 <- gt_table.0 %>%
      unite(gt, sam, sep="") %>%
      mutate(c_ref=str_count(gt, REF), c_alt=str_count(gt, ALT)) %>%
      mutate(heterozygosity=c_ref*c_alt) %>%
      summarize(pi=sum(heterozygosity)/total.0) %>%
      pull(pi)
    tmp.data <- data.frame(name=sam, pi0=pi0, pi4=pi4, pi0pi4=pi0/pi4)
    my.pi0pi4.spe <- rbind(my.pi0pi4.spe, tmp.data)
  }
  write_tsv(my.pi0pi4.spe, paste0(my.species,".sample_pi0pi4.txt"))
}

# Deleterious load in 20 random samples
my.pi0pi4.spe.ann <- read_tsv("ann.sample_pi0pi4.txt") %>% mutate(species="ANN")
my.pi0pi4.spe.arg <- read_tsv("arg.sample_pi0pi4.txt") %>% mutate(species="ARG")
my.pi0pi4.spe.pet <- read_tsv("pet.sample_pi0pi4.txt") %>% mutate(species="PET")
my.pi0pi4.spe.all <- rbind(my.pi0pi4.spe.ann,my.pi0pi4.spe.arg, my.pi0pi4.spe.pet)
my.pi0pi4.spe.all$pi0_s <- my.pi0pi4.spe.all$pi0/sd(my.pi0pi4.spe.all$pi0)
my.pi0pi4.spe.all$pi4_s <- my.pi0pi4.spe.all$pi4/sd(my.pi0pi4.spe.all$pi4)
cor.test((my.pi0pi4.spe.all$pi0_s/my.pi0pi4.spe.all$pi4_s), (my.pi0pi4.spe.all$pi4_s+my.pi0pi4.spe.all$pi0_s)/2)
