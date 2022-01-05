## Compare deleterious load between monomorphic and polymorphic (MAF>0.2) populations
## Kaichi Huang 2020 Dec

library(tidyverse)

inv_info <- read_tsv("inversion.list", col_names=F)
names(inv_info) <- c("inv","species","chr","start","end")

# Read in 0-fold site list
zero_sites <- read_tsv("Han412.0fold.2.txt",col_names=F)
names(zero_sites) <- c("chr","pos")
# Read in 4-fold site list
four_sites <- read_tsv("Han412.4fold.2.txt",col_names=F)
names(four_sites) <- c("chr","pos")

# Genotype table
gt_table.0 <- read_delim("pet.0fold.table",delim=" ",col_names=T)
gt_table.4 <- read_delim("pet.4fold.table",delim=" ",col_names=T)
gt_table.stop <- read_delim("pet.stop.table",delim=" ",col_names=T)

for (i in c(6,7,8)) {
  inv <- inv_info$inv[i]
  inv_chr <- inv_info$chr[i]
  inv_start <- inv_info$start[i]
  inv_end <- inv_info$end[i]
  # Sample lists
  sample_poly <- read_tsv(paste(inv,"sample_poly.txt",sep="."), col_names=T) %>% pull(Sample) # homozygous samples from polymorphic population
  sample_mono <- read_tsv(paste(inv,"sample_mono.txt",sep="."), col_names=T) %>% pull(Sample) # homozygous samples from monomorphic population
  # Subset to inversion region
  total.0 <- zero_sites %>%
    filter(chr==inv_chr & pos>inv_start & pos<inv_end) %>%
    summarize(count=n()) %>% pull(count)
  total.4 <- four_sites %>%
    filter(chr==inv_chr & pos>inv_start & pos<inv_end) %>%
    summarize(count=n()) %>% pull(count)
  gt_table.0.inv <- gt_table.0 %>% 
    filter(CHROM==inv_chr & POS>inv_start & POS<inv_end)
  gt_table.4.inv <- gt_table.4 %>%
    filter(CHROM==inv_chr & POS>inv_start & POS<inv_end)
  gt_table.stop.inv <- gt_table.stop %>%
    filter(CHROM==inv_chr & POS>inv_start & POS<inv_end)
  # Randomly choose a region of same size on Chr04 as control
  total.0.ctl <- zero_sites %>%
    filter(chr=="Ha412HOChr04" & pos>50000000 & pos<50000000+(inv_end-inv_start)) %>%
    summarize(count=n()) %>% pull(count)
  total.4.ctl <- four_sites %>%
    filter(chr=="Ha412HOChr04" & pos>50000000 & pos<50000000+(inv_end-inv_start)) %>%
    summarize(count=n()) %>% pull(count)
  gt_table.0.ctl <- gt_table.0 %>% 
    filter(CHROM=="Ha412HOChr04" & POS>50000000 & POS<50000000+(inv_end-inv_start))
  gt_table.4.ctl <- gt_table.4 %>%
    filter(CHROM=="Ha412HOChr04" & POS>50000000 & POS<50000000+(inv_end-inv_start))
  gt_table.stop.ctl <- gt_table.stop %>%
    filter(CHROM=="Ha412HOChr04" & POS>50000000 & POS<50000000+(inv_end-inv_start))
  # pi0/pi4
  calculate.pi0pi4 <- function(my.data.frame, samplelist, groupName) {
    for (sam in samplelist) {
      pi4 <- gt_table.4.inv %>%
        unite(gt, sam, sep="") %>%
        mutate(c_ref=str_count(gt, REF), c_alt=str_count(gt, ALT)) %>%
        mutate(heterozygosity=c_ref*c_alt) %>%
        summarize(pi=sum(heterozygosity)/total.4) %>%
        pull(pi)
      pi0 <- gt_table.0.inv %>%
        unite(gt, sam, sep="") %>%
        mutate(c_ref=str_count(gt, REF), c_alt=str_count(gt, ALT)) %>%
        mutate(heterozygosity=c_ref*c_alt) %>%
        summarize(pi=sum(heterozygosity)/total.0) %>%
        pull(pi)
      tmp.data <- data.frame(name=sam, pi0=pi0, pi4=pi4, pi0pi4=pi0/pi4, group=groupName, type="inversion")
      my.data.frame <- rbind(my.data.frame, tmp.data)
      # ctrl
      pi4 <- gt_table.4.ctl %>%
        unite(gt, sam, sep="") %>%
        mutate(c_ref=str_count(gt, REF), c_alt=str_count(gt, ALT)) %>%
        mutate(heterozygosity=c_ref*c_alt) %>%
        summarize(pi=sum(heterozygosity)/total.4) %>%
        pull(pi)
      pi0 <- gt_table.0.ctl %>%
        unite(gt, sam, sep="") %>%
        mutate(c_ref=str_count(gt, REF), c_alt=str_count(gt, ALT)) %>%
        mutate(heterozygosity=c_ref*c_alt) %>%
        summarize(pi=sum(heterozygosity)/total.0) %>%
        pull(pi)
      tmp.data <- data.frame(name=sam, pi0=pi0, pi4=pi4, pi0pi4=pi0/pi4, group=groupName, type="control")
      my.data.frame <- rbind(my.data.frame, tmp.data)
    }
    my.data.frame
  }
  my.pi0pi4.pop <- data.frame(sam=character(), pi0=numeric(), pi4=numeric(), pi0pi4=numeric(), group=factor(), type=factor())
  my.pi0pi4.pop <- calculate.pi0pi4(my.pi0pi4.pop, sample_poly, "poly")
  my.pi0pi4.pop <- calculate.pi0pi4(my.pi0pi4.pop, sample_mono, "mono")
  write_tsv(my.pi0pi4.pop, paste(inv,"pi0pi4_pop.txt",sep="."))
  # Normalize and compare
  control.poly <- mean(my.pi0pi4.pop3$pi0pi4[which(my.pi0pi4.pop3$type=="control"&my.pi0pi4.pop3$group=="poly")])
  control.mono <- mean(my.pi0pi4.pop3$pi0pi4[which(my.pi0pi4.pop3$type=="control"&my.pi0pi4.pop3$group=="mono")])
  pi0pi4.ctrl.poly <- data.frame(pi0pi4.c=my.pi0pi4.pop3$pi0pi4[which(my.pi0pi4.pop3$type=="inversion"&my.pi0pi4.pop3$group=="poly")] / control.poly, group="poly")
  pi0pi4.ctrl.mono <- data.frame(pi0pi4.c=my.pi0pi4.pop3$pi0pi4[which(my.pi0pi4.pop3$type=="inversion"&my.pi0pi4.pop3$group=="mono")] / control.mono, group="mono")
  my.pi0pi4.pop3.new <- rbind( pi0pi4.ctrl.poly, pi0pi4.ctrl.mono)
  t.test(my.pi0pi4.pop3.new$pi0pi4.c[which(my.pi0pi4.pop3.new$group=="poly")], my.pi0pi4.pop3.new$pi0pi4.c[which(my.pi0pi4.pop3.new$group=="mono")])
  # Pnonsense/P0
  calculate.pp <- function(my.data.frame, samplelist, groupName) {
    for (sam in samplelist) {
      total0 <- gt_table.0.inv %>%
        unite(gt, sam, sep="") %>%
        mutate(c_ref=str_count(gt, REF), c_alt=str_count(gt, ALT)) %>%
        mutate(poly=c_ref*c_alt) %>%
        summarize(count=sum(poly)) %>%
        pull(count)
      totals <- gt_table.stop.inv %>%
        unite(gt, sam, sep="") %>%
        mutate(c_ref=str_count(gt, REF), c_alt=str_count(gt, ALT)) %>%
        mutate(poly=c_ref*c_alt) %>%
        summarize(count=sum(poly)) %>%
        pull(count)
      tmp.data <- data.frame(name=sam, t0=total0, ts=totals, pp=totals/total0, group=groupName, type="inversion")
      my.data.frame <- rbind(my.data.frame, tmp.data)
      # ctrl
      total0 <- gt_table.0.ctl %>%
        unite(gt, sam, sep="") %>%
        mutate(c_ref=str_count(gt, REF), c_alt=str_count(gt, ALT)) %>%
        mutate(poly=c_ref*c_alt) %>%
        summarize(count=sum(poly)) %>%
        pull(count)
      totals <- gt_table.stop.ctl %>%
        unite(gt, sam, sep="") %>%
        mutate(c_ref=str_count(gt, REF), c_alt=str_count(gt, ALT)) %>%
        mutate(poly=c_ref*c_alt) %>%
        summarize(count=sum(poly)) %>%
        pull(count)
      tmp.data <- data.frame(name=sam, t0=total0, ts=totals, pp=totals/total0, group=groupName, type="control")
      my.data.frame <- rbind(my.data.frame, tmp.data)
    }
    my.data.frame
  }
  my.pp.pop <- data.frame(name=character(), t0=numeric(), ts=numeric(), pp=numeric(), group=factor(), type=factor())
  my.pp.pop <- calculate.pi0pi4(my.pp.pop, sample_poly, "poly")
  my.pp.pop <- calculate.pi0pi4(my.pp.pop, sample_mono, "mono")
  write_tsv(my.pi0pi4.pop, paste(inv,"pp_pop.txt",sep="."))
  # Normalize and compare
  control.poly <- mean(my.pp.pop3$pp[which(my.pp.pop3$type=="control"&my.pp.pop3$group=="poly")])
  control.mono <- mean(my.pp.pop3$pp[which(my.pp.pop3$type=="control"&my.pp.pop3$group=="mono")])
  pp.ctrl.poly <- data.frame(pp.c=my.pp.pop3$pp[which(my.pp.pop3$type=="inversion"&my.pp.pop3$group=="poly")] / control.poly, group="poly")
  pp.ctrl.mono <- data.frame(pp.c=my.pp.pop3$pp[which(my.pp.pop3$type=="inversion"&my.pp.pop3$group=="mono")] / control.mono, group="mono")
  my.pp.pop3.new <- rbind(pp.ctrl.poly, pp.ctrl.mono)
  t.test(my.pp.pop3.new$pp.c[which(my.pp.pop3.new$group=="poly")], my.pp.pop3.new$pp.c[which(my.pp.pop3.new$group=="mono")])
}
