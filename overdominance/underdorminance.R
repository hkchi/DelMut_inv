## Calculate U-scores for the inversions
## Kaichi Huang 2021 Oct

library(tidyverse)
library(HWxtest)

inv_info <- read_tsv("inversion_info.tsv", col_names=F)
names(inv_info) <- c("inv","chr","start","end")

sample_info <- read_tsv("sample_info_apr_2018.tsv", col_names=T) %>% select(name,population) %>% rename(Sample=name, Population=population)

# Population U-scores
my.Us <- data.frame()
for (i in 1:nrow(inv_info)) {
  inv <- inv_info$inv[i]
  species <- inv_info$species[i]
  samples <- read_tsv( paste0(species,".samplelist.txt"), col_names=F)
  names(samples) <- c("sample")
  genotype_info <- read_tsv(paste(species,inv,"genotypes.txt",sep="."), col_names=T) %>% select(sample,triangle_genotype) %>%
    inner_join(.,sample_list) %>% inner_join(.,sample_info)
  genotype_info2 <- genotype_info %>% mutate(locus1 = case_when(triangle_genotype == 0 ~ '0/0',
                                                                triangle_genotype == 1 ~ "0/1",
                                                                triangle_genotype == 2 ~ "1/1",
                                                                TRUE ~ "NA")) %>%
    select(-Sample, -triangle_genotype) %>% rename(pop=Population) # Genotype table for hwx.test
  hw_result <- hwx.test(genotype_info2, statName="U")
  hw_utest <- hwdf(hw_result)
  hw_utest$tmp <- rownames(hw_utest)
  hw_utest <- hw_utest %>% separate(tmp, c("Population","locus"), "\\.")
  
  tmp.Us <- hw_utest %>% select(p.value="P-val(U)", u.score="obs-U", n="N") %>% mutate(locus=inv)
  my.Us <- rbind(my.Us, tmp.Us)
}
my.Us <- my.Us %>% mutate(p.type=case_when(p.value<0.01~"significant",TRUE~"nonsignificant"))
write_tsv(my.Us, "my.Us.txt")

# Species-wide U-scores
my.Us.all.inv <- data.frame()
my.Us.all.snps <- data.frame()
for (i in 1:nrow(inv_info)) {
  inv <- inv_info$inv[i]
  species <- inv_info$species[i]
  samples <- read_tsv( paste0(species,".samplelist.txt"), col_names=F)
  names(samples) <- c("sample")
  genotype_info <- read_tsv(paste(species,inv,"genotypes.txt",sep="."), col_names=T) %>% select(sample,triangle_genotype) %>%
    inner_join(.,sample_list) %>% inner_join(.,sample_info)
  genotype_info2 <- genotype_info %>% mutate(locus1 = case_when(triangle_genotype == 0 ~ '0/0',
                                                                triangle_genotype == 1 ~ "0/1",
                                                                triangle_genotype == 2 ~ "1/1",
                                                                TRUE ~ "NA")) %>%
    select(-Sample, -triangle_genotype) %>% rename(pop=Population)
  # Permutate the samples
  tmp.Us.all.inv <- c()
  for (j in 1:100) {
    genotype_info2.tmp <- genotype_info2 %>% group_by(pop) %>% sample_n(1) %>% ungroup()
    genotype_info2.tmp$pop <- "All"
    hw_result <- hwx.test(genotype_info2.tmp, statName="U")
    hw_utest <- hwdf(hw_result)
    tmp.Us.all.inv.tmp <- hw_utest %>% select(p.value="P-val(U)", u.score="obs-U", n="N") %>% mutate(locus=inv)
    tmp.Us.all.inv <- rbind(tmp.Us.all.inv, tmp.Us.all.inv.tmp)
  }
  my.Us.all.inv <- rbind(my.Us.all.inv,tmp.Us.all.inv)
  # 100 random SNPs
  genotype_info <- read.table(paste(inv,"snps.genotypes.txt",sep="."), stringsAsFactors=F) %>% as.matrix() %>% t() %>% as.data.frame(.,stringsAsFactors=F)
  names(genotype_info) <- paste0("locus",1:ncol(genotype_info))
  row.names(genotype_info) <- sample_list$Sample
  genotype_info2 <- cbind(sample_list, genotype_info) %>% inner_join(sample_info) %>% rename(pop=Population)
  genotype_info2.tmp <- genotype_info2 %>% group_by(pop) %>% sample_n(1) %>% ungroup()
  genotype_info2.tmp$pop <- "All"
  hw_result <- hwx.test(genotype_info2.tmp, statName="U")
  hw_utest <- hwdf(hw_result)
  tmp.Us.all.snps <- hw_utest %>% select(p.value="P-val(U)", u.score="obs-U", n="N") %>% mutate(locus=inv)
  my.Us.all.snps <- rbind(my.Us.all.snps, tmp.Us.all.snps)
}
write_tsv(my.Us.all.inv, "my.Us.all.inv.txt")
write_tsv(my.Us.all.snps, "my.Us.all.snps.txt")
