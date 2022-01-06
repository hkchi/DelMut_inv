## Examine inversion genotypes and phenotypes
## Kaichi Huang 2021 Nov

library(tidyverse)

inv_info <- read_tsv("inversion_info.tsv", col_names=F)
names(inv_info) <- c("inv","chr","start","end")

my.data <- data.frame()

for (i in 1:nrow(inv_info)) {
  inv <- inv_info$inv[i]
  species <- inv_info$species[i]
  samples <- read_tsv( paste0(species,".samplelist.txt"), col_names=F)
  names(samples) <- c("sample")
  genotype_info <- read_tsv(paste(species,inv,"genotypes.txt",sep="."), col_names=T) %>% select(sample,triangle_genotype) %>%
    inner_join(.,sample_list)
  my.traits <- system(paste0('ls gwas/',species),intern=T)
  my.traits <- gsub(".znorm.txt", "", my.traits[grep("znorm",my.traits)])
  for (trait in my.traits) {
    tmp.data <- read_tsv(paste0("gwas/",species,"/",trait,".znorm.txt"), col_names=F) %>% select(Sample=X1, value=X3) %>%
      inner_join(.,genotype_info)
    mean.00 <- mean(tmp.data$value[which(tmp.data$triangle_genotype==0)],na.rm=T)
    mean.01 <- mean(tmp.data$value[which(tmp.data$triangle_genotype==1)],na.rm=T)
    mean.11 <- mean(tmp.data$value[which(tmp.data$triangle_genotype==2)],na.rm=T)
    if (mean.01 > max(mean.00, mean.11)) {
      # possible overdominance
      if (mean.00 > mean.11) {
        pvalue <- t.test(tmp.data$value[which(tmp.data$triangle_genotype==1)], tmp.data$value[which(tmp.data$triangle_genotype==0)])$p.value
      } else {
        pvalue <- t.test(tmp.data$value[which(tmp.data$triangle_genotype==1)], tmp.data$value[which(tmp.data$triangle_genotype==2)])$p.value
      }
      my.data <- rbind(my.data, data.frame(species=species, inversion=inv, trait=trait, type="overdominance", p.value=pvalue))
    } else if (mean.01 < min(mean.00, mean.11)) {
      # possible underdominance
      if (mean.00 < mean.11) {
        pvalue <- t.test(tmp.data$value[which(tmp.data$triangle_genotype==1)], tmp.data$value[which(tmp.data$triangle_genotype==0)])$p.value
      } else {
        pvalue <- t.test(tmp.data$value[which(tmp.data$triangle_genotype==1)], tmp.data$value[which(tmp.data$triangle_genotype==2)])$p.value
      }
      my.data <- rbind(my.data, data.frame(species=species, inversion=inv, trait=trait, type="underdominance", p.value=pvalue))
    }
  }
}

write_tsv(my.data, "gt_trait.inv.txt")
