## Examine TE density with recombination rate and inversion
## Kaichi Huang 2020 Dec

library(tidyverse)
library(gridExtra)
library(grid)

win <- 500000
logit <- function(x) {log(x/(1-x))}

# Recombination rate
my.rate <- read_tsv("recombination_rate.tsv", col_names=F) %>% rename(chr=X1, bp=X2, recombination_rate=X3) 

# TE density
my.density <- read_tsv(paste("prop.LTR-RT",win,"tsv",sep="."), col_names=T)
#mean(my.density$TE_density)

# Correlation with recombination rate
my.density <- inner_join(my.density, my.rate)
summary(lm(data=my.density, logit(TE_density)~recombination_rate))

# Compare inversions with genome-wide
inv_info <- read_tsv("inversion_info.tsv", col_names=F)
names(inv_info) <- c("inv","chr","start","end")

my.density <- my.density %>%
  mutate(category=case_when(recombination_rate*1e6 < 0.01 ~ "null",
                            (recombination_rate*1e6 >= 0.01) & (recombination_rate*1e6 < 1) ~ "reduced",
                            TRUE ~ "high")
  ) %>%
  mutate(category=factor(category, levels=c("null","reduced","high"))
  )

i=1 # ann01.01
  inv_chr <- inv_info$chr[i]
  inv_start <- inv_info$start[i]
  inv_end <- inv_info$end[i]
  my.density.inv <- my.density %>%
    filter(chr==inv_chr & bp %in% seq((inv_start%/%win)*win,(inv_end%/%win)*win,win))
  t.test(logit(my.density$TE_density[which(my.density$category=="reduced")]), logit(my.density.inv$TE_density[which(my.density.inv$category=="reduced")]))

i=2 # ann05.01
  inv_chr <- inv_info$chr[i]
  inv_start <- inv_info$start[i]
  inv_end <- inv_info$end[i]
  my.density.inv <- my.density %>%
    filter(chr==inv_chr & bp %in% seq((inv_start%/%win)*win,(inv_end%/%win)*win,win))
  t.test(logit(my.density$TE_density[which(my.density$category=="null")]), logit(my.density.inv$TE_density[which(my.density.inv$category=="null")]))
  t.test(logit(my.density$TE_density[which(my.density$category=="reduced")]), logit(my.density.inv$TE_density[which(my.density.inv$category=="reduced")]))
