# TE

Summarize density of transposable elements (LTR-retrotransposons) across the chromosomes

```
win=500000
gff="EDTA.TEanno.gff" # Output from EDTA
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }' $gff | awk '{print $1"\t"length($2)}' > chr_len.txt
grep -v "^#" $gff | awk '$9~/^ID=/{print}' | awk '$3~/LTR/{print $1"\t"$3"\t"$4"\t"$5}' | sort -k 1,1 -k 3n > LTR-RT.tsv

perl prop.EDTA.pl LTR-RT.tsv chr_len.txt $win > prop.LTR-RT.$win.tsv
```

Examine TE density with recombination rate and inversion

```
Rscript TE_RR.R
```

