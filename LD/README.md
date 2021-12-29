# LD

Plot the distribution of recombination rate: 

```
Rscript RR_plot.R
```


Convert vcf to plink format: 

```
for species in ann arg pet
do
	plink --vcf $species.vcf.gz --allow-extra-chr --double-id --out $species --threads 5
done
```

Calculate LD for pairwise windows across the chromosome: 

```
while read inv species chr start end
do
	plink --r2 --bfile $species --allow-extra-chr --keep $inv.homo.keeplist.txt --chr $chr --maf 0.1 --bp-space 5000 --ld-window 999999999 --ld-window-kb 999999 --ld-window-r2 0 --out LD_${inv}_homo --threads 48
	plink --r2 --bfile $species --allow-extra-chr --keep $inv.poly.keeplist.txt --chr $chr --maf 0.1 --bp-space 5000 --ld-window 999999999 --ld-window-kb 999999 --ld-window-r2 0 --out LD_${inv}_poly --threads 48
	awk -F "[ ]+" 'NR>1{print $1"\t"$2"\t"$4"\t"$5"\t"$7}' LD_${inv}_homo.ld | perl LD_window.pl > LD_${inv}_homo.window.ld
	awk -F "[ ]+" 'NR>1{print $1"\t"$2"\t"$4"\t"$5"\t"$7}' LD_${inv}_poly.ld | perl LD_window.pl > LD_${inv}_poly.window.ld
done < inversion.list
```

Plot LD: 

```
Rscript LD_plot.R
```


Calculate Fst between arrangement for each inversion: 

```
while read inv species chr start end
do
	from=$((start - 20000000))
	to=$((end + 20000000))
	vcftools --gzvcf $species.vcf --chr $chr --from-bp $from --to-bp $to --weir-fst-pop $inv.group_0.txt --weir-fst-pop $inv.group_2.txt --fst-window-size 1000000 --fst-window-step 200000 --out $inv.w1000000s200000
	vcftools --gzvcf $species.vcf --chr $chr --from-bp $from --to-bp $to --weir-fst-pop $inv.group_0.txt --weir-fst-pop $inv.group_2.txt --fst-window-size 100000 --fst-window-step 20000 --out $inv.w100000s20000
done < inversion.list
```

Plot Fst: 

```
Rscript Fst_plot.R
```