# overdominance

Examine inversion genotypes and phenotypes to test for overdorminance/underdorminance

[gt_trait.inv.R](https://github.com/hkchi/DelMut_inv/tree/master/overdominance/gt_trait.inv.R)

Calculate U-scores for the inversions

[underdorminance.R](https://github.com/hkchi/DelMut_inv/tree/master/overdominance/underdorminance.R)

Generate genotype tables for 100 random SNPs: 
```
awk '{print $2":"$3-100000"-"$4+100000}' inversion_info.tsv > inversions.intervals
for species in ann arg pet
do
  gatk SelectVariants \
    -V $species.vcf.gz \
    -XL inversions.intervals \
    --exclude-non-variants --exclude-filtered \
    --selectExpressions "vc.isNotFiltered()" \
    --max-nocall-fraction 0.0 \
    -O $species.subset.vcf.gz
done
while read species inv low high
do
  bcftools view -q $low -Q $species.subset.vcf.gz | bcftools query -f '[ %TGT]\n' - | shuf | head -100 | sed 's/|/\//g' > $inv.snps.genotypes.txt
done < inversions.frq # AF file for the inversions
```
