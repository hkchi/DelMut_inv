# pi0pi4

Prepare annotation files

```
# Tables of all 4-fold and 0-fold degenerated sites
# Annotation files downloaded from https://sunflowergenome.org/
perl degeneratedSites.pl CDS_gene_strand.txt CDS.fasta Han412
awk '$1~/^Ha412HOChr[0-9][0-9]$/' Han412.4fold.txt | sort -k1,1 -k2,2n > Han412.4fold.2.txt
awk '$1~/^Ha412HOChr[0-9][0-9]$/' Han412.0fold.txt | sort -k1,1 -k2,2n > Han412.0fold.2.txt

# Tables of sample genotypes
for species in ann arg pet
do
  vcf="$species.vcf.gz"
  bcftools query -l $vcf > $species.samplelist.txt
  header="CHROM POS REF ALT"
  while read sample
  do
  	header="$header $sample"
  done < $species.samplelist.txt
  echo $header > $species.0fold.table
  echo $header > $species.4fold.table
  bcftools query -f '%CHROM %POS %REF %ALT[ %TGT]\n' -R Han412.0fold.2.txt $vcf >> $species.0fold.table
	bcftools query -f '%CHROM %POS %REF %ALT[ %TGT]\n' -R Han412.4fold.2.txt $vcf >> $species.4fold.table
	echo $header > $species.stop.table
	java -Xmx20G -jar snpEff.jar -c snpEff.config -no-utr -no-downstream -no-upstream -no-intergenic HA412 $vcf | grep -E '^#|stop_gained|stop_lost' | bcftools query -f '%CHROM %POS %REF %ALT[ %TGT]\n' >> $species.stop.table
```

Calculate and analyze pi0/pi4
[pi0pi4.R](https://github.com/hkchi/DelMut_inv/tree/master/pi0pi4/pi0pi4.R)

Calculate and analyze Pnonsense/P0
[stop_codon.R](https://github.com/hkchi/DelMut_inv/tree/master/pi0pi4/stop_codon.R)

Compare monomorphic and polymorphic populations
[pop.R](https://github.com/hkchi/DelMut_inv/tree/master/pi0pi4/pop.R)
