## Generate genomic positions of all 4-fold and 0-fold degenerated sites on the genome
## Kaichi Huang 2020 Nov

use warnings;
use strict;

my $cds = $ARGV[0];
my $fasta = $ARGV[1];
my $out_prefix = $ARGV[2];

# First two positions ~ 4-fold
my %hash1 = (
  'TC'=>1,
  'CT'=>1,
  'CC'=>1,
  'CG'=>1,
  'AC'=>1,
  'GT'=>1,
  'GC'=>1,
  'GG'=>1,
  );
# Last two positions ~ 0-fold
my %hash2 = (
  'GA'=>1,
  'GG'=>1,
  'TA'=>1,
  'TG'=>1,
  );

# CDS position and gene info
my %cds;
open CDS, "<$cds" or die $!;
while (<CDS>) {
	chomp;
	my ($chr, $start, $end, $strand, $gene) = (split)[0, 1, 2, 3, 4];
	if (!exists $cds{$gene}) {
		$cds{$gene}{"chr"} = $chr;
		$cds{$gene}{"strand"} = $strand; # strand stay consistent within a gene
		$cds{$gene}{"pos"} = [];
	}
	push @{$cds{$gene}{"pos"}}, ($start, $end);
}
close CDS;

# Loop through CDS sequences
$/=">";
open FASTA, "<$fasta" or die $!;
open OUT_FOUR, ">$out_prefix.4fold.txt";
open OUT_ZERO, ">$out_prefix.0fold.txt";
while (<FASTA>) {
	# Get $name and $seq
	chomp;
	my @a = split(/\n/,$_);
	my $name = $a[0];
	my $seq = join("", splice(@a,1));
	# Exclude sequences that do not begin with "ATG" or end with "TAG/TGA/TAA" or the length of which not multiple of 3
	if ($seq !~ /^ATG/i || $seq !~ /(TAG|TGA|TAA)$/i || length($seq)%3!=0) {next;}
	# Get $gene and $len
	my $gene = (split(/:/,$name))[1];
	my $len = length($seq);
	# Get the chr and pos info
	my $gene_chr = $cds{$gene}{"chr"};
	my @gene_pos = @{$cds{$gene}{"pos"}};
	# Set up
	my ($seg_start,$seg_end);
	if ($cds{$gene}{"strand"} eq "+") {
		$seg_start = shift @gene_pos;
		$seg_end = shift @gene_pos;
	} else {
		# reverse
		$seg_start = pop @gene_pos;
		$seg_end = pop @gene_pos;
	}
	my $delta = 0;
	# Go through the sequence
	for (my $i=0; $i<$len; $i+=3) {
		my ($pos1, $pos2, $pos3);
		if ($cds{$gene}{"strand"} eq "+") {
			$pos1 = $seg_start + ($i - $delta);
			if ($pos1 > $seg_end) {
				$delta = $delta + ($seg_end - $seg_start + 1);
				$seg_start = shift @gene_pos;
				$seg_end = shift @gene_pos;
				$pos1 = $seg_start + ($i - $delta);
			}
			$pos2 = $seg_start + ($i - $delta) + 1;
			if ($pos2 > $seg_end) {
				$delta = $delta + ($seg_end - $seg_start + 1);
				$seg_start = shift @gene_pos;
				$seg_end = shift @gene_pos;
				$pos2 = $seg_start + ($i - $delta) + 1;
			}
			$pos3 = $seg_start + ($i - $delta) + 2;
			if ($pos3 > $seg_end) {
				$delta = $delta + ($seg_end - $seg_start + 1);
				$seg_start = shift @gene_pos;
				$seg_end = shift @gene_pos;
				$pos3 = $seg_start + ($i - $delta) + 2;
			}
		} else {
			# reverse
			$pos1 = $seg_start - ($i - $delta);
			if ($pos1 < $seg_end) {
				$delta = $delta + ($seg_start - $seg_end + 1);
				$seg_start = pop @gene_pos;
				$seg_end = pop @gene_pos;
				$pos1 = $seg_start - ($i - $delta);
			}
			$pos2 = $seg_start - ($i - $delta) - 1;
			if ($pos2 < $seg_end) {
				$delta = $delta + ($seg_start - $seg_end + 1);
				$seg_start = pop @gene_pos;
				$seg_end = pop @gene_pos;
				$pos2 = $seg_start - ($i - $delta) - 1;
			}
			$pos3 = $seg_start - ($i - $delta) - 2;
			if ($pos3 < $seg_end) {
				$delta = $delta + ($seg_start - $seg_end + 1);
				$seg_start = pop @gene_pos;
				$seg_end = pop @gene_pos;
				$pos3 = $seg_start - ($i - $delta) - 2;
			}
		}
		$seq =~ m/(\w)(\w)(\w)/g;
		my $one_two = $1.$2;
		my $two_three = $2.$3;
		# most the 1st position are 0-fold
		if (!exists $hash2{$two_three}) {
			print OUT_ZERO "$gene_chr\t$pos1\n";
		}
		# all the 2nd position are 0-fold
		print OUT_ZERO "$gene_chr\t$pos2\n";
		# only the 3rd position can be 4-fold
		if (exists $hash1{$one_two}) {
			print OUT_FOUR "$gene_chr\t$pos3\n";
		}
	}
}
close FASTA;
close OUT_FOUR;
close OUT_ZERO;
