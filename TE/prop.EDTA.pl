## Calculate TE density from EDTA .GFF
## Kaichi Huang 2020 Jun

use warnings;
use strict;

my $te = $ARGV[0];
my $chr_len = $ARGV[1];
my $win = $ARGV[2];

my %chr_end;
open CHR, $chr_len;
while(<CHR>){
	chomp;
	my @a = split(/\t/, $_);
	$a[0] =~ s/^>//;
	$chr_end{$a[0]} = $a[1];
}
close CHR;

my %info;
foreach my $chr (keys %chr_end) {
	my $bin = $win;
	$info{$chr}{"bp"} = [];
	while ($bin < $chr_end{$chr}) {
		push @{$info{$chr}{"bp"}}, $bin;
		$bin += $win;
	}
	push @{$info{$chr}{"bp"}}, $chr_end{$chr};
}

my $chr = "";
my ($bp1, $bp2, $bin, $sum);
my $end;
open TE, $te;
while(<TE>){
	chomp;
	my @a = split(/\t/, $_);
	if ($a[0] ne $chr) {
		if ($chr ne "") {
			while ($chr_end{$chr} > $bin) {
				push @{$info{$chr}{"prop"}}, $sum/$win;
				$bin += $win;
				$sum = 0;
			}
			push @{$info{$chr}{"prop"}}, $sum/($chr_end{$chr}-($bin-$win));
		}
		$chr = $a[0];
		$bin = $win;
		$sum = 0;
		$info{$chr}{"prop"} = [];
		$end=0;
	}
	$bp1 = $a[2];
	$bp2 = $a[3];
	if ($bp1 < $end) {
		$bp1 = $end;
	}
	while ($bp2 > $bin) {
		if ($bp1 < $bin) {
			$sum += $bin - $bp1;
			$bp1 = $bin;
		}
		push @{$info{$chr}{"prop"}}, $sum/$win;
		$bin += $win;
		$sum = 0;
	}
	$sum += $bp2 - $bp1;
	$end=$bp2;
}
while ($chr_end{$chr} > $bin) {
	push @{$info{$chr}{"prop"}}, $sum/$win;
	$bin += $win;
	$sum = 0;
}
push @{$info{$chr}{"prop"}}, $sum/($chr_end{$chr}-($bin-$win));
close TE;

$chr = "";
my ($pos, $bin, $sum);
open TE, $te;
while(<TE>){
	chomp;
	my @a = split(/\t/, $_);
	if ($a[0] ne $chr) {
		if ($chr ne "") {
			while ($chr_end{$chr} > $bin) {
				push @{$info{$chr}{"count"}}, $sum;
				$bin += $win;
				$sum = 0;
			}
			push @{$info{$chr}{"count"}}, $sum*$win/($chr_end{$chr}-($bin-$win));
		}
		$chr = $a[0];
		$bin = $win;
		$sum = 0;
		$info{$chr}{"count"} = [];
	}
	$pos = ($a[2] + $a[3]) / 2;
	while ($pos > $bin) {
		push @{$info{$chr}{"count"}}, $sum;
		$bin += $win;
		$sum = 0;
	}
	$sum += 1;
}
while ($chr_end{$chr} > $bin) {
	push @{$info{$chr}{"count"}}, $sum;
	$bin += $win;
	$sum = 0;
}
push @{$info{$chr}{"count"}}, $sum*$win/($chr_end{$chr}-($bin-$win));
close TE;

print "chr\tbp\tTE_density\tTE_count\n";
foreach my $chrom (sort {$a cmp $b} keys %info){
	foreach my $i (0..$#{$info{$chrom}{"bp"}}) {
		print "$chrom\t$info{$chrom}{'bp'}[$i]\t$info{$chrom}{'prop'}[$i]\t$info{$chrom}{'count'}[$i]\n";
	}
}
