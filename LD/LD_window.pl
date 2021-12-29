## Summarize LD in windows
## Kaichi Huang 2020 Dec

use warnings;
use strict;
use POSIX;
use List::Util qw(sum);

my $window_size = 500000; # 500-kbp windows
my %count_hash;
my %sum_hash;
my $current_window = 0;
print "win1\twin2\tn\tmean_r2\n";

while(<STDIN>){
	chomp;
	my @a = split(/\t/,$_);
	my $pos1 = $a[1];
	my $pos2 = $a[3];
	my $r2 = $a[4];
	my $window1 = floor($pos1/$window_size)*$window_size;
	my $window2 = floor($pos2/$window_size)*$window_size;
	$count_hash{$window1}{$window2}++;
	$sum_hash{$window1}{$window2} += $r2;
	if ($window1 != $current_window){
		foreach my $current_window2 (sort {$a<=>$b} keys %{$count_hash{$current_window}}){
			my $mean_r2 = $sum_hash{$current_window}{$current_window2}/$count_hash{$current_window}{$current_window2};
			print "$current_window\t$current_window2\t$count_hash{$current_window}{$current_window2}\t$mean_r2\n";
		}
	}
	$current_window = $window1;
}
