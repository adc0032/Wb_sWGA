#!/usr/bin/env perl
 
use warnings;use strict;
 
use Getopt::Long;
 
my %args = (window => 10, step=>1);
 
die "Usage: $0 [--window XX] [--step XX] [file1 file2...]" unless GetOptions(\%args, 'window=i', 'step=i');
 
local $/ = ">";
readline;
 
while(<>){
    chomp;
    my($id, $sequence) = split /\n/, $_, 2;
    $sequence =~ tr/\n//d;
    $sequence = uc($sequence);
  
    my $end = length($sequence) - $args{window};
    my $start = 0;
    while($start <= $end){
	my $window = substr($sequence, $start, $args{window});
	my $CG = $window =~ tr/CG/CG/;
	print join("\t", $id, $start, $start + $args{window}, $CG/$args{window}*100),"\n";
	$start += $args{step};
    }
}
