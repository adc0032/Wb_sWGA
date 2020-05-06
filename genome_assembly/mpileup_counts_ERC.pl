#!/usr/bin/perl

#snpfinderBT2Mpileup.pl by erchan

#Description: Calculates SNPS based on MPILEUP.

#Usage:  perl snpfinderBT2Mpileup.pl <*.mpileup> <PREFIX>

use strict; use warnings;



################SET FILETERS#################

my $filt_mincoverage=0;

#my $filt_maxcoverage=$ARGV[2]; #set based on genome base coverage distribution

#my $majorallelemin=$ARGV[2];  #set after looking at allele frequency MAYBE SET AFTER ANALYZING ALL SAMPLES TO MAKE SURE SNP IS NOT EXCLUDED.  (ie. 80% in one sample is now called reference call)

#my $strandbias=$ARGV[3];  #1/4 of the reads must be on the opposite strand

my $quallimit=0;

################SET FILETERS#################



open (PILE, "<$ARGV[0]") or die "Could not open pileup file.";

#open (INDELS, ">$ARGV[1]INDELS.txt");
#print INDELS join("\t","chromosome","position","INS/DEL","length","\n");

open (RAW, ">$ARGV[1]RAWSNPSPileup.txt");
print RAW join("\t","chromosome","position","reference","major","count_major","count_major_plus","alt","count_alt","count_alt_plus","major_freq","minor_freq","\n");

open (STATS, ">$ARGV[1]RAFSTATSPileup.txt");

open (COVDIS, ">$ARGV[1]COVDISTPileup.txt");





my %reffreqdist;

my %covdist;

my $linect=0;

while (<PILE>) {
    
    chomp;
    
    my $prefilt=0;  #flag filters  0=pass 1=fail
    
    my @ATGC=("A","T","G","C","R","a","t","g","c","r");
    
    my @init=(0,0,0,0,0,0,0,0,0,0);
    
    my %base_ct;
    
    @base_ct{@ATGC} = @init;
    
    my ($chrom, $base, $ref, $coverage, $unfiltcalls, $quals, undef)=split("\t", $_);
    
    ###remove indel and map quality info from base calls
    
    my @genotype=split(//,$unfiltcalls);
    
    my @qual=split(//,$quals);
    
    my $calls="";
    
    for (my $i=0; $i<@genotype; $i++) {
        
        if ($genotype[$i] eq "+") {
            
	#    print INDELS "$chrom\t$base\tINS\t$genotype[$i+1]\n";
            
            my $stringtest=$genotype[$i+1];
            
            if (exists $genotype[$i+2]) {
                
                $stringtest=$genotype[$i+1].$genotype[$i+2];
                
            }
            
            if (exists $genotype[$i+3]) {
                
                $stringtest=$genotype[$i+1].$genotype[$i+2].$genotype[$i+3];
                
            }
            
            my $adjust=0;
            
            if ($stringtest =~ m/(\d+)/) {
                
                $adjust=$1;
                
                if ($adjust>9) {
                    
                    $adjust++;
                    
                }
                
                if ($adjust>99) {
                    
                    $adjust++;
                    
                }
                
            }
            
            $i=($i+1)+$adjust;
            
            next;
            
        }
        
        if ($genotype[$i] eq "-") {
            
	   # print INDELS "$chrom\t$base\tDEL\t$genotype[$i+1]\n";
            
            my $stringtest=$genotype[$i+1];
            
            if (exists $genotype[$i+2]) {
                
                $stringtest=$genotype[$i+1].$genotype[$i+2];
                
            }
            
            if (exists $genotype[$i+3]) {
                
                $stringtest=$genotype[$i+1].$genotype[$i+2].$genotype[$i+3];
                
            }
            
            my $adjust=0;
            
            if ($stringtest =~ m/(\d+)/) {
                
                $adjust=$1;
                
                if ($adjust>9) {
                    
                    $adjust++;
                    
                }
                
                if ($adjust>99) {
                    
                    $adjust++;
                    
                }
                
            }
            
            $i=($i+1)+$adjust;
            
            next;
            
        }
        
        if ($genotype[$i] eq "^") {
            
            $i=$i+1;
            
            next
            
	    }
        
        if ($genotype[$i] eq "\$") {
            
            next;
            
        }
        
        else {
            
            $calls=$calls.$genotype[$i];
            
        }
        
    }
    
    if (length($calls) != length($quals)) {
        
        print "$_\n\n$calls\n$quals\nerror\n";
        
        exit;
        
    }
    
    my @qualcalls=split(//,$calls);
    
    $calls="";
    
    for (my $i=0; $i<@qualcalls; $i++) {
        
        if ((ord($qual[$i])-33) >=$quallimit) {  #ord changes character to ascii value.  the quality scores are given in phred+64.
            
            if ($qualcalls[$i] ne "*") {
                
                $calls=$calls.$qualcalls[$i];
                
            }
            
        }
        
        else {
            
            next;
            
        }
        
    }
    
    if (length($calls) <= $filt_mincoverage) {
        
        next; # $prefilt=1;   #differ from the non all pileup file.
        
    }
    
    $coverage=length($calls);
    
    $base_ct{A}++ while ($calls =~ m/A/g);
    
    $base_ct{T}++ while ($calls =~ m/T/g);
    
    $base_ct{G}++ while ($calls =~ m/G/g);
    
    $base_ct{C}++ while ($calls =~ m/C/g);
    
    $base_ct{a}++ while ($calls =~ m/a/g);
    
    $base_ct{t}++ while ($calls =~ m/t/g);
    
    $base_ct{g}++ while ($calls =~ m/g/g);
    
    $base_ct{c}++ while ($calls =~ m/c/g);
    
    $base_ct{R}++ while ($calls =~ m/\./g);
    
    $base_ct{r}++ while ($calls =~ m/,/g);
    
    
    
    my %base_tot;
    
    $base_tot{A}=$base_ct{A}+$base_ct{a};
    
    $base_tot{T}=$base_ct{T}+$base_ct{t};
    
    $base_tot{G}=$base_ct{G}+$base_ct{g};
    
    $base_tot{C}=$base_ct{C}+$base_ct{c};
    
    $base_tot{R}=$base_ct{R}+$base_ct{r};
    
    my $rank=0;
    
    
    
    print RAW "$chrom\t$base\t$ref";
    
    my $refcalls=$base_tot{R};
    
    my $snpcalls=0;
    
    my $maf=0;
    
    foreach my $check (sort {hashValueAscendingNum (%base_tot)} (keys(%base_tot))) {
        
        $rank++;
        
        if ($rank<=2) {  #controls how many alternate alleles to report
            
            if ($rank==2) {
                
                $maf=$base_tot{$check}/($coverage)*100;
                
            }
            
            if ($base_tot{$check} ==0) {
                
                print RAW "\tNA\t0\t0";
                
                next;
                
            }
            
            my $smallcheck = $check;
            
            $smallcheck =~ tr/A-Z/a-z/;
            
            my $allele=$check;
            
            if ($check eq "R") {
                
                $allele=$ref;
                
            }
            
            if ($check ne "R") {
                
                $snpcalls=$base_tot{$check};
                
            }
            
            my $clonecheck=(($base_ct{$check}+$base_ct{$smallcheck})-abs($base_ct{$check}-$base_ct{$smallcheck}))/(2*($base_ct{$check}+$base_ct{$smallcheck}));
            
            print RAW "\t$allele\t$base_tot{$check}\t$base_ct{$check}";
            
        }
        
        else {
            
            next;
            
        }
        
    }
    
    my $reffreq=$refcalls/($coverage)*100;
    
    print RAW "\t$reffreq\t$maf\n";
    
    if ($coverage >= 20) {
        
        $reffreqdist{int($reffreq)}++;
        
    }
    
    $covdist{$coverage}++;
    
    
    
}



print STATS "$ARGV[1]%\tRef\n";

foreach my $perc (sort {$a<=>$b} keys %reffreqdist) {
    
    print STATS "$perc\t$reffreqdist{$perc}\n";
    
}





print COVDIS "$ARGV[1]\tCOVDIST\n";

foreach my $cov (sort {$a<=>$b} keys %covdist) {
    
    print COVDIS "$cov\t$covdist{$cov}\n";
    
}



sub hashValueAscendingNum {
    
    my (%hash) = @_;
    
    return $hash{$b} <=> $hash{$a};
    
}

close PILE;

#close INDELS;

close STATS;

close RAW;

close COVDIS;

exit;

