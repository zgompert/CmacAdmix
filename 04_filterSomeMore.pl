#!/usr/bin/perl

use warnings;
use strict;

# filter vcf files

# usage: vcfFilter.pl infile.vcf

# change the marked variables below to adjust settings


############# stringency variables, edit as desired #########################

# maximum depth to avoid repeats (if there are 16 copies of each SNP in a single individual while there are only 2 copies (on average) of every other SNP, this could indicate that you've actually captured a SNP and a paralog or two (in other words, they do NOT actually align together)
#You need to calculate what you want this to be by looking at the distribution of depth coverage in your samples and choosing the mean + 3 SD or some other value that appears to capture what you want
my $maxCoverage =  48000; 

# don't call if SNPs are closer together than this
# can set to a negative number if you want to skip this step altogether
my $mind = 2; 

my $outPrefix = "output_filtered_files/morefilter_maxcov48000_mind2";

#################### input files #######################################

#snps.txt should have one row per snp, column 1 = scaffold (or LG), col. 2 = position in bp
my $SNPtextFile = "output_filtered_files/SNPS_filtered_minCov2x_maxND20_cmac_admix_variants.txt";

# your VCF file that was passed from the shell script: 
my $in = shift(@ARGV);



##### this set is for whole genome shotgun data#

#################### checking SNP distances  ###########################################
my $tooclose;
my $i;

my @line;
my @scaf;
my @pos;
my @keep;


### get set of SNPs and precompute distnaces ##
#snps.txt should have one row per snp, column 1 = scaffold (or LG), col. 2 = position in bp
open(IN, $SNPtextFile) or die "failed to read the locus info file\n";
while(<IN>){
	#MAKE SURE THE LINE BELOW MATCHES THE FORMAT OF YOUR TEXT FILE:
	#This is set up to match a string of digits, hyphens and capital letters (scaffold name) plus a space (data delimiter) plus another string of digits (Scaffold position)
	#Zach's original version is set up to take ONLY a numeric name for scaffold (no characters)
	#The code below is set up for a contig name containing uppercase letters, digits, dots, and vertical bars, then some spaces before the scaffold position, which is comprised solely of digits (\d+)
	m/([A-Z0-9\.\|]+)\s+(\d+)/ or die "failed to match $_\n";
	push(@scaf,$1);
	push(@pos,$2);
	push(@keep,1);
}
close(IN);
my $n = @scaf; 
print "checking neighbors for $n SNPs\n";

#This for-loop is set up in a C-like style, not perl-style
for($i=0; $i<$n; $i++){
	$tooclose = 0;
	unless($i==0){
		#changed == to eq below since my scaffold names are strings, not numerics
		if($scaf[$i] eq $scaf[$i-1] and ($pos[$i] - $pos[$i-1]) < $mind){
			$tooclose = 1;
		}
	}
	unless($i==($n-1)){
		#same here... == to eq
		if($scaf[$i] eq $scaf[$i+1] and ($pos[$i+1] - $pos[$i]) < $mind){
			$tooclose = 1;
		}
	}
	if($tooclose==1){
		$keep[$i] = 0;
	}
}


################### checking read depth ###################################

open (IN, $in) or die "Could not read the infile = $in\n";

# MAKE SURE THE LINE BELOW MATCHES YOUR VCF FILE NAME: 
$in =~ m/output_filtered_files\/([a-zA-Z0-9_]+\.vcf)$/ or die "Failed to match the variant file\n";
open (OUT, "> $outPrefix.$1") or die "Could not write the outfile\n";


my $flag = 0;
my $cnt = 0;

while (<IN>){
	chomp;
	$flag = 1;
	print "\n";
	if (m/^\#/){ ## header row, always write
		$flag = 1;
	}

	#MAKE SURE TO EDIT THE REG. EX. BELOW TO MATCH THE START OF SEQUENCE LINES IN YOUR VCF FILE: 
	#Here I've changed it to match A-Z
	elsif (m/[A-Z]/){ 
		$flag = 1;
		m/DP=(\d+)/ or die "Syntax error, DP not found\n";
		if ($1 > $maxCoverage){
			$flag = 0;
			print "fail DP\n";
		}
		$tooclose = shift(@keep);
		if ($tooclose == 0){
			$flag = 0;
			print "fail too close\n";				
		}
		if ($flag == 1){
			$cnt++; ## this is a good SNV
		}
	}
	else{
		print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
		$flag = 0;
	}
	if ($flag == 1){
		print OUT "$_\n";
	}
}
close (IN);
close (OUT);

print "Finished filtering $in\nRetained $cnt variable loci\n";
