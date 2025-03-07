#!/usr/bin/perl

use warnings;
use strict;

# this program filters a vcf file based on overall sequence coverage, number of non-reference reads, number of alleles, and reverse orientation reads

# usage vcfFilter.pl infile.vcf
#
# change the marked variables below to adjust settings
#

#### stringency variables, edits as desired
my $minCoverage = 3072; # minimum number of sequences; DP, 2x = 1536*2 = 3072 (2x number of samples)
my $minAltRds = 10; # minimum number of sequences with the alternative allele; AC
my $notFixed = 1.0; # removes loci fixed for alt; AF

## tests are p-values here
my $bqrs = 0.001; # maximum absolute value of the base quality rank sum test; BQB
my $mqrs = 0.0001; # maximum absolute value of the mapping quality rank sum test; MQB
my $rprs = 0.001; # minimum absolute value of the read position rank sum test; RPB

my $mq = 30; #mapping quality; typically 20-30
my $miss = 614; # maximum number of individuals with no data, 80+ = 3072*0.2 = 614.4

##### this set is for whole genome shotgun data
my $d;

my @line;

my $in = shift(@ARGV); #input .vcf file must be passed to perl script here from the shell script
open (IN, $in) or die "Could not read the infile = $in\n";
#make sure the regular expression matching below will capture the name of your $in file:
$in =~ m/input_vcf_file\/([a-z\_]+\.vcf)$/ or die "Failed to match the variant file\n";
open (OUT, "> output_filtered_files\/filtered_minCov2x_maxND20_$1") or die "Could not write the outfile\n";

my $flag = 0;
my $cnt = 0;

while (<IN>){
	chomp;
	$flag = 1;
	if (m/^\#/){ ## header row, always write (matches any row starting with a "#")
		$flag = 1;
	}
	elsif (m/^[A-Z]/){ ## this is a sequence line, you might need to edit this reg. expr.
			   ##the '^' up carot mark means "start of line" for matching
			   ## in my case, all the contig IDs for beetles start with a capital letter
		$flag = 1;
		# $d matches the following: digit/digit:0,0,0:digits
		$d = () = (m/\d\/\d:0,0,0:\d+/g); ## for bcftools call
		if ($d >= $miss){
			$flag = 0;
			##print "fail missing : ";
		}
		if (m/[ACTGN]\,[ACTGN]/){ ## two alternative alleles identified
			$flag = 0;
			#print "fail allele : ";
		}
		@line = split(/\s+/,$_);
		if(length($line[3]) > 1 or length($line[4]) > 1){
			$flag = 0;
			#print "fail INDEL : ";
		}
		m/DP=(\d+)/ or die "Syntax error, DP not found\n";
		if ($1 < $minCoverage){
			$flag = 0;
			#print "fail DP : ";
		}
		m/AC1*=(\d+)/ or die "Syntax error, DP not found\n";
		if ($1 < $minAltRds){
			$flag = 0;
			#print "fail DP : ";
		}
## bcftools call version
	
		m/DP4=\d+,\d+,(\d+),(\d+)/ or die "Syntax error DP4 not found\n";
		if(($1 + $2) < $minAltRds){
			$flag = 0;
		}
		m/AF1*=([0-9\.e\-]+)/ or die "Syntax error, AF not found\n";
		if ($1 == $notFixed){
			$flag = 0;
		#	print "fail AF : ";
		}

		if(m/MQ=([0-9\.]+)/){
			if ($1 < $mq){
				$flag = 0;
#				print "fail MQ : ";
			}
		}
		else{
			$flag = 0;
			print "faile no MQ : ";
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
