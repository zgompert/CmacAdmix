#!/usr/bin/perl

use warnings;
use strict;

# this script takes an input vcf file and pulls out just the scaffold name and position 
# this produces a new output SNP text file for use in the FilterSomeMore.pl script

# usage make_SNP_file.pl infile.vcf



#input .vcf file must be passed to perl script here from the shell script:
my $in = shift(@ARGV); 

#Open the vcf file ($in) and assign it a handle (IN):
open (IN, $in) or die "Could not read the infile = $in\n";

#make sure the regular expression matching below will capture the name of your $in file:
#The purpose of this is to pull out the part of the matched name in parenthes for use in naming the outfile
$in =~ m/output_filtered_files\/([a-zA-Z0-9\_]+)\.vcf/ or die "Failed to match the variant file\n";

#Create and name your future outfile: 
open (OUT, "> output_filtered_files\/SNPS_$1.txt") or die "Could not write the outfile\n";


#The <word> carots are a shortcut command to read in the next line of a text file
#so while you've read in the next line of the file whose handle is called IN, do: 
while (<IN>){
	#chomp is a command that removes the newline (\n) from the end of a line:
	chomp;

	#If a row starts with a "#", it's a header row, and thus isn't needed in our output:
	if (m/^\#/){ 
		#print "#was a header\n"; #test to see if everything is working
	}

	#If the row that was read in is NOT a header row, then:
	# double-check this a sequence line, you might need to edit this reg. expr
	# example sequence line: ENA|CAACVG010015778|CAACVG010015778.1	55916  .  T  C  47.4887	
	# in my case, all the beetle contigs from goran ref. genome start with a capital letter 
	elsif (m/(^[A-Z0-9\.\|]+).{0,1}([0-9]+)/){ 

		#the above checks to see if the line starts with a string of A-Z, 0-9, |, or dot
		##the parentheses around ([A-Z0-9-]+) tells perl to save that bit to $1
		#the .{0,1} says you want to match whatever character comes after that string as well
		#this is to account for the utf8 tab delimiters in vcf files
		#most text is encoded as ascII so you could FIND tabs using [\t]. Not vcfs. 
		#finally, the ([0-9]+) matches the string of digits AFTER the character that came after the first string anchored (^) to the beginning of that line and saves that string of digits as $2...

		#print "matching worked $1 $2";
		print OUT "$1 $2\n";
	}
	else{
		print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
		
	}
}
close (IN);
close (OUT);

print "Finished using $in to make SNP text file\n";
