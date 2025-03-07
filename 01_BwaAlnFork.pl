#!/usr/bin/perl
#
# run bwa alignment on each individual  
#


use Parallel::ForkManager;
my $max = 36;
my $pm = Parallel::ForkManager->new($max);
#defining a variable named  $genome as your chosen reference genome: 
my $genome = "input_ref_genome_files/CAACVG01.fasta";  


FILES:
##for each .fastq file passed to this perl script from the shell script: 
foreach $file (@ARGV){ 
	##fork: 
	$pm->start and next FILES; 
	##If the next file passed through @ARGV matches any string containing some combination of characters A through Z, 0 through 9, dot, or underscore followed by .fastq, then:
        if ($file =~ m/([A-Z0-9_\.]+)\.fastq/){  
		##save the matched name as a real variable called $ind (perl auto-saved what was matched above in a temporary variable called $1) Thus, $ind contains the name of the beetle whose sequences are going to be aligned to the reference:
	        $ind = $1; 
    	}
        else {
        	die "Failed to match $file\n"; #if the name doesn't match, something went wrong
	}
	
	#The input for bwa aln is: bwa aln ref.fa short_read.fq > aln_sa.sai
	#second alignment: 
       	system "bwa aln -n 5 -l 20 -k 2 -t 1 $genome $file > output_sai_files/aln_$ind.sai\n"; 
	
	##system writes everything in the quotes to the command line. 
	#Thus, this line calls bwa and runs the alignment.
	 
	#-n is the total number of mismatches allowed per sequence. 
	#Since my sequences are about 86 bp long, 5 mismatches/86 bp = 0.058% mismatches
	#default for -n is 0.04 (percent mismatches), so I've set this slightly higher
	#Since all my sequences are almost exactly the same length, it's fine to set a value instead of a percentage for this. If you have variable sequence lengths, a percentage might be necessary
	
	#-t is number of threads, set here to one since I'm parallelizing with fork instead
	
	#-k is the maximum edit distance in the seed, 
	#in other words, the maximum number of point mutations (insertions, deletions, or substitutions) required to turn sequence "a" into sequence "b".
	#For example, "ATTGG" can be turned into "AGGTT" with a minimum of FOUR substitutions
	#In this case, I set k to 2, which means if the seed sequence (which I set to 20 bp) has more than 2 mismatches, it is not considered an alignment and the program moves on
	#-k is DIFFERENT than -n. -n is the maximum differences allowed in the WHOLE SEQUENCE, while -k is the maximum differences allowed in the first seeded subsequence
	
	#-l is "take the first INT subsequence as the seed." 
	#In other words, this variable is specifying seed length (in my case, I specified it as 20). 
	#Thus, bwa aln will look at the first 20 bp of each sequence, and if the mismatch is greater than 2 (as specified in -k) in that first part of the sequence, it will drop it and move on
	
	#-q trims the read down to a smaller piece


	$pm->finish;
}

$pm->wait_all_children;



