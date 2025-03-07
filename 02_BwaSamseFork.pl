#!/usr/bin/perl
#
# run bwa conveting sam to bam on each individual  
#


use Parallel::ForkManager;
my $max = 36;
my $pm = Parallel::ForkManager->new($max);
#defining a variable named $genome which contains the name of the reference genome fasta file: 
my $genome = "input_ref_genome_files/CAACVG01.fasta"; 

FILES:

#for each $file passed to this perl script from the shell script as @ARGV: 
foreach $file (@ARGV){ 
	##fork: 
	$pm->start and next FILES; 
	#MAKE SURE THE MATCHING LINE ACCOUNTS FOR FILE LOCATION IF YOU USED THE FULL PATH!!!
	#REMEMBER TO USE \ TO ESCAPE / AND . CHARACTERS
	#if the file name matches some combination of A-Z, a-z, 0-9, underscore, and/or dots, followed by .fastq, then: 
        if ($file =~ m/input_fastq_files\/([A-Za-z0-9_\.]+)\.fastq/){  
		#save the matched name as a real variable called $ind 
		#perl auto-saved the stuff in square brackets as a temporary variable called $1) 
		#Thus, $ind contains the name of an individual beetle
		
	        $ind = $1; 
    	}
        else {
        	die "Failed to match $file\n"; #if the name doesn't match, something went wrong
	}
	
	#The input for bwa aln is: bwa aln ref.fa short_read.fq > aln_sa.sai
	#The input for bwa samse is: bwa samse ref.fa aln_sa.sai short_read.fq > aln-se.sam

#convert .sai to .sam: 
	system "bwa samse $genome output_sai_files/aln_$ind.sai $file > output_sam_files/aln_$ind.sam\n";     


##system writes everything in the quotes to the command line. 
# Thus, the line above calls bwa and runs samse (converting alignment from sai to sam) using my reference genome ($genome), my sai files that I am making dead sure match my fastq file (aln_$ind.sai), my fastq file ($file) and outputs the .sam version as aln_name-of-butterfly.sam. 

	$pm->finish;
}

$pm->wait_all_children;



