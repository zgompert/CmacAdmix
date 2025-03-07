#!/usr/bin/perl
#
# manage running estpEM on multiple split .gl files at once by using the fork manager 
#


use Parallel::ForkManager;
my $max = 24; #maximum number of threads/cores to fork to
my $pm = Parallel::ForkManager->new($max);

FILES:
# for each $fq (fastq file) passed from the shell script, do: 
foreach $fq (@ARGV){  
	$pm->start and next FILES; ## fork
	
	# example file name: 
	#CA.3.F19.L_morefilter_maxcov48000_mind2.filtered_minCov2x_maxND20_cmac_admix_variants.gl
	# matching the beetle ID for use naming the output file: 
	#$fq =~ m/([A-Z\_]+\.[A-Z0-9\.]+)\_[a-zA-Z0-9\_\.]+\.gl/ or die "failed to match $fq\n";

	#NEED TO ALTER MATCHING FOR THE NEW PREFIX BRIAN ADDED FOR DENOTING WHO RAN THE EXPERIMENT (AS = A SPRINGER, BK = BRIAN KISSMER)
	$fq =~ m/(AS\.[A-Z\_]+\.[A-Z0-9\.]+)\_[a-zA-Z0-9\_\.]+\.gl/ or die "failed to match $fq\n";
	

        #Perl automatically saves everything in parentheses as $1 
        #Here we assign $1 to a "real" variable name
        #so it doesn't get lost or replaced with a new $1:
        
	$id = $1;
        
	##printing to the slurmout file:
	print "Running estpEM on $id\n"; 

	##running estpEM version 0.1 on each $fq file: 
	## recall this is being run from output_split_gl_files
	system "/uufs/chpc.utah.edu/common/home/u6000989/bin/estpEM -i $fq -h 2 -o ../output_allele_freq_files/$id\_estpEM_allele_freqs.txt\n";
	
	#options in estpEM: 
	
	#-i Infile with genetic data for the population
	#-o Outfile for allele frequency estimates [default = out.txt]
	#-e Tolerance for EM convergence [default = 0.001]
	#-m Maximum number of EM iterations [default = 20]
	#-h Number of header lines [default = 0]
	

        $pm->finish;
}

$pm->wait_all_children;



