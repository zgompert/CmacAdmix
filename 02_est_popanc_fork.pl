#!/usr/bin/perl
#
# manage running popanc est on multiple treatment group files at once by using the fork manager 
#


use Parallel::ForkManager;
my $max = 24; #maximum number of threads/cores to fork to
my $pm = Parallel::ForkManager->new($max);

FILES:
# for each $fq (fastq file) passed from the shell script, do: 
foreach $popanc_file (@ARGV){  
	$pm->start and next FILES; ## fork
	
	# example popanc file name: 
	#BF_BZ.5.F7.C_d_10_n_3.hdf5   

	$popanc_file =~ m/output_hdf5_files\/([A-Z]+_[A-Z]+\.[A-Z0-9\.]+)_(d_[0-9]+_n_[0-9]+)/ or die "failed to match $popanc_file\n";
		

        #Perl automatically saves everything in parentheses as $1 
        #Here we assign $1 to a "real" variable name
        #so it doesn't get lost or replaced with a new $1:
        
	$id = $1;
	$popanc_params = $2;
	
	$popanc_file2 = $popanc_file;
	$popanc_file2 =~ s/ch1/ch2/;

	##printing to the slurmout file:
	print "Running est_popanc on $id\n"; 

	##running est_popanc on each $admix_file file: 
	system "/uufs/chpc.utah.edu/common/home/u6000989/bin/est_popanc -o output_summary_files/${id}_${popanc_params}_popQ.txt -p popQ -c 0.95 -s 0 -w 1 $popanc_file $popanc_file2\n";
	
	#options in popanc: 
	
	##./est_popanc version 0.1 - 9 December 2014

	##Usage:   estpost [options] infile1.hdf5 infile2.hdf5
	##-o     Outfile [default = postout.txt]
	##-p     Name of parameter to summarize: scale, popQ, samQ
	##-c     Credible interval to calculate [default = 0.95]
	##-b     Number of additinal MCMC samples to discard for burn-in [default = 0]
	##-h     Number of bins for posterior sample histogram [default = 20]
	##-s     Which summary to perform: 0 = posterior estimates and credible intervals
	##                       1 = histogram of posterior samples
	##                       2 = convert to plain text
	##                       3 = MCMC diagnostics
	##-w     Write parameter identification to file, boolean [default = 1]
	##-v     Display estpost software version
	

        $pm->finish;
}

$pm->wait_all_children;



