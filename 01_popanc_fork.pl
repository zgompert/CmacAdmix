#!/usr/bin/perl
#
# manage running popanc on multiple treatment group files at once by using the fork manager 
#


use Parallel::ForkManager;
my $max = 24; #maximum number of threads/cores to fork to
my $pm = Parallel::ForkManager->new($max);

FILES:
# for each $fq (fastq file) passed from the shell script, do: 
foreach $admix_file (@ARGV){  
	$pm->start and next FILES; ## fork
	
	# example admixed file name: 
	#BF_BZ.5.F7.C_for_popanc.txt   
	
	#and example parental file name: 
	#BF_for_BFxBZ_for_popanc.txt
	#
	# matching the beetle ID for use naming the output file: 
	#$fq =~ m/([A-Z\_]+\.[A-Z0-9\.]+)\_[a-zA-Z0-9\_\.]+\.gl/ or die "failed to match $fq\n";

	$admix_file =~ m/input_popanc_files\/(([A-Z]+)_([A-Z]+)\.[A-Z0-9\.]+)_for_popanc\.txt/ or die "failed to match $admix_file\n";
	

        #Perl automatically saves everything in parentheses as $1 
        #Here we assign $1 to a "real" variable name
        #so it doesn't get lost or replaced with a new $1:
        
	$id = $1;
        $parent1 = $2;
	$parent2 = $3;

	##printing to the slurmout file:
	print "Running popanc on $id with parents $parent1 and $parent2\n"; 

	##running popanc on each $admix_file file: 
	#recall my parent files follow the structure: BF_for_BFxBZ_for_popanc.txt
	
	#system "/uufs/chpc.utah.edu/common/home/u6000989/bin/popanc input_popanc_files/${parent1}_for_${parent1}x${parent2}_for_popanc.txt input_popanc_files/${parent2}_for_${parent1}x${parent2}_for_popanc.txt $admix_file -o output_hdf5_files/${id}_d_10_n_3.hdf5 -n 3 -d 10 -m 10000 -b 1000 -t 5 -c 1\n";

	#second run with different parameters: 
	
	#system "/uufs/chpc.utah.edu/common/home/u6000989/bin/popanc input_popanc_files/${parent1}_for_${parent1}x${parent2}_for_popanc.txt input_popanc_files/${parent2}_for_${parent1}x${parent2}_for_popanc.txt $admix_file -o output_hdf5_files/${id}_d_50_n_15.hdf5 -n 15 -d 50 -m 10000 -b 1000 -t 5 -c 1\n";
	
	#ZG run with different parameters: 
	
        system "/uufs/chpc.utah.edu/common/home/u6000989/bin/popanc input_popanc_files/${parent1}_for_${parent1}x${parent2}_for_popanc.txt input_popanc_files/${parent2}_for_${parent1}x${parent2}_for_popanc.txt $admix_file -o output_hdf5_files/${id}_d_1000000_n_15_ch1.hdf5 -n 15 -d 1000000 -m 30000 -b 10000 -t 5 -c 1 -l 1 -u 100000 -a 1000 -s 1\n";
        system "/uufs/chpc.utah.edu/common/home/u6000989/bin/popanc input_popanc_files/${parent1}_for_${parent1}x${parent2}_for_popanc.txt input_popanc_files/${parent2}_for_${parent1}x${parent2}_for_popanc.txt $admix_file -o output_hdf5_files/${id}_d_1000000_n_15_ch2.hdf5 -n 15 -d 1000000 -m 30000 -b 10000 -t 5 -c 1 -l 1 -u 100000 -a 1000 -s 1\n";
	
	#options in popanc: 
	
	##./popanc version 0.1 -- 25 Aug 2015

	##Usage:   popAnc [options] infile_p0 infile_p1 infile_Admix
	##-o HDF5 outfile [default = output.hdf5]
	##-n Maximum number of SNPs in each direction to define a locus for LDA [default = 3]
	##-d Maximum distance in each direction to define a locus for LDA [default = 10]
	##-m Number of MCMC steps for the analysis [default = 10000]
	##-b Discard the first n MCMC samples as a burn-in [default = 1000]
	##-t Thin MCMC samples by recording every nth value [default = 1]
	##-c Scale for correlated beta process; ignored if -s 1 [default = 1]
	##-s Estimate scale parameter for beta process model [default = 0]
	##-q Draw samples from simple beta-binomial for CIs on anc. freqs. [default = 0]
	##-u Upper bound for uniform prior on scale [default = 10]
	##-l Lower bound for uniform prior on scale [default = 0.1]
	##-a Proposal +- for metropolis update of scale [default = 0.005]
	##-z Write local ancestry posterior probabilities [default = 0]
	

        $pm->finish;
}

$pm->wait_all_children;



