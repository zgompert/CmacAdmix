#!/usr/bin/perl
#
# run entropy jobs
#
# usage: 
# perl RunEntropyFork.pl files_of_initial_admixture_proportions 
#

use Parallel::ForkManager;
my $max = 20;
my $pm = Parallel::ForkManager->new($max);

#recall that I passed the initial -q value files to this perl script as the only @ARGV
foreach $in (@ARGV){ 
	#match the text of the infile (should be cmac_admix_lda_k followed by a digit followed by .txt):
	$in =~ m/input\_files\/cmac\_admix\_lda\_k(\d)\.txt/ or die "failed to match the ldak file";
	#save the part of the string in parentheses above as $k (auto-saved as $1 in perl):
	$k = $1;
	
	#the line below sets the number of chains. In this case I want 10, so 0 through 9:
	foreach $ch (0..9){
		#fork across the CHAINS (not the ldak* files):
		$pm->start and next; 

		#run entropy (first run, 10,000 burn in, 20,000 steps):
#		system "/uufs/chpc.utah.edu/common/home/u6010038/cmac_admix/09_entropy/buerklelab-mixedploidy-entropy-9c1c9c33b477/entropy -i input_files/morefilter_maxcov48000_mind2.filtered_minCov2x_maxND20_cmac_admix_variants.gl -l 20000 -b 10000 -t 5 -k $k -Q 0 -s 50 -q $in -o output_MCMC_files/cmac_admix_MCMC_k$k.ch$ch.hdf5 -w 1 -m 1 -n 2\n";

		## second entropy run for the sake of speed (small, 1,500 burn in, 4,000 steps): 
		system "/uufs/chpc.utah.edu/common/home/u6010038/cmac_admix/09_entropy/buerklelab-mixedploidy-entropy-9c1c9c33b477/entropy -i input_files/morefilter_maxcov48000_mind2.filtered_minCov2x_maxND20_cmac_admix_variants.gl -l 4000 -b 1500 -t 5 -k $k -Q 0 -s 50 -q $in -o output_short_MCMC_files/cmac_admix_short_MCMC_k$k.ch$ch.hdf5 -w 1 -m 1 -n 2\n";

		# -i is the .gl infile; it is hard-coded here
		# -l is the number of MCMC steps for the sampling iterations, here 20,000 steps
		# -b is the number of burn-in MCMC steps (i.e. the number to discard), here 10,000
		# -t is the amount of MCMC iterations to store; here keeping every 5th step
		# -k is the putative number of demes/populations; this is matched from the ldak* files
		
		# -Q is intra- (0) or inter- (1) specific ancestry
		# in other words, are your individuals from within (0) or between (1) species?
		# here, all my individuals are of the same species, so I chose 0

		# -s is the Dirichlet initialization value; between 20-50 is appropriate, here 50
		# the higher the number you choose here, the more similar your intial values will be		
		# -q is a file with the initial admixture proportions for each individual
		# these are the ldak* files passed to the script in @ARGV (hence $in)
		# these estimates are used to ensure all the chains match for easy concatination
		# if this is not done, label swapping might happen where a pop is named A the first
		# time but B the second time, making stacking of chains hard since pops don't match
		
		# -o is the name of the output file, make sure to note which k and chain for each
		# -w is whether the output should (1) or should not (0) contain pop. allele freqs\

		# -m is whether the input file (the .gl file) is in read counts (0) 
		# or in genotype likelihoods (1). Mine is in genotype likelihood format, so 1
		
		# -n is the ploidy, either a number or a matrix, mine are all diploid, so 2

		$pm->finish;
		}
}


$pm->wait_all_children;



