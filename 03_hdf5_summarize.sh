#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=mcmc_summarize_cmac_admix
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amy.springer@aggiemail.usu.edu


##this is a batch script that summarizes mcmc output of entropy across chains


#Usage:   estpost.entropy [options] infile1.hdf5 infile2.hdf5
#-o     Outfile [default = postout.txt]
#-p     Name of parameter to summarize, e.g., 'q'
#-c     Credible interval to calculate [default = 0.95]
#-b     Number of additinal MCMC samples to discard for burn-in [default = 0]
#-h     Number of bins for posterior sample histogram [default = 20]
#-s     Which summary to perform: 0 = posterior estimates and credible intervals
#                                 1 = histogram of posterior samples
#                                 2 = convert to plain text
#                                 3 = calculate DIC or WAIC (based on input file)
#                                 ('-p deviance' for DIC and '-p likdata' for WAIC)
#                                 4 = MCMC diagnostics
#-v     Display estpost.entropy software version

#The known parameters are: gprob, zprob, p, q, pi, fst, alpha, gamma, deviance, likdata


#change working directory:
cd /uufs/chpc.utah.edu/common/home/u6010038/cmac_admix/09_entropy


#change the numbers after "in" below to match the k-values you used, i.e. 2 3 4, or 4 5 6, etc.
#in the case of my admixed beetles, I only ran k = 3 because I have exactly 3 parental populations 
#for the values listed below in my variable kvals, do:
for kvals in 3
do 
#first, assign a local variable called FILES to save the string that will call all the 10 chain files:
FILES=output_short_MCMC_files/*k$kvals*hdf5
#use echo to print what's in FILES to make sure it's sucking up all the right files
echo my_files_are_$FILES
echo hello_$kvals
#run entropy for admixture proportion summaries: 
/uufs/chpc.utah.edu/common/home/u6010038/cmac_admix/09_entropy/buerklelab-mixedploidy-entropy-9c1c9c33b477/estpost.entropy -o output_MCMC_summary_files/cmac_admix_MCMC_summaries_q_k$kvals -p q -s 0 $FILES 
done
 


#run entropy for genotype summaries:
#genotype summaries are taken across k to account for uncertainty in k:  
#again, in this case (admixed beetles), I am only running k = 3 because I had 3 parental pops
/uufs/chpc.utah.edu/common/home/u6010038/cmac_admix/09_entropy/buerklelab-mixedploidy-entropy-9c1c9c33b477/estpost.entropy -o output_MCMC_summary_files/cmac_admix_MCMC_summaries_gprob_k3 -p gprob -s 0 output_short_MCMC_files/*hdf5

 
