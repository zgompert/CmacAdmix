#!/bin/sh 
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=run_varne_cmac_admix
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amy.mo.springer@gmail.com

#load perl:
module load perl

#change working directory to location of my perl script:
cd /uufs/chpc.utah.edu/common/home/gompert-group5/projects/cmac_admix/11_estimating_var_Ne

#run the perl script 02_run_varne.pl and pass FINAL allele freq files to it:

#running varne for gen 1  to gen 18/19/20 only for all pops and reps:
perl 03_run_varne.pl input_allele_freq_files/AS.*.F[12][890].*txt

#running varne for gen 1 to gen 5/6/7 for all pops and reps: 
#perl 03_run_varne.pl input_allele_freq_files/AS.*.F[567].*txt

#running varne for gen 5/6/7 to gen 18/19/20 for all pops and reps:
#perl 03_run_varne.pl input_allele_freq_files/AS.*_*.F[12][890].*txt

#print a message if something is messed up:
early()
{
 echo ' '
 echo ' ############ WARNING:  EARLY TERMINATION ############# '
 echo ' '
}
exit
