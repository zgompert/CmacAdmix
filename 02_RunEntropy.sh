#!/bin/sh 
#SBATCH --time=140:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=entropy_cmac_admix
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amy.springer@aggiemail.usu.edu

#load modules necessary to run entropy:
module load gsl
module load hdf5
module load perl

#change working directory:
cd /uufs/chpc.utah.edu/common/home/u6010038/cmac_admix/09_entropy

#run the perl script RunEntropyFork.pl and pass all the ldak* files to it 
#not passing the .gl file here because I don't want to screw up the forking; will just hard-code it
perl 02_RunEntropyFork.pl input_files/cmac_admix_lda_*.txt 

#print a message if something is messed up:
early()
{
 echo ' '
 echo ' ############ WARNING:  EARLY TERMINATION ############# '
 echo ' '
}
exit
