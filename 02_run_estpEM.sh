#!/bin/sh 
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=estpEM_cmac_admix
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amy.springer@aggiemail.usu.edu

#load perl:
module load perl

#change working directory to location of my .gl file:
cd /uufs/chpc.utah.edu/common/home/u6010038/cmac_admix/10_estpEM/output_split_gl_files

#run the perl script 02_fork_estpEM.pl (from up one directory) and pass all the split .gl files to it
perl ../02_estpEM_fork.pl *gl 

#print a message if something is messed up:
early()
{
 echo ' '
 echo ' ############ WARNING:  EARLY TERMINATION ############# '
 echo ' '
}
exit
