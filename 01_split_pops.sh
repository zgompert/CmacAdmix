#!/bin/sh 
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=split_pops_cmac_admix
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amy.springer@aggiemail.usu.edu

#load perl:
module load perl

#change working directory to location of my .gl file:
#cd /uufs/chpc.utah.edu/common/home/u6010038/cmac_admix/10_estpEM/input_gl_file
#re-run with new genome alignmeent: 
cd /uufs/chpc.utah.edu/common/home/u6010038/cmac_admix/10_estpEM/input_gl_file_ng

#run the perl script 01_split_pops.pl (from up one directory) and pass the .gl file to it
perl ../01_split_pops.pl *gl 

#print a message if something is messed up:
early()
{
 echo ' '
 echo ' ############ WARNING:  EARLY TERMINATION ############# '
 echo ' '
}
exit
