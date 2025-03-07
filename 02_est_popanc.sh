#!/bin/sh 
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=popanc_est_cmac_admix
#SBATCH --mail-user=zach.gompert@usu.edu

#presumably purge all previously opened modules
#then load perl, gcc, and hdf5:
module purge
module load gcc/8.5.0 hdf5/1.10.7
module load perl


#change working directory to location of my perl file:
cd /uufs/chpc.utah.edu/common/home/gompert-group5/projects/cmac_admix/12_popanc



#run the perl script 02_est_popanc_fork.pl and pass all .hdf5 files from the output_popanc_files directory to it:

#perl 02_est_popanc_fork.pl output_hdf5_files/*d_10*hdf5
#perl 02_est_popanc_fork.pl output_hdf5_files/*d_1000*hdf5 
perl 02_est_popanc_fork.pl output_hdf5_files/*d_1000000_n_15_ch1*

