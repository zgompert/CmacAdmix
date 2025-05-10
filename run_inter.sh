#!/bin/bash 
#SBATCH --time=140:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=braker
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

source ~/.bashrc

module load interproscan

cd /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/t_crist_hyw154_green_h2/braker

cat braker.aa | perl -p -i -e 's/\*//g' > braker_aa.fasta

/uufs/chpc.utah.edu/sys/installdir/r8/interproscan/5.63/interproscan.sh -cpu 48 -i braker_aa.fasta -goterms -b interpro_aa
