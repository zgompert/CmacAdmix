#!/bin/sh
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=BCF_call_haydenii
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=amy.springer@aggiemail.usu.edu

#load bcftools module:
module load bcftools/1.9

#change working directory:
cd /uufs/chpc.utah.edu/common/home/u6010038/cmac_admix/06_variant_calling

#run bcftools to call variants: 
#BCFTOOLs call version 1.9 options used

bcftools call -v -c -f GQ -p 0.01 -O v -o output_vcf_file/cmac_admix_variants.vcf output_bcf_file/cmac_admix.bcf  
#-v option lets you output variant sites only (thus, the variant calling bit: only keeps the parts of each sequence that are variable sites
#-c/m = the original calling method (conflicts with -m) or alternative model for multiallelic and rare-variant calling (conflicts with -c)
#f = format to output: GQ or GP
#p = p-value threshold: when using -c, accept variants with a p-value less than value specified
#P = the prior for expected substitution rate (use bigger for more sensitivity). Use only with -m
#O = output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' (here it is 'v') 
#o = write output to the following file name [file name]



