# Summary of scripts

## Scripts and data for population size

- [popsize.stan](popsize.stan)
- [cmac_admix_lentil_pop_growth_forzach.R](cmac_admix_lentil_pop_growth_forzach.R)
- [Lentil_pop_growth_data.csv](Lentil_pop_growth_data.csv)

## Scripts for processing the raw sequence data

- [parse_barcodes768.pl](parse_barcodes768.pl)
- [PhiXFilterFork.pl](PhiXFilterFork.pl)

## Scripts for alignment, variant calling and filtering

- [01_BwaAlnFork.pl](01_BwaAlnFork.pl)
- [01_vcfFilter.pl](01_vcfFilter.pl)
- [02_BwaSamseFork.pl](02_BwaSamseFork.pl)
- [02_call_variants.sh](02_call_variants.sh)
- [03_check_read_depth.R](03_check_read_depth.R)
- [03_Sam2BamFork.pl](03_Sam2BamFork.pl)
- [04_filterSomeMore.pl](04_filterSomeMore.pl)



## Scripts for genotype and allele frequency estimation

- [01_make_lda_file.R](01_make_lda_file.R)
- [02_estpEM_fork.pl](02_estpEM_fork.pl)
- [01_split_pops.pl](01_split_pops.pl)
- [01_split_pops.sh](01_split_pops.sh)
- [02_make_SNP_file.pl](02_make_SNP_file.pl)
- [02_RunEntropyFork.pl](02_RunEntropyFork.pl)
- [02_RunEntropy.sh](02_RunEntropy.sh)
- [02_run_estpEM.sh](02_run_estpEM.sh)
- [04_remove_headers_for_pca.txt](04_remove_headers_for_pca.txt)
- [mkDpPlots.R](mkDpPlots.R)

## Scripts for ancestry inference

- [01_popanc_fork.pl](01_popanc_fork.pl)
- [01_run_popanc.sh](01_run_popanc.sh)
- [02_est_popanc_fork.pl](02_est_popanc_fork.pl)
- [02_est_popanc.sh](02_est_popanc.sh)
- [plotAncestryEarly.R](plotAncestryEarly.R)
- [plotAncestry.R](plotAncestry.R)
- [allele_freqs_for_popanc2.R](allele_freqs_for_popanc2.R)


## Scripts for Ne, the null model and tests of parallelism

- [03_run_varne.pl](03_run_varne.pl)
- [03_run_varne.sh](03_run_varne.sh)
- [nullModel.R](nullModel.R)




