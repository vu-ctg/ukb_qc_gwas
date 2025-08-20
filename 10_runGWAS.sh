# Example plink input file used to run GWAS on UKB data that has been quality controlled using the prior pipeline steps.
# On Snellius, this analysis is run across 1000 genomic chunks (filtered by QC and split by ancestry group) with automated parallelization and resource optimization
# On UKB-RAP, this analysis is run per chromosome with input lists of QC-filtered SNPs generated locally and ancestry-specific individuals
# Assumes files definining the phenotypes, covariates, variant and subject exclusions have already been created

plink2 --bfile <chunk-prefix> 
  --out <chunk-prefix> \
  --remove defaultSubjectExclusions.txt \
  --exclude defaultVariantExclusions.txt \
  --pheno plink_phenos_cont.txt \
  --covar plink_covars.txt \
  --covar-name sex,array01,age,pop_pc1-pop_pc20 \
  --maf 0.0001 \
  --geno 0.05 \
  --freq \
  --ci 0.95 \
  --covar-variance-standardize \
  --glm hide-covar no-x-sex \
  --memory 1700 
