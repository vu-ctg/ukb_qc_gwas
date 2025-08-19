# ukb_qc_gwas

Author: Jeanne Savage (j.e.savage@vu.nl) & Sven Stringer

This repository describes the analysis pipeline used by the CTG lab for quality control of imputed array data and GWAS in the UK Biobank. Scripts assume that the data can be downloaded locally (no longer possible in most cases since the transition to the online Research Analysis Platform). This pipeline has been used for publications on the UKB data carried out under application #16406.

#Quality control
Order of operations:
1) Download and decrypt bgen and sample files (data field 22828)
2) ['chunkChromosomes'](chunkChromosomes) Split bgen files into smaller chunks (10k SNPs) for parallel processing 
3) ['convertBgenToPlink'](convertBgenToPlink) Filter SNPs and convert bgen chunk files to plink hardcall format  
- MAF >= 0.00001
- INFO score  >= 0.9
- Missing geno < 0.05 after hard call 
- If multiple snps have the same chr position, take the one with the highest info score (multiple can remain if they have the same INFO value). Unfortunately the --bgen option in plink does not work well with variants that have multiple variants on the same position. This resulted in ~60K variants that should have been filtered not being filtered, which will be excluded in a later step.
4) ['assignAncestry'](assignAncestry) Project ancestry information from 1000 Genomes to assign UKB individuals to continental ancestry superpopulations and subgroups
5) ['sampleQC'](sampleQC) Identify quality control indicators for samples such as relatedness and poor genotyping quality
6) ['splitPopulations.sh'](splitPopulations.sh) Split full sample into continental ancestry superpopulations and including vs. excluding relatives 
7) ['computePCs.sh'](computePCs.sh) Compute within-ancestry principal components to use as analysis covariates
8) ['prepareSexChr'](prepareSexChr)
