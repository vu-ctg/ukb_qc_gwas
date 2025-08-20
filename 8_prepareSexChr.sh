# Steps and scripts to prepare imputed sex chromosome genotype files and integrate with autosomes
# Imputed sex chromosomes (X and XY) were released added to CTG pipeline in February 2019


### Prep
# create list 'ids_excl_XY_phasing.txt' from UKB field noting samples excluded from phasing (these don't appear to be included in bgen files anyway)

#Link Sex Chromosome genotype data from download dir
ln -s /download/genotype_release_2/imp/*_chrX*.sample .
ln -s /download/genotype_release_2/imp/decrypted/*_chrX*.txt .
ln -s /download/genotype_release_2/imp/decrypted/*_chrX*.bgen .


### Step 1 - convert bgen chromosome files to plink
mkdir qc
cd qc
plink2 --bgen ../_004_ukb_imp_chrX_v3.bgen --sample ../ukbXXXXX_imp_chrX_v3_sXXXXXX.sample --geno .05 --maf 0.0001 --remove ids_excl_XY_phasing.txt --hard-call-threshold 0.1 --make-bed --out chrX_raw
plink2 --bgen ../_004_ukb_imp_chrXY_v3.bgen --sample ../ukbXXXXX_imp_chrXY_v3_sXXXXXX.sample --geno .05 --maf 0.0001 --remove ids_excl_XY_phasing.txt --hard-call-threshold 0.1 --make-bed --memory 60000 --out chrXY_raw


### Step 2 - check QC metrics (maf, hwe, missingness)
# For X chr, these metrics should be compared between males and females and large differences removed
plink2 --bfile chrX_raw --freq --missing --out chrX_snpqc --memory 2000000
plink2 --bfile chrX_raw --hardy --out chrX_snpqc_hwe_MF
plink2 --bfile chrX_raw --freq --geno-counts --missing --hardy --keep-females --out chrX_snpqc_female --memory 2000000
plink2 --bfile chrX_raw --freq --geno-counts --missing --keep-males --out chrX_snpqc_male --memory 2000000

plink2 --bfile chrXY_raw --freq --missing --hardy --out chrXY_snpqc --memory 60000
plink2 --bfile chrXY_raw --freq --geno-counts --missing --hardy --keep-females --out chrXY_snpqc_female --memory 60000
plink2 --bfile chrXY_raw --freq --geno-counts --missing --keep-males --out chrXY_snpqc_male --memory 60000

# Check basic qc (frequency, missingness, hwe) for males and females
# Variant missingness rates >.05
awk '$5 > .05 {print $2}' chrX_snpqc.vmiss | grep -v "ID" > genolt95_all
awk '$5 > .05 {print $2}' chrX_snpqc_female.vmiss | grep -v "ID" > genolt95_female
awk '$5 > .05 {print $2}' chrX_snpqc_male.vmiss | grep -v "ID" > genolt95_male
awk '$5 > .05 {print $2}' chrXY_snpqc.vmiss | grep -v "ID" >> genolt95_all
awk '$5 > .05 {print $2}' chrXY_snpqc_female.vmiss | grep -v "ID" >> genolt95_female
awk '$5 > .05 {print $2}' chrXY_snpqc_male.vmiss | grep -v "ID" >> genolt95_male
wc -l genolt*

# Subject missingness
# Probably not important to include since call rate is lower overall on X; whole-genome based metrics are better
awk '$5 > .03 {print $2}' chrX_snpqc.smiss > imisslt97
wc -l imisslt97

# HWE in females
awk '$14 < .000000005 {print $2}' chrX_snpqc_hwe_MF.hardy.x | grep -v "ID" > hwe5e6
awk '$10 < .0000005 {print $2}' chrXY_snpqc.hardy | grep -v "ID" >> hwe5e6
wc -l hwe5e6

# Heterozygous calls in male X chrom (should be none)
awk '$6 > 0 {print $2}' chrX_snpqc_male.gcount



### Step 3: Create list of variants to filter
Rscript - << 'END'
require(data.table)
all <- fread("genolt95_all",h=F,sep="\t")
f <- fread("genolt95_female",h=F,sep="\t")
m <- fread("genolt95_male",h=F,sep="\t")
table(m$V1 %in% f$V1)
exclsnpmiss <- c(m$V1,f$V1)
table(duplicated(exclsnpmiss))

exclhwe <- fread("hwe5e6",h=F,sep="\t")

#Test difference in missingness & allele freq
fafx <- fread("chrX_snpqc_female.afreq")
fafxy <- fread("chrXY_snpqc_female.afreq")
mafx <- fread("chrX_snpqc_male.afreq")
mafxy <- fread("chrXY_snpqc_male.afreq")
faf <- rbind(fafx,fafxy)
maf <- rbind(mafx,mafxy)

af <- merge(faf,maf,by=c("#CHROM","ID", "REF", "ALT"))
af$diffmaf <- abs(af$ALT_FREQS.x - af$ALT_FREQS.y)
table(af$diffmaf > .02)                         #122
max(af$diffmaf)                                         #0.049458
exclafdiff <- af$ID[af$diffmaf > .02]

fmissx <- fread("chrX_snpqc_female.vmiss")
mmissx <- fread("chrX_snpqc_male.vmiss")
fmissxy <- fread("chrXY_snpqc_female.vmiss")
mmissxy <- fread("chrXY_snpqc_male.vmiss")
fmiss <- rbind(fmissx,fmissxy)
mmiss <- rbind(mmissx,mmissxy)

miss <- merge(fmiss,mmiss,by=c("#CHROM","ID"))
miss$diffmiss <- abs(miss$F_MISS.x - miss$F_MISS.y)
table(miss$diffmiss > .02)                              # 5329
max(miss$diffmiss)                                      #0.0749576
exclmissdiff <- miss$ID[miss$diffmiss > .02]

write.table(excl,"X_excludesnps_freq-hwe-mfdiffs",quote=F,row.names=F,col.names=F)
END


### Step 4: Create list of samples to filter
# Good samples from post-QC autosomal files
awk '{print $1,$1}' /qc/final/phenotypes/geno_samplevars/geno_ukbv2b_samplevars.txt > ids_keep_ukb_qc
# Ambiguous sex individuals to exclude from contributing sex chr genotypes
awk '$5==0 {print $1,$1}' chrX_raw.fam > ids_ambig_sex



### Step 5: Create filtered files
plink2 --bfile chrX_raw --keep ids_keep_ukb_qc --remove ids_ambig_sex --extract ../../convertBgenToPlink/ukbrsids_maf1e_5_info0_9_uniquepos.snplist --exclude X_excludesnps_freq-hwe-mfdiffs --make-bed --out ../chrX_clean
plink2 --bfile chrXY_raw --keep ids_keep_ukb_qc --remove ids_ambig_sex --extract ../../convertBgenToPlink/ukbrsids_maf1e_5_info0_9_uniquepos.snplist --exclude X_excludesnps_freq-hwe-mfdiffs --make-bed --out ../chrXY_clean
cd ../

### Step 6: Check stats of clean files
# Change snp id format to match autosomes: CHR:POS:A1_A2 whera A1 and A2 are alphabetically ordered
awk '{print $1,$5<$6?$1":"$4":"$5"_"$6:$1":"$4":"$6"_"$5,$3,$4,$5,$6}' chrX_clean.bim > chrX_clean.bim.new
mv chrX_clean.bim.new chrX_clean.bim
awk '{print $1,$5<$6?$1":"$4":"$5"_"$6:$1":"$4":"$6"_"$5,$3,$4,$5,$6}' chrXY_clean.bim > chrXY_clean.bim.new
mv chrXY_clean.bim.new chrXY_clean.bim

# Check for duplicate snpids
echo "Check for duplicate snpids.."
cut -d' ' -f2 chrX_clean.bim | sort | uniq -d
cut -d' ' -f2 chrXY_clean.bim | sort | uniq -d

# Compute allele frequencies and missingness
mkdir stats
plink2 --bfile chrX_clean --freq --missing --out stats/chrX_clean_stats
plink2 --bfile chrXY_clean --freq --missing --out stats/chrXY_clean_stats


# 15 SNPs are all alone on the PAR2 region and cause weird manhattan plots (possibly also indicate problems genotyping this region) - select these to remove (add default in filterOutVariants.R script


### Step 7: Split ancestry populations

for POP in AFR AMR EAS SAS EUR; do
for chr in X XY; do
echo $POP
plink2 --bfile chr${chr}_clean --keep ../splitPopulations/extract_${POP}.plinkids --freq --missing --make-bed --out ${TEMP}/chr${chr}_${POP}
done
done




