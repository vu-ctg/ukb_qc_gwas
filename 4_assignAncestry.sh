#STEP 1a - preprocess 1kg genotype data
# Script preprocess_1kg.sh

### Prep ###

#Download 1kG reference files
#Set file locations, including file with ukb_maf and ukb_info for filtered SNPs (maf>1e-5 and info>=0.9) 
SNP_INFO_FILE=/qc/scripts/geno_release2/convertBgenToPlink/filteredsnps_maf1e_5_info0_9.txt
SCRIPT_DIR=/qc/scripts/geno_release2/assignAncestry
UKB_PLINK_DIR=/qc/scripts/geno_release2/convertBgenToPlink/createPlinkChunks
KG_CHR_DIR=/datasets/1kgenomes/phase3v5/intermediate
TEMP_INPUT=/scratch/$USER/input
TEMP_OUTPUT=/scratch/$USER/output

cd $SCRIPT_DIR

module load R/3.3.1
module load plink2/1.90b3y

mkdir -p genos_ukb
mkdir -p genos_1kg
mkdir -p merged_ukb_1kg
mkdir -p $TEMP_INPUT
mkdir -p $TEMP_OUTPUT

# # The following files should be located in main directory:
# "1KGPpopulations.key" - plain text file with 2 columns listing 1kG subpopulation names and associated superpopulation
# "LongRangeLDregions_Price2008_hg19_20141005.txt" - plain text file listing range of high LD regions in plink set format

### Extract and reformat genotype files to merge UKB and 1kG datasets ###
# Create list of called SNPs in UKB sample from imputed files
echo "Create allchr.info.."
awk '{print $4 " " $1 " " $5 " " $6 " " $7 " " $8 " " $10 " " $2}' $SNP_INFO_FILE > allchr.info
cd $SCRIPT_DIR

echo "Perform snp selection for merge with 1000 genomes and ancesetry assignment (proper SNPs with INFO=1)..."
Rscript - <<'END' 
require(data.table)
called <- fread("allchr.info")
#Only keep SNPs with info=1, and A,C,T,G alleles (ie. remove indels since they are much more difficult to reliably match to other files)
called <- subset(called,INFO_UKB==1 & nchar(A1_UKB)==1 & nchar(A2_UKB)==1)
called[, CHR_POS := paste0(CHR,":",POS_HG19)]
#Remove indels

setcolorder(called,c("RSID_UKB","CHR","POS_HG19","A1_UKB","A2_UKB","MAF_UKB","INFO_UKB","SNPID_UNIQUE","CHR_POS"))
fwrite(called,"ukb_called_snps",quote=F,row.names=F,sep=" ")
fwrite(called[,.(CHR_POS)],"ukb_called_snps.chr_bp",quote=F,row.names=F,col.names=F)
fwrite(called[,.(RSID_UKB)],"ukb_called_snps.rsids",quote=F,row.names=F,col.names=F)
fwrite(called[,.(SNPID_UNIQUE)],"ukb_called_snps.plinkid",quote=F,row.names=F,col.names=F)
END

## Change 1kG SNP IDs to match format of ukb
echo "Preprocess 1kg files..."
awk '{print $2 "\t" $3 "\t" $3 "\t" $9}' ukb_called_snps > ukb_called_snps.range
sed -i '1d' ukb_called_snps.range

cd genos_1kg
echo "Extract selected ukb snps from 1000 genome files..."
# Extract called ukb snps (do this in parallel on compute node)
for chr in {1..22}; do
(
cd $SCRIPT_DIR/genos_1kg
plink --memory 2000 --allow-no-sex --bfile ${KG_CHR_DIR}/1kg3v5_chr${chr} --extract range ../ukb_called_snps.range --freq --make-bed --out 1kg3v5_chr${chr}_filtered
)&
done
wait

#Rename snpid to match plink ids of UKB
cd $SCRIPT_DIR/genos_1kg
echo "Rename snpid to plink id in 1000 genome files..."
for chr in {1..22}; do
mv 1kg3v5_chr${chr}_filtered.bim 1kg3v5_chr${chr}_filtered.temp
awk '{alleleOrder= ($5<$6 ? $5"_"$6 : $6 "_" $5); print $1 "\t" $1 ":" $4 ":" alleleOrder "\t" $3 "\t" $4 "\t" $5 "\t" $6}' 1kg3v5_chr${chr}_filtered.temp > 1kg3v5_chr${chr}_filtered.bim
done

wc -l 1kg3v5_chr*_filtered.bim
head 1kg3v5_chr1_filtered.bim
rm *.temp

cat *.bim > ../all_1kg_snps.bim
cd ..

#Identify overlappping snps between UKB and 1kg (to reduce the number of extracted SNPs from UKB before processing)
echo "Identify overlapping snps UKB 1kg..."
Rscript - <<'END'
require(data.table)
kgsnps <- fread("all_1kg_snps.bim",header=F)
ukb_called_snps <- fread("ukb_called_snps.plinkid",sep=" ",header=F)
overlap_snps <- intersect(kgsnps$V2,ukb_called_snps$V1)
message(length(overlap_snps), " snps overlap between 1kg and ukb snp selection")

write(overlap_snps,"overlapping_ukb_1kg_snps.plinkid")
END

#Merge 1kg chr files
echo "Merge 1kg chr files..."
rm kg_filelist
for chr in {1..22}; do
echo "./genos_1kg/1kg3v5_chr${chr}_filtered" >> kg_filelist
done

plink --merge-list kg_filelist --allow-no-sex --make-bed --out $TEMP_INPUT/kg_file
cp $TEMP_INPUT/kg_file.* ./genos_1kg





#!#!#!#!#!

#STEP 1b - Preprocess UKB genotype data: merge and prune to prepare for PCA analysis with smartpca
# Script preprocess_ukb.sh 

#Extract UKB called SNPs from mergedPlink file
cd $SCRIPT_DIR
echo "copy UKB plink fileset to scratch..."
cp ${UKB_PLINK_DIR}/mergedPlinkChunks_final.* $TEMP_INPUT

echo "process UKB plink set..."
plink --bfile $TEMP_INPUT/mergedPlinkChunks_final --extract overlapping_ukb_1kg_snps.plinkid --allow-no-sex --maf 0.001 --hwe 0.000001 --exclude range LongRangeLDregions_Price2008_hg19_20141005.txt --make-bed --out $TEMP_INPUT/ukb_file

echo "copy back UKB 1000genomes SNP subset to cache..."
cp $TEMP_INPUT/ukb_file.* ./genos_ukb


echo "Merge ukb and 1kg file and only retain non-missing snps..."
#merge ukb and 1kg file
plink --bfile $TEMP_INPUT/ukb_file --bmerge ./genos_1kg/kg_file --geno 0 --make-bed  --out $TEMP_INPUT/ukb_1kg_merge
#turned out to be unnecessary - use in case of merging errors:
#plink --bfile $TEMP_INPUT/ukb_1kg_merge --exclude $TEMP_INPUT/ukb_1kg_merge.missnp --allow-no-sex --freq --missing --make-bed --out $TEMP_INPUT/ukb_1kg_all


echo "Copy merged fileset from scratch..."
cp $TEMP_INPUT/ukb_1kg_merge.* ./merged_ukb_1kg

wc -l $TEMP_INPUT/ukb_1kg_merge.bim
wc -l $TEMP_INPUT/ukb_1kg_merge.fam


#!#!#!#!#!

#STEP 2 - Compute PCs in 1kg and project on UKB (using smartpca) 
# Script project_1kgpcs_on_ukb.sh 

cd $SCRIPT_DIR
module load eb
module load EIGENSOFT/7.2.1-foss-2016b
mkdir -p pca

# Make smartPCA ped file with ancestry groups
Rscript - <<'END'
require(data.table)
fam <- fread("./merged_ukb_1kg/ukb_1kg_all_pruned.fam")
pop <- fread("/datasets/1kgenomes/phase3v5/final/1kg3v5.ALL.panel",h=F)
out <- merge(fam[,1:5],pop[,c(1,3)],by="V1",all.x=T,sort=F)

stopifnot(all(fam[["V1"]]==out[["V1"]]))   # Make sure order is not changed

out[is.na(out$V3.y),V3.y:= "UKB"]
table(out$V3.y)
fwrite(out,"ukb_1kg_all_pruned.pedind",quote=F,row.names=F,col.names=F,sep=" ")
END

cp  ./merged_ukb_1kg/ukb_1kg_all_pruned.* $TEMP_INPUT

# Create the smartPCA parameter text file in pca folder using (uncommented) text below:
cd pca

awk '{print $1}' ../1KGPpopulations.key > ../1KG_populations.list

cat << 'END' > ukb_1kg_smartpca.param
genotypename: ${TEMP_INPUT}/ukb_1kg_all_pruned.bed
snpname: ${TEMP_INPUT}/ukb_1kg_all_pruned.bim
indivname: ../ukb_1kg_all_pruned.pedind
evecoutname: ukb_1kg.evec
evaloutname: ukb_1kg.eval
altnormstyle: YES
numoutevec: 30
numoutlieriter: 5
numthreads: 15
poplistname: ../1KG_populations.list
snpweightoutname: ukb_1kg.snpweights
END

# Run smartPCA
cat ukb_1kg_smartpca.param
#/programs/eigensoft/EIG-6.1.3/bin/smartpca -p ukb_1kg_smartpca.param
smartpca -p ukb_1kg_smartpca.param
                                           

#!#!#!#!#!

#STEP 3 - Assign ancestry to UKB subjects
# Script assign_ancestry.sh

cd $SCRIPT_DIR
mkdir -p ibd
mkdir -p output

Rscript - <<'END'
# Calculate mahalanobis distance for each sample and assign to best fitting population
require(data.table)
require(ggplot2)
evec <- fread("./pca/ukb_1kg.evec",skip=1)
names(evec) <- c("id",paste0("pc",1:30),"pop")
head(evec)
dim(evec)

ukb <- evec[pop=="UKB"]
evec <- evec[!pop=="UKB"]
evec$pop <- factor(evec$pop)
pops <- paste(levels(evec$pop))

for (i in pops) {
  print(i)
  eval(parse(text=paste0( i, " <- evec[pop=='", i, "']")))
  eval(parse(text=paste0( "m.", i, " <- colMeans(", i, "[,2:31])")))
  eval(parse(text=paste0( "c.", i, " <- cov(", i, "[,2:31])")))
  eval(parse(text=paste0( "ukb$", i, " <- mahalanobis(ukb[,2:31],m.", i, ",c.",i,")")))
}
head(ukb)
ukb$minMahal <- apply(ukb[,..pops],1,min)
ukb$popMatch <- pops[apply(ukb[,..pops],1,which.min)]
table(ukb$popMatch)

# Identify outliers who were +/-6SD from best matching subpopulation
ukb$mahalSD <- NA
for (i in pops) {
  ukb$mahalSD[ukb$popMatch==i] <- scale(ukb[popMatch==i,minMahal])
  }  
capture.output(table(round(ukb$mahalSD)>6,ukb$popMatch),file="./output/ukb_pop_assignment_outliers.txt")

key <- read.table("1KGPpopulations.key",strings=F)
names(key) <- c("pop","super_pop")
out <- merge(ukb,key,by.x="popMatch",by.y="pop",all.x=T)
names(out)[1] <- "sub_pop"

capture.output(table(round(out$mahalSD)>6,out$super_pop),file="./output/ukb_superpop_assignment_outliers.txt")
capture.output(table(out$super_pop),file="./output/ukb_superpopulations.txt")
capture.output(table(out$sub_pop,out$super_pop),file="./output/ukb_populationmatrix.txt")

ggsave("./output/ukb_1kg_mahalSD_from_assigned_superpop.png", plot = 
         ggplot(data = out, aes(x=super_pop, y=mahalSD)) + geom_boxplot(aes(color=sub_pop)))


anc <- out
anc <- anc[anc$mahalSD<=6,]

# Check against other measures of ancestry & qc metrics
#anc.ukb <- fread("anc_vars.txt")
#out$f.eid <- as.numeric(substr(out$id,1,7))
#require(car)
#anc.ukb$ethnicity <- recode(anc.ukb$f.21000.0.0,"c(1,1001,1002,1003)='White';c(2,2001,2002,2003,2004)='Mixed';c(3,3001,3002,3003)='SAS';c(3004,5)='EAS';c(4,4001,4002,4003)='Black';6='Other';else=NA")
#anc <- merge(out,anc.ukb,by="f.eid",all.x=T)
#table(anc$super_pop,anc$f.21000.0.0)
#table(anc$super_pop,anc$f.22006.0.0)
#table(anc$super_pop,anc$ethnicity)
#write.table(table(anc$super_pop,anc$ethnicity),"./output/selfrpt_ethn_vs_assigned_superpop.table")
#table(is.na(anc$f.22006.0.0),anc$super_pop)   # 25861 additional EUR, 2027 AFR, 893 AMR, 605 EAS, 2577 SAS
# correlation between UKB provided PCs and our own
#cor(anc$pc1,anc$f.22009.0.1,use="p")  # 0.9786144
#cor(anc$pc2,anc$f.22009.0.2,use="p")  # 0.6513453
#cor(anc$pc3,anc$f.22009.0.3,use="p")  # 0.7123841
#cor(anc$pc4,anc$f.22009.0.4,use="p")  # 0.1922447
#cor(anc$pc5,anc$f.22009.0.5,use="p")  # 0.2932173
#cor(anc$pc6,anc$f.22009.0.5,use="p")  # 0.0006709


# Visualize principal components 
ggsave("./output/ukb_1kgPC1PC2_superpop.png", plot=
         ggplot(data=anc, aes(pc1, pc2)) + geom_point(aes(color=super_pop, shape=super_pop)))
ggsave("./output/ukb_1kgPC3PC4_superpop.png", plot=
         ggplot(data=anc, aes(pc3, pc4)) + geom_point(aes(color=super_pop, shape=super_pop)))
ggsave("./output/ukb_1kgPC5PC6_superpop.png", plot=
         ggplot(data=anc, aes(pc5, pc6)) + geom_point(aes(color=super_pop, shape=super_pop)))

ggsave("./output/ukb_1kgPC1PC2_facet_superpop.png", plot=
         ggplot(data=anc, aes(pc1, pc2)) + geom_point(aes(color=super_pop, shape=super_pop)) + facet_wrap(~super_pop))
ggsave("./output/ukb_1kgPC3PC4_facet_superpop.png", plot=
         ggplot(data=anc, aes(pc3, pc4)) + geom_point(aes(color=super_pop, shape=super_pop)) + facet_wrap(~super_pop))
ggsave("./output/ukb_1kgPC5PC6_facet_superpop.png", plot=
         ggplot(data=anc, aes(pc5, pc6)) + geom_point(aes(color=super_pop, shape=super_pop)) + facet_wrap(~super_pop))

anclong <- data.frame(super_pop=anc$super_pop,pc=rep(1:6,each=(dim(anc)[1])),value=c(anc$pc1,anc$pc2,anc$pc3,anc$pc4,anc$pc5,anc$pc6))
ggsave("./output/ukb_1kgPC1-6_superpop_hist.png", plot=ggplot(anclong, aes(value))+geom_density(aes(color=super_pop))+facet_wrap(~pc,nrow=2))

eur <- anc[anc$super_pop=="EUR",]
dim(eur)
ggsave("ukb_1kgPC1PC2_subpop_EUR.png", plot=
         ggplot(data=eur, aes(pc1, pc2)) + geom_point(aes(color=sub_pop, shape=sub_pop)))
ggsave("ukb_1kgPC3PC4_subpop_EUR.png", plot=
         ggplot(data=eur, aes(pc3, pc4)) + geom_point(aes(color=sub_pop, shape=sub_pop)))
ggsave("ukb_1kgPC5PC6_subpop_EUR.png", plot=
         ggplot(data=eur, aes(pc5, pc6)) + geom_point(aes(color=sub_pop, shape=sub_pop)))

# Compare with PC plots for the 1kG reference genomes
evec <- merge(evec,key,by.x="pop",by.y="sub_pop",all.x=T)
evec$sub_pop <- evec$pop
table(evec$super_pop)
table(evec$sub_pop,evec$super_pop)
ggsave("./output/1kgPC1PC2_superpop.png", plot=
         ggplot(data=evec, aes(pc1, pc2)) + geom_point(aes(color=pop)))
ggsave("./output/1kgPC3PC4_superpop.png", plot=
         ggplot(data=evec, aes(pc3, pc4)) + geom_point(aes(color=pop)))
ggsave("./output/1kgPC5PC6_superpop.png", plot=
         ggplot(data=evec, aes(pc5, pc6)) + geom_point(aes(color=pop)))

eur <- evec[evec$super_pop=="EUR",]
dim(eur)
ggsave("1kgPC1PC2_subpop_EUR.png", plot=
         ggplot(data=eur, aes(pc1, pc2)) + geom_point(aes(color=sub_pop, shape=sub_pop)))
ggsave("1kgPC3PC4_subpop_EUR.png", plot=
         ggplot(data=eur, aes(pc3, pc4)) + geom_point(aes(color=sub_pop, shape=sub_pop)))
ggsave("1kgPC5PC6_subpop_EUR.png", plot=
         ggplot(data=eur, aes(pc5, pc6)) + geom_point(aes(color=sub_pop, shape=sub_pop)))

anc[,imprankid:=sapply(strsplit(anc$id,split=":"),first)]
fwrite(anc,file="./output/kgpop_assignment.csv")

END

