#Requires mergedPlinkChunks file from step 3 and "geno_ukbv2b_samplevars.txt" from step #5
#for ancestry assignment and to determine filtered out subjects
SCRIPT_DIR=/qc/scripts/geno_release2/splitPopulations
cd $SCRIPT_DIR

#Create 1kg superpopulation specific file
# only includes samples with ibd info, 1kg pop assignment, no sex discrepancies, no bad ukb qc flags 
Rscript - <<'END'
library(data.table)
samplevarfile <- "/qc/final/phenotypes/geno_samplevars/geno_ukbv2b_samplevars.txt"

d <- fread(samplevarfile)
pops <- unique(d$super_pop) #AFR,AMR,EAS,EUR,SAS ancestry super populations
for(p in pops){
fwrite(d[super_pop==p & relative_exclusion==F,.(f.eid,f.eid)],file=paste0("unrelated_subjects_",p,".plinkids"),sep=" ",col.names=F)
fwrite(d[super_pop==p,.(f.eid,f.eid)],file=paste0("extract_",p,".plinkids"),sep=" ",col.names=F)
}
fwrite(d[relative_exclusion==F,.(f.eid,f.eid)],file=paste0("unrelated_subjects.plinkids"),sep=" ",col.names=F)

END

FINALDIR=/qc/final/genotypes/release2b
TEMP=/scratch/$USER/input
mkdir -p $TEMP

cp /qc/scripts/geno_release2/convertBgenToPlink/createPlinkChunks/mergedPlinkChunks_final.* $TEMP

MERGED_PLINK="${TEMP}/mergedPlinkChunks_final"

for POP in AFR AMR EAS SAS EUR; do
echo $POP

plink --bfile $MERGED_PLINK --update-ids ../sampleQC/ukb_linked_sampleid_update.txt --keep ./extract_${POP}.plinkids --freq --missing --make-bed --out ${TEMP}/full_${POP}

mkdir -p ${FINALDIR}/${POP}/full_$POP
cp $TEMP/full_${POP}.* ${FINALDIR}/${POP}/full_$POP

plink --bfile $TEMP/full_${POP} --keep-allele-order --maf 0.01 --keep ./unrelated_subjects_${POP}.plinkids --extract ../assignAncestry/ukb_1kg-1500-150-r10.prune.in --make-bed --freq --missing --out $TEMP/pruned_${POP}

plink --bfile $TEMP/full_${POP} --keep-allele-order --maf 0.01  --extract ../assignAncestry/ukb_1kg-1500-150-r10.prune.in --make-bed --freq --missing --out $TEMP/pruned_all_${POP}


mkdir -p ${FINALDIR}/${POP}/pruned_r2_0_1_maf_0_01_unrel_$POP
mkdir -p ${FINALDIR}/${POP}/pruned_r2_0_1_maf_0_01_all_$POP
cp $TEMP/pruned_all_${POP}.* ${FINALDIR}/${POP}/pruned_r2_0_1_maf_0_01_all_$POP
cp $TEMP/pruned_${POP}.* ${FINALDIR}/${POP}/pruned_r2_0_1_maf_0_01_unrel_$POP

done
