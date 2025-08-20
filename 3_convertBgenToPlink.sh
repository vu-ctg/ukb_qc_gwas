#STEP 0
#Make sure to move a mergedPhenoData.sample (a bgen .sample file) with the correct number of subjects into the directory where bgen (chunk) files are stored

#!#!#!#!# 

#STEP 1
#To allow merging all plink files it is important to remove variants with non-unique position,also some SNP filtering is applied (e.g. maf and info, missing).
#In particular snps are included if:
#*maf_ukb >= 0.00001
#*info_ukb >= 0.9
#*if multiple snps have the same chr position, take the one with the highest info score (multiple can remain if they have the same INFO value)
#*missing geno < 0.05 after hard call [--geno 0.05 plink option]
#In all these files unique_id is the plink id used within CTG in the format chr:pos:A1_A2 where A1 and A2 are in lexicographical order. 

#Script performSNPselection.sh:

#Maf and info per snp files per chromosome downloaded from UKB
UKB_VERSION="v3"
UKB_PREFIX="_004_ukb"

MFI_DIR=/download/genotype_release_2/imp/decrypted
META_INFO_DIR="/qc/final/genotypes/release2b/meta_info/"
mkdir -p $META_INFO_DIR
FILTERED_SNP_FILE=/qc/scripts/geno_release2/convertBgenToPlink/filteredsnps_maf1e_5_info0_9.txt
FULL_SNP_FILE=${META_INFO_DIR}/full_ukb_snplist.txt

#Select snps with MAF>=0.00001 and info>=0.9
#Create unique id (chr:bp_hg19:A1_A2 where A1,A2 in alphabetical order)
echo "Create ${FILTERED_SNP_FILE} and ${FULL_SNP_FILE}"

#Write headers
echo "" | awk 'BEGIN{FS="\t";OFS="\t"}{print "CHR","SNPID_UNIQUE","SNPID_UKB","RSID_UKB","POS_HG19","A1_UKB","A2_UKB","MAF_UKB","MIN_ALLELE_UKB","INFO_UKB"}' - > $FILTERED_SNP_FILE
cp $FILTERED_SNP_FILE $FULL_SNP_FILE

#Write (filtered) SNP info
for CHR in {1..22} X XY; do
    MFI_FILE=${MFI_DIR}/${UKB_PREFIX}_mfi_chr${CHR}_${UKB_VERSION}.txt
    echo $MFI_FILE
    head $MFI_FILE
    awk -v CHR=$CHR 'BEGIN{FS="\t";OFS="\t"}$6>=0.00001 && $8>=0.9 {uniqueID=($4<$5?$4"_"$5:$5"_"$4);print CHR,CHR":"$3":"uniqueID,$0}' $MFI_FILE  >> $FILTERED_SNP_FILE
    awk -v CHR=$CHR 'BEGIN{FS="\t";OFS="\t"}{uniqueID=($4<$5?$4"_"$5:$5"_"$4);print CHR,CHR":"$3":"uniqueID,$0}' ${MFI_FILE}  >> $FULL_SNP_FILE
done


module load R/3.3.1
Rscript - << END
message("Excluding SNPs on same position, variants with the highest info score is selected on each position with duplicates")
message("This step can take a long time")
library(data.table)
inputfile <- "/qc/scripts/geno_release2/convertBgenToPlink/filteredsnps_maf1e_5_info0_9.txt" 
outputfile <- "/qc/scripts/geno_release2/convertBgenToPlink/filteredsnps_maf1e_5_info0_9_uniquepos.txt"
snpoutputfile <- "/qc/scripts/geno_release2/convertBgenToPlink/ukbrsids_maf1e_5_info0_9_uniquepos.snplist"
snpoutputfile2 <- "/qc/scripts/geno_release2/convertBgenToPlink/ukbrsids_maf1e_5_info0_9.snplist"
df <- fread(inputfile)
fwrite(df[,.(RSID_UKB)],file=snpoutputfile2)

subdf<- df[, .SD[which.max(abs(INFO_UKB))], by=.(CHR,POS_HG19)]
fwrite(subdf[,.(RSID_UKB)],file=snpoutputfile)
fwrite(subdf,file=outputfile)
END


#!#!#!#!# 

# STEP 2
# Run plink conversion on bed chunks
# Script startChunkAnalysis.sh with input createPlinkChunks.txt analysis file
plink --bgen <chunk-prefix>.bgen ref-first \
  --sample mergedPhenoData.sample \
  --geno 0.05 --hard-call-threshold 0.1 \
  --extract ukbrsids_maf1e_5_info0_9_uniquepos.snplist \
  --make-bed \
  --memory 60000 --threads 15 \
  --out <chunk-prefix>_filtered


#!#!#!#!# 

# STEP 3:
# Merge all plink chunks into one large plink file set run. The final step of this script will create a file hardcalled_ukb_snplist.txt. This file includes all SNPs in the plink files and a column 'exclude_snp' flagging SNPs that need to be excluded from analysis since they do not meet the filtering criteria.
# Script mergeFullPlinkFiles.sh 

SCRIPT_DIR="/qc/scripts/geno_release2/convertBgenToPlink"
CHUNK_DIR="${SCRIPT_DIR}/createPlinkChunks/results_bgen"
OUTPUT_DIR="${SCRIPT_DIR}/createPlinkChunks"
TEMP_DIR="/scratch/"
mkdir -p $TEMP_DIR
mkdir -p $TEMP_DIR/input
mkdir -p $TEMP_DIR/output
MERGE_FILE="mergedPlinkChunks"

BIM_DIR=$TEMP_DIR/input
MERGE_DIR=$TEMP_DIR/output

echo "Copy chunks to scratch..."
for f in $CHUNK_DIR/*.bim $CHUNK_DIR/*.fam $CHUNK_DIR/*.bed; do
 cp ${f} ${BIM_DIR};
done

echo "Create merge file.."
rm ${MERGE_FILE}.files
for f in $BIM_DIR/*.bim
do
  echo $f | cut -d'.' -f1 >> ${MERGE_FILE}.files
done

head -n3 ${MERGE_FILE}.files

echo "Change snp id into plink snpid (with spaces)"
#Change snp id into plink snpid: CHR:POS:A1_A2 whera A1 and A2 are alphabetically ordered
for f in $BIM_DIR/*.bim
  do
    awk '{print $1,$5<$6?$1":"$4":"$5"_"$6:$1":"$4":"$6"_"$5,$3,$4,$5,$6}' $f > ${f}.new
    mv ${f}.new ${f}
  done

module load plink2

echo "Run plink with --merge-list to merge chunks..."
plink --merge-list ${MERGE_FILE}.files --out $MERGE_DIR/${MERGE_FILE}

#Identify duplicate snpids
echo "Identify duplicate snpids.."
cut -d' ' -f2 ${MERGE_DIR}/${MERGE_FILE}.bim | sort | uniq -d > ${MERGE_DIR}/excluded_duplicate_ids.txt

#Remove duplicate snpids while copying back
echo "Remove duplicate snpids while copying back..."
plink --bfile ${MERGE_DIR}/${MERGE_FILE} --exclude ${MERGE_DIR}/excluded_duplicate_ids.txt --make-bed --out ${MERGE_DIR}/${MERGE_FILE}_final

echo "Compute allele frequencies and missingness"
plink --bfile ${MERGE_DIR}/${MERGE_FILE}_final --freq --missing --out ${MERGE_DIR}/${MERGE_FILE}_final

#Replace tabs with spaces in bim file
echo "Replace tabs with spaces in bim file..."
sed 's/\t/ /g' ${MERGE_DIR/${MERGE_FILE}_final.bim > ${MERGE_DIR}/${MERGE_FILE}.bim.space
mv ${MERGE_DIR}/${MERGE_FILE}_final.bim ${MERGE_DIR}/${MERGE_FILE}.bim.original
mv ${MERGE_DIR}/${MERGE_FILE}.bim.space ${MERGE_DIR}/${MERGE_FILE}_final.bim

echo "Copy merged plink files back..."
cp $MERGE_DIR/${MERGE_FILE}_final.* $OUTPUT_DIR
cp $MERGE_DIR/excluded_duplicate_plinkids.txt $OUTPUT_DIR
cp ${MERGE_FILE}.files $OUTPUT_DIR


#!#!#!#!# 

# STEP 4:
# check whether the number of SNPs and subjects is as expected (.bim and .fam files)
