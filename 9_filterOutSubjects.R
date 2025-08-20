# Saves a list of all subject IDs for exclusion from analysis that are:
#1) related (kinship coefficient >0.04 based on UKB provided info)
#2) discordant for self-reported and genetic sex
#3) have no genetic data
#4) recommended by UKBIO to exclude due to poor heterozygosity/missingness
#5) withdrawn from UKBIO project (assumes environmental variable 'WITHDRAWNSUBJECT_FILE' defining the file of withdrawn subject IDs provided by UKB)
# ancestry-specific filtering happens at the level of the plink files, which are already split into 5 super_pop groups


#File name of output file with subject IDs to exclude
outFile <- args[1]

#File name of UKB subjec QC information
phenotypeRFile <- "/qc/final/phenotypes/geno_samplevars/geno_ukbv2b_samplevars_poppcs.txt"

#Default output filename
if (is.na(outFile)){
  outFile <- "defaultSubjectExclusions.txt"
}

#Create default
message("Creating list of exclusion subjects...")
message("output file: ", outFile)
stopifnot(!is.na(outFile))

# UKBIO variables used for exclusion criteria
#f.eid: subject id
#f.31.0.0: self-reported sex
#f.22001.0.0: genetic sex 
#f.22010.0.0: recommended genomic analysis exclusions (1=poor heterzygosity/missing)
#f.22011.0.0-4: genetic relatedness pairing (if any)
#f.22012.0.0-4: genetic relatedness factor (estimated kinship coefficient [approx 0.5 x values in GRM matrix])
#relative_exclusion: subjects excluded until no 3rd degree or higher relatives remaining in dataset, defined in earlier QC pipeline step
### data qc steps removed sex/heterozygosity/missingness/ancestry/relatedness abnormalities from genetic analysis dataset
requiredVariables <- c("f.eid","relative_exclusion")

stopifnot(file.exists(phenotypeRFile))
message(phenotypeRFile)
message("Load phenotype data..")

bd <- fread(phenotypeRFile)
stopifnot(all(requiredVariables %in% names(bd)))

message("Filter out withdrawn subjects...")

withdrawnSubjects <- fread(Sys.getenv("WITHDRAWNSUBJECT_FILE"))[[1]]
withdrawnSubjects <- c(withdrawnSubjects,bd[f.eid<0,.(f.eid)]$f.eid)

message(paste(length(withdrawnSubjects),"withdrawn subjects removed"))

message("Filter out subjects that did not surivive geno sample QC...")

noGoodGenoAvailable <- bd[is.na(relative_exclusion),.(f.eid)]$f.eid

#Remove relationship > kinship coefficient threshold
message("Filter out related subjects...")

relatedSubjectIDs<-bd[relative_exclusion==T,.(f.eid)]$f.eid

message(length(noGoodGenoAvailable)," subjects that did not survive sample QC excluded (inconsistent sex, excess heterozygosity, missingness etc..)")
message(paste(length(relatedSubjectIDs)," related subjects excluded"))


excludedSubjects <- unique(c(noGoodGenoAvailable,withdrawnSubjects,relatedSubjectIDs))
message(length(excludedSubjects)," subjects filtered out")
message(nrow(bd)-length(excludedSubjects)," subjects remain (no ancestry selection applied yet)")
message("write to ",outFile)
d <- data.table(FID=excludedSubjects,IID=excludedSubjects)
fwrite(d,file=outFile,row.names=F,col.names=F,quote=F,sep=" ")

message("Done")
