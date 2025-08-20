# Identify quality control indicators for samples such as relatedness and poor genotyping quality (mostly from information provided by UKB). Generate lists of relatives including twins and trios for use in specific analyses.

Rscript - <<'END'
library(data.table)
impdownloaddir <- "/download/genotype_release_2/imp"
caldownloaddir <- "/download/genotype_release_2/cal"

# files for called/imputed genotypes and QC info downloaded from UKB; filenames will be project-specific
samplefile <- file.path(impdownloaddir,"ukbXXXX_imp_chr22_v3_sXXXXXX.sample")
samplecalfile <- file.path(impdownloaddir,"ukbXXXXX_cal_chr22_v2_sXXXXXX.fam")
relfile <- file.path(impdownloaddir,"ukbXXXXX_rel_chr22_sXXXXXX.dat")
samplecalqcfile <- file.path(caldownloaddir,"_ukb_sqc_v2.txt")
ancestryFile <- file.path("../assignAncestry/output","kgpop_assignment.csv")
outputFile <- "/qc/final/phenotypes/geno_samplevars/geno_ukbv2b_samplevars.txt"

#prepare 1kg ancestry assignment file
anc <- fread(ancestryFile)
anc[,id:=NULL]

#prepare relatedness pairs ukb provided
rel_excl_file="rel_exclusions.RData"

if(!file.exists(rel_excl_file)){

rel <- fread(relfile)
rel <- rel[Kinship>0.04,]

relids <- unique(c(rel$ID1,rel$ID2))

message("Filter out related subjects based on kinship coefficient >0.04")   #0.06 is cousins but quite a few 0.04 are reported as well

# Rank individuals for exclusion based on highest number of relatives in the dataset
relations <- sapply(relids,function(i)nrow(rel[ID1==i | ID2==i,]))

relfrq <- data.table(id=relids,nrel=relations)
setorder(relfrq,-nrel)

relleft <- rel
excludesubj <- c()
message(nrow(relleft))
for(i in 1:nrow(relfrq))
{
if(i %% 1000 ==0) message(i," ",nrow(relleft)," ",length(excludesubj))
prevrels <- nrow(relleft)
relleft <- relleft[ID1 != relfrq[i,id] & ID2 != relfrq[i,id],]
if(prevrels>nrow(relleft)) excludesubj <- c(excludesubj,relfrq[i,id])
}
save(excludesubj, file=rel_excl_file)
}else{
load(rel_excl_file)
}

#link ranked sample ids to ukb sample ids and save to file to use in --update-ids plink option
impsubj <- fread(samplefile)
calsubj <- fread(samplecalfile)

impsubj[,old_FID:=0:(nrow(impsubj)-1)]
impsubj[,old_IID:=old_FID]
impsubj[,missing:=NULL]

impsubj <- impsubj[-1,]
setcolorder(impsubj,c("old_FID","old_IID","ID_1","ID_2","sex"))
fwrite(impsubj,file="ukb_linked_sampleid_update.txt",sep=" ")

#sample qc cal file
samplecalqc <- fread(samplecalqcfile)
setnames(samplecalqc,c("sampleid","id","array","batch",
"platename",
"well",
"cluster_CR",
"dQC",
"internal_pico",
"reported_gender",
"inferred_gender",
"x_int",
"y_int",
"subm_platename",
"subm_well",
"sample_qc_miss_rate",
"heterozygosity",
"het_pc_corrected",
"het_missing_outliers",
"sex_aneuploidy",
"in_kinship_table",
"excluded_kinship_inference",
"excess_relatives",
"in_white_brit_subset",
"used_in_pc_comp",
paste0("PC_ukb",1:40),
"in_phasing_input_autosome",
"in_phasing_input_X",
"in_phasing_input_XY"))

samplecalqc[,calrankid:= 1:nrow(samplecalqc)]
setnames(calsubj,c("FID","IID","P1","P2","sex","batch"))
calsubj[,calrankid:= 1:nrow(calsubj)]

stopifnot(all(nrow(samplecalqc)==nrow(calsubj)))
d <- merge(calsubj,samplecalqc,by=c("calrankid","batch"))
message("Merge cal with imp subjects...")
d <- merge(impsubj,d,by.x="ID_1",by.y="FID")
message(nrow(d), " nr of subjects to start with")
message("Merge with ancestry info...")
d <- merge(d,anc,by.x="old_FID",by.y="imprankid",x.all=T)

setnames(d,old="ID_1",new="f.eid")
head(d)
d <- d[,-c(1,3:5,7:10)]
setnames(d,old="sex.y",new="sex")
d[,relative_exclusion:=f.eid %in% excludesubj]



#Apply sample filters

subjqclog <- "sampleqcreport.csv"
logsteps <- data.table(step=character(),removed=numeric(),left=numeric())

remove_subjects <- function(dat,label,condition){
res <- list(label,sum(condition),nrow(dat)-sum(condition))
list(log=res,d=d[!condition,])
}
message("Before sample qc ",nrow(d), "subjects")

res <- remove_subjects(d,"inconsistent_gender",with(d,reported_gender!=inferred_gender))
d <- res$d
logsteps <- rbind(logsteps,res$log)
message("After inconsistent sex removal ", nrow(d), " subjects remaining")

res <- remove_subjects(d,"sex_aneuploidy",with(d,sex_aneuploidy==1))
d <- res$d
logsteps <- rbind(logsteps,res$log)
message("After aneuploidy removal ", nrow(d), " subjects remaining")

res <- remove_subjects(d,"het_missing_outliers",with(d,het_missing_outliers==1))
d <- res$d
logsteps <- rbind(logsteps,res$log)
message("After het_miss_outlier removal ", nrow(d), " subjects remaining")

res <- remove_subjects(d,"excluded_in_kinship_inference",with(d,excluded_kinship_inference==1))
d <- res$d
logsteps <- rbind(logsteps,res$log)
message("After excluded in kinship inf removal ", nrow(d), " subjects remaining")

res <- remove_subjects(d,"no_assigned_pop",with(d,is.na(d$super_pop)))
d <- res$d
logsteps <- rbind(logsteps,res$log)
message("After no pop assignment removal ", nrow(d), " subjects remaining")

res <- remove_subjects(d,"excess_relatives",with(d,excess_relatives==1))
d <- res$d
logsteps <- rbind(logsteps,res$log)
message("After excess rel removal ", nrow(d), " subjects remaining")

fwrite(logsteps,file=subjqclog)
fwrite(d,file=outputFile,sep=" ")



## Examine types of relationships in UKB based on derived genetic relatedness values
## & create plink-style fam files for trios and MZs to use in genetic analyses

rel <- fread(relfile)
rel <- rel[Kinship > 0,]
table(round(rel$Kinship,1))


rel$degree <- integer()
rel$type <- character()
# set relationships based on KING paper recommendations http://people.virginia.edu/~wc9c/KING/manual.html
rel[Kinship >= .354,'degree'] <- 0
rel[Kinship >= .177 & Kinship < .354,'degree'] <- 1
rel[Kinship >= .0884 & Kinship < .177,'degree'] <- 2
rel[Kinship >= .0442 & Kinship < .0884,'degree'] <- 3
table(rel$degree)

rel[degree==0,'type'] <- "MZ"
rel[IBS0 < .001 & degree==1,'type'] <- "PO"
rel[IBS0 > .001 & degree==1,'type'] <- "FS"
rel[IBS0 > .001 & degree==2,'type'] <- "2nd"
rel[IBS0 > .001 & degree==3,'type'] <- "3rd"
table(rel$type)


pdf("ukb_relatives.pdf")
plot(rel$IBS0, rel$Kinship, type="p",
  col = rel$degree+1, pch=as.numeric(factor(rel$type)),
  main = "UKB relative pairs",
  xlab="Estimated Pr(IBS=0)", ylab = "Estimated Kinship Coefficient")
dev.off()

fwrite(rel,"ukb_relative_pair_types.txt",sep="\t")



# Create subsets of relations of interest to annotate in fam files
po <- rel[type=="PO",]
mz <- rel[type=="MZ",]

table(duplicated(mz$ID1))
table(duplicated(mz$ID2))       #all mzs are paired together, no dups
table(mz$ID1 %in% c(po$ID1,po$ID2))
table(mz$ID2 %in% c(po$ID1,po$ID2))     #7 pairs have parent/child

table(duplicated(po$ID1))
table(duplicated(po$ID2))

mz$FID <- paste(mz$ID1,mz$ID2,sep="_")
mz$PID <- paste(mz$ID1,mz$ID2,"1",sep="_")
mz$MID <- paste(mz$ID1,mz$ID2,"2",sep="_")
mz.fam <- data.frame(FID=rep(mz$FID,2),IID=c(mz$ID1,mz$ID2),PID=rep(mz$PID,2),
                        MID=rep(mz$MID,2),sex=0,pheno=-9)

# trios have index individual that must be related to two other individuals
## (but this can also include one parent with two children)
## according to Bycroft 2019 article, there are no 3-generation families in UKB
## should be 1,066 trios from 1029 parents, 37 quartets (possibly slightly less due to withdrawal)
allids <- c(po$ID1,po$ID2)
dupids <- allids[duplicated(allids)]
dup1 <- po[po$ID1 %in% dupids,]
dup2 <- po[po$ID2 %in% dupids,c(2,1,3:7)]
names(dup2)[1:2] <- c("ID1","ID2")              #order dataset so that each dup individual is an index
dups <- rbind(dup1,dup2)
dups$pairid <- paste(apply(dups[,1:2],1,min),apply(dups[,1:2],1,max),sep="_")
table(duplicated(dups$pairid))                  #pairids duplicated when both individuals have multiple relatives


dems <- fread("/gwa_covs.txt")       #created covariate file to add age/sex info to infer parents
dems$sex <- dems$sex-1
dems$sex[dems$sex==0] <- 2              #switch m/f coding (ukb is 1=f/2=m)
dems1 <- dems[,2:4]
names(dems1) <- c("ID1","sex1","age1")
dems2 <- dems1
names(dems2) <- c("ID2","sex2","age2")
dups <- merge(dups,dems2,by="ID2",all.x=T)
dups <- merge(dups,dems1,by="ID1",all.x=T)

reorder <- dups[dups$age1 > dups$age2,c(2,1,3:8,11,12,9,10)]            #reorder so child is index id
names(reorder) <- names(dups)
dups[dups$age1 > dups$age2,] <- reorder
dups$rel1 <- dups$rel2 <- character()           #check that reorder worked
for (i in 1:dim(dups)[1]) {
        if (dups$age1[i] < dups$age2[i]) { dups$rel1[i] <- "child" }
        if (dups$age2[i] < dups$age1[i]) { dups$rel2[i] <- "child" }
        if ((dups$age1[i] > dups$age2[i]) & dups$sex1[i]==2) { dups$rel1[i] <- "mother" }
        if ((dups$age2[i] > dups$age1[i]) & dups$sex2[i]==2) { dups$rel2[i] <- "mother" }
        if ((dups$age1[i] > dups$age2[i]) & dups$sex1[i]==1)  { dups$rel1[i] <- "father" }
        if ((dups$age2[i] > dups$age1[i]) & dups$sex2[i]==1)  { dups$rel2[i] <- "father" }
}
table(dups$rel1, dups$rel2)

dups <- dups[!duplicated(dups$pairid),]         #remove redundant pairs (exactly the same P-O ids)
##!!#! Note: two individuals are offspring of an MZ parent so there appear to be three parents. Difficult to resolve which one is the actual parent (but doesn't matter genetically - so remove based on clues in the data e.g. more similar home environment).

table(dups$ID1 %in% dups$ID1[duplicated(dups$ID1)])
trio <- dups[dups$ID1 %in% dups$ID1[duplicated(dups$ID1)],]             #remove non-dup children (e.g. children tagged by dup ids of parents that had two children in the dataset)
table(table(trio$ID1))                          #should be 2 of each ID
table(table(trio$ID2))
table(trio$rel2)

#Create ped file for trios
fa <- trio[trio$rel2=="father",]
mo <- trio[trio$rel2=="mother",]
table(mo$ID1 %in% fa$ID1)
table(fa$ID1 %in% mo$ID1)
table(mo$ID2[mo$ID1 %in% fa$ID1==F] %in% c(mz$ID1,mz$ID2))      #dups w/2 mo mostly MZs, otherwise error? remove
mo <- mo[mo$ID1 %in% fa$ID1==T,]
mo.o <- mo[,c("ID1","ID2","sex1")]
fa.o <- fa[,c("ID1","ID2")]
names(mo.o) <- c("IID","MID","SEX")
names(fa.o) <- c("IID","PID")
out <- merge(mo.o,fa.o,by="IID")
dim(out)                                        #should be the same as mo/fa objects!
out$FID <- paste(out$IID,out$PID,out$MID,sep="_")
out$PHENO <- -9
out <- out[,c("FID","IID","PID","MID","SEX","PHENO")]
fwrite(out,"trios.fam",sep="\t",col.names=F)

#Adjust ped file for mzs
mz.fam$FID <- as.character(mz.fam$FID)
mz.fam$PID <- as.character(mz.fam$PID)
mz.fam$MID <- as.character(mz.fam$MID)
table(mz.fam$IID %in% out$IID)
for (i in mz.fam$IID[which(mz.fam$IID %in% out$IID)]) {
        mz.fam[mz.fam$IID==i,c("FID","PID","MID")] <- out[out$IID==i,c("FID","PID","MID")]
}
for (i in mz.fam$IID) {
        mz.fam$sex[mz.fam$IID==i] <- dems$sex[dems$IID==i]
}

fwrite(mz.fam,"mz_twins.fam",sep="\t",col.names=F)


END

