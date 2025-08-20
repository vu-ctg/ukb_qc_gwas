# Compute within-ancestry principal components to use as analysis covariates

# STEP 1: compute population pcs
# Script computePopPCs.sh
SCRIPT_DIR=/qc/scripts/geno_release2/computePopPCs
TEMP_INPUT=/scratch/$USER/input
mkdir -p $TEMP_INPUT
cd $TEMP_INPUT

#Run flashpca2
NR_PCS=30
MAX_MEM=100000
NPROCS=47

#Project PCs to UKB sample
#cd flaspca_1kg
for POP in AFR EUR EAS SAS AMR
do
echo $POP
PRUNED_FILE=/qc/final/genotypes/release2b/$POP/pruned_r2_0_1_maf_0_01_all_$POP/pruned_all_$POP
cp ${PRUNED_FILE}.* $TEMP_INPUT
SUFFIX=_${POP}_flashpcs
/programs/flashpca2/flashpca --bfile pruned_all_$POP -n $NPROCS -d $NR_PCS -m $MAX_MEM --suffix ${SUFFIX} > log${SUFFIX}
cp *${SUFFIX}* $SCRIPT_DIR
done
echo DONE


# STEP 2: add derived population PCS to the covariate file geno_ukbv2b_samplevars.txt
# Script integratePopPCs.sh

Rscript - <<'END'
library(data.table)
genovarprefix="/qc/final/phenotypes/geno_samplevars/geno_ukbv2b_samplevars"
inputfile <- paste0(genovarprefix,".txt")
stopifnot(file.exists(inputfile))

d  <- fread(inputfile)
pops <- unique(d$super_pop)
npcs <- 30
pcnames <- paste0("pop_pc",1:npcs)

#initialize pop_pcs
d[,eval(pcnames) := as.numeric(NA)]
setkey(d,f.eid)

for (p in pops)
{
poppcs <- fread(paste0("pcs_",p,"_flashpcs"))

d2<- merge(d,poppcs,by.x="f.eid",by.y="FID",all.x=T)
setkey(d,f.eid)
stopifnot(d2$f.eid==d$f.eid)
for(i in 1:npcs)
{
d[super_pop==p,eval(paste0("pop_pc",i)) := d2[super_pop==p,eval(paste0("PC",i)),with=F]]
stopifnot(all(d2[d2$super_pop==p,eval(paste0("PC",i)),with=F] == d[d$super_pop==p,eval(paste0("pop_pc",i)),with=F],na.rm=T))
}
}
fwrite(d,file=paste0(genovarprefix,"_poppcs.txt"),sep=" ",na="NA")
END
