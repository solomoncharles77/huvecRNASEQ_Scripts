#!/bin/bash

#SBATCH --job-name=gctaAssoc
#SBATCH --account=vasccell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cs806@leicester.ac.uk
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --export=NONE

cd $SLURM_SUBMIT_DIR


module load gcta
module load plink

# gcta --bfile ../genotypeData/rawGeno/newImpute2025/exprGeno --make-grm-bin --out genoFiles/vsmcGRM
# gcta --reml --grm-bin genoFiles/vsmcGRM --pheno phenoFiles/rnaEditingIndices_Normalized_AssocReady.txt --mpheno 2 --out heritability/A2G_heritability

# HUVEC Random AEI GWAS with GCTA
gcta --mlma \
--bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--pheno phenoFiles/huvecRandomAEI_clean_AssocReady.txt \
--mpheno 2 \
--out assocRes/huvecRandomAEIGWAS \
--qcovar covFiles/huvec_Sex_3gpc_covariates.txt \
--thread-num 8

# clump with plink
plink --bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--clump assocRes/huvecRandomAEIGWAS_suggestiveHits.txt --clump-p1 0.00005 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/huvecRandomAEIGWAS_clumped_suggestive_loci

################

# HUVEC RefSeq AEI GWAS with GCTA
gcta --mlma \
--bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--pheno phenoFiles/huvecRefSeqAEI_clean_AssocReady.txt \
--mpheno 2 \
--out assocRes/huvecRefSeqAEIGWAS \
--qcovar covFiles/huvec_Sex_3gpc_covariates.txt \
--thread-num 8

# clump with plink
plink --bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--clump assocRes/huvecRefSeqAEIGWAS_suggestiveHits.txt --clump-p1 0.00005 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/huvecRefSeqAEIGWAS_clumped_suggestive_loci

########################

# HUVEC MMSites AEI GWAS with GCTA
gcta --mlma \
--bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--pheno phenoFiles/huvecMMSitesAEI_clean_AssocReady.txt \
--mpheno 2 \
--out assocRes/huvecMMSitesAEIGWAS \
--qcovar covFiles/huvec_Sex_3gpc_covariates.txt \
--thread-num 8

# clump with plink
plink --bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--clump assocRes/huvecMMSitesAEIGWAS_suggestiveHits.txt --clump-p1 0.00005 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/huvecMMSitesAEIGWAS_clumped_suggestive_loci

###################################################################
###################################################################

# HUVEC Mean Editing Level GWAS with GCTA
gcta --mlma \
--bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--pheno phenoFiles/huvecMeanEditingLevel_clean_AssocReady.txt \
--mpheno 1 \
--out assocRes/huvecMELgwas \
--qcovar covFiles/huvec_Sex_3gpc_covariates.txt \
--thread-num 8

##################################################

# HUVEC Mean Editing Level GWAS with GCTA
gcta --mlma \
--bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--pheno phenoFiles/huvecCodingSequenceEditingLevel_clean_AssocReady.txt \
--mpheno 1 \
--out assocRes/huvecMcELgwas \
--qcovar covFiles/huvec_Sex_3gpc_covariates.txt \
--thread-num 8

