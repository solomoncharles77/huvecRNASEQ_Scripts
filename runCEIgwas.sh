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

# HUVEC Random CEI GWAS with GCTA
gcta --mlma \
--bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--pheno phenoFiles/huvecRandomCEI_clean_AssocReady.txt \
--mpheno 2 \
--out assocRes/huvecRandomCEIGWAS \
--qcovar covFiles/huvec_Sex_3gpc_covariates.txt \
--thread-num 8

# clump with plink
plink --bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--clump assocRes/huvecRandomCEIGWAS_suggestiveHits.txt --clump-p1 0.00005 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/huvecRandomCEIGWAS_clumped_suggestive_loci

################

# HUVEC RefSeq CEI GWAS with GCTA
gcta --mlma \
--bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--pheno phenoFiles/huvecRefSeqCEI_clean_AssocReady.txt \
--mpheno 2 \
--out assocRes/huvecRefSeqCEIGWAS \
--qcovar covFiles/huvec_Sex_3gpc_covariates.txt \
--thread-num 8

# clump with plink
plink --bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--clump assocRes/huvecRefSeqCEIGWAS_suggestiveHits.txt --clump-p1 0.00005 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/huvecRefSeqCEIGWAS_clumped_suggestive_loci

########################

# HUVEC MMSites CEI GWAS with GCTA
gcta --mlma \
--bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--pheno phenoFiles/huvecMMSitesCEI_clean_AssocReady.txt \
--mpheno 2 \
--out assocRes/huvecMMSitesCEIGWAS \
--qcovar covFiles/huvec_Sex_3gpc_covariates.txt \
--thread-num 8

# clump with plink
plink --bfile ../genotypeData/rawGeno/newImpute2025/huvecGeno \
--clump assocRes/huvecMMSitesCEIGWAS_suggestiveHits.txt --clump-p1 0.00005 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/huvecMMSitesCEIGWAS_clumped_suggestive_loci