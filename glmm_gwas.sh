#SCRIPT WILL PERFORM glmm-GWAS (with covariates) USING GCTA
#removes any individuals in the test case
#


PHENO=$1

echo "Running GLMM GWAS for phenotype: $PHENO"

########################################
#module load GCTA/1.94.1-gfbf-2023a
##########################################

BFILE="all_maf_geno_mind_unrel"
REMOVALS="test_${PHENO}_sibs_nonsibs.ids.no_pc"
SPARSE="./GRMs/${PHENO}_grm_sparse05"
PHENOTYPE_FILE="${PHENO}.pheno"

#Keep single sex for breast and prostate cancer

GRMs/gcta-1.94.4-linux-kernel-3-x86_64/gcta64 --bfile all_maf_geno_mind_unrel \
       --fastGWA-mlm-binary \
       --grm-sparse "$SPARSE" \
       --pheno "$PHENOTYPE_FILE"  \
       --threads 8 \
       --out  "$PHENO"_fastgwa_16pc \
       --qcovar ukbb_pc1_16.txt.subset \
       --exclude  hla_snps.txt \
       --remove "$REMOVALS"
