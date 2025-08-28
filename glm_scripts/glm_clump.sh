#SCRIPT WILL PERFORM CLUMPING (FROM GWAS RESULTS) USING PLINK


PHENO=$1
BFILE="all_maf_geno_mind_unrel"

########################################
module load PLINK/1.9b_6.21-x86_64
#######################################

#Formatting of summary statistics
awk '{if ($8=="ADD") print}' "$PHENO"_ml_gwas_w_"$2"_pc_results.PHENO1.glm.logistic.hybrid > "$PHENO"_ml_gwas_w_"$2"_pc_results.PHENO1.glm.logistic

sed -i '1s/^/#CHROM  POS     SNP     REF     ALT     A1      FIRTH?  TEST    OBS_CT  OR      LOG(OR)_SE      Z_STAT  P       ERRCODE\n/' gwas_results/"$PHENO"_ml_gwas_w_"$2"_pc_results.PHENO1.glm.logistic

plink --out gwas_results/"$PHENO"_ml_w_"$2"_pc \
 --clump-kb 500 \
 --clump-r2 0.1 \
 --clump-p1 0.5 \
 --clump-p2 1 \
 --clump "$PHENO"_ml_gwas_w_"$2"_pc_results.PHENO1.glm.logistic \
 --bfile "$BFILE"
