#SCRIPT WILL PERFORM CLUMPING (FROM GWAS RESULTS) USING PLINK

PHENO=$1
BFILE="all_maf_geno_mind_unrel"

########################################
module load PLINK/1.9b_6.21-x86_64
#######################################

plink --out "$PHENO"_lmm_w_16_pc \
 --clump-kb 500 \
 --clump-r2 0.1 \
 --clump-p1 0.5 \
 --clump-p2 1 \
 --clump "$PHENO"_fastgwa_16pc.fastGWA \
 --bfile "$BFILE"

