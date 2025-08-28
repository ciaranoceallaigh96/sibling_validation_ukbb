#SCRIPT WILL PERFORM GWAS (with variable number of PCs) USING PLINK
#removes any individuals in the test sets
#
#
#Runs PLINK-GLM GWAS with Covariates on FULL train set

PHENO=$1

echo "Running GWASs for phenotype: $PHENO"


########################################
module load PLINK/2.00-alpha2-x86_64_avx2
##########################################

BFILE="all_maf_geno_mind_unrel"
REMOVALS="test_${PHENO}_sibs_nonsibs.ids.no_pc"


#--keep-males for prostate cancer
#--keep-females for breast cancer


#0 PCs
plink2 --bfile "$BFILE" --pheno "$PHENO".full.pheno --remove "$REMOVALS" --threads 24 --exclude range extended_hla.txt --glm allow-no-covars --out "$PHENO"_ml_gwas_w_0_pc_results

for i in 2 4 6 8 10 12 14 16; do
  end=$((i + 2))
  cut -f 1-"$end" ukbb_sex_age_pc1_16.txt.subset > ukbb_pc1_"$i".txt.subset
  plink2 --bfile "$BFILE" --pheno "$PHENO".full.pheno --remove "$REMOVALS" --covar ukbb_pc1_"$i".txt.subset --threads 24 --exclude range extended_hla.txt --glm --out "$PHENO"_ml_gwas_w_"$i"_pc_results
done


