#
#SCRIPT WILL PERFORM PRS USING PLINK on both SIBS and NON-SIBS
#

echo "Current date and time: $(date)"


PHENO=$1
NUM_SNPS=$2
PC=$3

module load PLINK/2.00-alpha2-x86_64_avx2

echo "Running PRS for phenotype: $PHENO"
echo "Number of SNPs to include: $NUM_SNPS"

tail -n +2  "$PHENO"_ml_w_"$PC"_pc.clumped | awk '{print $3}' | head -n $NUM_SNPS >"$PHENO"_ml_w_"$PC"_pc.clumped.top."$NUM_SNPS".snps



SNPS_INCLUDED2=$(wc -l < "$PHENO"_ml_w_"$PC"_pc.clumped.top."$NUM_SNPS".snps)

echo "Using A1 column as effect allele."
echo "Converting OR to ln(OR) i.e. natural log"

#Reformat for PLINK --score input
awk 'NR==FNR { snps[$1]; next } $3 in snps {print $3, $6, log($10)}' "$PHENO"_ml_w_"$PC"_pc.clumped.top."$NUM_SNPS".snps "$PHENO"_ml_gwas_w_"$PC"_pc_results.PHENO1.glm.logistic > gwas_results/"$PHENO"_ml_w_"$PC"_pc.clumped.top."$NUM_SNPS".prs



#TEST ON SIBS
plink2 --bfile all_maf_geno_mind \
 --pheno "$PHENO".full.pheno \
 --score "$PHENO"_ml_w_"$PC"_pc.clumped.top."$NUM_SNPS".prs 1 2 3 \
 --out sibs_prs_w_"$PC"_pc_"$PHENO"_"$NUM_SNPS" \
 --allow-no-sex \
 --keep individuals_in_pairs_with_"$PHENO".txt


#TEST ON NON-SIBS
plink2 --bfile all_maf_geno_mind \
 --pheno "$PHENO".full.pheno \
 --score "$PHENO"_ml_w_"$PC"_pc.clumped.top."$NUM_SNPS".prs 1 2 3 \
 --out non_sibs_prs_w_"$PC"_pc_"$PHENO"_"$NUM_SNPS" \
 --allow-no-sex \
 --keep "$PHENO"_included_individuals_age.txt


sed -i 's/#FID/FID/g' prs_results/non_sibs_prs_w_"$PC"_pc_"$PHENO"_"$NUM_SNPS".sscore
sed -i 's/#FID/FID/g' prs_results/sibs_prs_w_"$PC"_pc_"$PHENO"_"$NUM_SNPS".sscore


