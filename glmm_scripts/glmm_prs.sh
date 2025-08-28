#SCRIPT WILL PERFORM PRS USING PLINK on both SIBS and NON-SIBS
#

echo "Current date and time: $(date)"

PHENO=$1
NUM_SNPS=$2
PC=$3

module load PLINK/2.00-alpha2-x86_64_avx2

echo "Running PRS for phenotype: $PHENO"
echo "Number of SNPs to include: $NUM_SNPS"

tail -n +2  "$PHENO"_lmm_w_"$PC"_pc.clumped | awk '{print $3}' | head -n $NUM_SNPS > "$PHENO"_lmm_w_"$PC"_pc.clumped.top."$NUM_SNPS".snps


SNPS_INCLUDED2=$(wc -l < "$PHENO"_lmm_w_"$PC"_pc.clumped.top."$NUM_SNPS".snps)

echo "Using A1 column as effect allele."
echo "No need to convert from OR to ln(OR) as --fastGWA-mlm-binary already reports log(OR) as BETA"

#Reformat for PLINK --score input [SNP A1 ln(OR)]
awk 'NR==FNR { snps[$1]; next } $2 in snps {print $2, $4, $11}' "$PHENO"_lmm_w_"$PC"_pc.clumped.top."$NUM_SNPS".snps "$PHENO"_fastgwa_16pc.fastGWA > lmm_results/"$PHENO"_lmm_w_"$PC"_pc.clumped.top."$NUM_SNPS".snps.prs




#TEST ON SIBS
plink2 --bfile all_maf_geno_mind \
 --pheno "$PHENO".full.pheno \
 --score "$PHENO"_lmm_w_"$PC"_pc.clumped.top."$NUM_SNPS".snps.prs 1 2 3 \
 --out sibs_prs_lmm_"$PC"_pc_"$PHENO"_"$NUM_SNPS" \
 --allow-no-sex \
 --keep individuals_in_pairs_with_"$PHENO".txt
 

#TEST ON NON-SIBS
plink2 --bfile all_maf_geno_mind \
 --pheno ."$PHENO".full.pheno \
 --score "$PHENO"_lmm_w_"$PC"_pc.clumped.top."$NUM_SNPS".snps.prs 1 2 3 \
 --out non_sibs_prs_lmm_"$PC"_pc_"$PHENO"_"$NUM_SNPS" \
 --allow-no-sex \
 --keep "$PHENO"_included_individuals_age.txt

 ######################################################################

sed -i 's/#FID/FID/g' non_sibs_prs_lmm_"$PC"_pc_"$PHENO"_"$NUM_SNPS".sscore
sed -i 's/#FID/FID/g' sibs_prs_lmm_"$PC"_pc_"$PHENO"_"$NUM_SNPS".sscore


