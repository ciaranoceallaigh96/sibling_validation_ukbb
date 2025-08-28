#Recalculated PC generation

#Load PLINK
########################################
module load PLINK/2.00-alpha2-x86_64_avx2

##########################################
BFILE="all_maf_geno_mind"


plink2 --bfile "$BFILE" \
       --extract UKBB_pruned_2.prune.in \
       --pca 50 approx \
       --threads 8 \
       --exclude range extended_hla.txt \
       --out UKBB_pca


