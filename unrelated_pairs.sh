#Creates non-sibling pairs based on sibling pair characteristics

PHENO=$1
FAM_FILE="all_maf_geno_mind_unrel.fam"
SIBLING_FILE="unique_filtered_sibling_pairs_with_${PHENO}.csv"
OG_FAM_FILE="all_maf_geno_mind_unrel.fam"

#Number of sibling pairs
n_pairs=$(tail -n +2 "$SIBLING_FILE" | wc -l)

#Extract IDs from OG_FAM_FILE
awk '{print $2}' "$OG_FAM_FILE" > og_fam_ids.txt

#Separate cases and controls
cp ../phenotype_files/"$PHENO".ids cases.txt
awk '{print $2, $5}' "$FAM_FILE" > sex_info.txt
awk '{print $1}' "$FAM_FILE" > all_ids.txt
grep -vFf cases.txt all_ids.txt > controls.txt

#Filter cases and controls by OG_FAM_FILE presence
grep -Ff og_fam_ids.txt cases.txt > cases_in_og_fam.txt
grep -Ff og_fam_ids.txt controls.txt > controls_in_og_fam.txt

#Filter cases and controls by sex
if [[ "$PHENO" == "prostate_cancer_self_report" ]]; then
    grep -Ff cases_in_og_fam.txt sex_info.txt | awk '$2 == 1 {print $1}' > filtered_cases.txt
    grep -Ff controls_in_og_fam.txt sex_info.txt | awk '$2 == 1 {print $1}' > filtered_controls.txt
elif [[ "$PHENO" == "breast_cancer_self_report" ]]; then
    grep -Ff cases_in_og_fam.txt sex_info.txt | awk '$2 == 2 {print $1}' > filtered_cases.txt
    grep -Ff controls_in_og_fam.txt sex_info.txt | awk '$2 == 2 {print $1}' > filtered_controls.txt
else
    cp cases_in_og_fam.txt filtered_cases.txt
    cp controls_in_og_fam.txt filtered_controls.txt
fi

#Extract sibling IDs
awk -F, '{print $1; print $2}' "$SIBLING_FILE" | sort | uniq | grep -Ff og_fam_ids.txt > sibling_ids.txt

#Exclude sibling IDs
grep -vFf sibling_ids.txt filtered_cases.txt > cases_without_siblings.txt
grep -vFf sibling_ids.txt filtered_controls.txt > controls_without_siblings.txt

#Calculate available number of cases and controls
num_cases=$(wc -l < cases_without_siblings.txt)
num_controls=$(wc -l < controls_without_siblings.txt)
max_pairs=$((num_cases<num_controls ? num_cases : num_controls))

oversample_factor=10
target_sample_size=$((n_pairs * oversample_factor))
target_sample_size=$((target_sample_size < max_pairs ? target_sample_size : max_pairs))

#Randomly sample larger pools
shuf -n "$target_sample_size" --random-source=<(yes 42) cases_without_siblings.txt > oversample_cases.txt
shuf -n "$target_sample_size" --random-source=<(yes 42) controls_without_siblings.txt > oversample_controls.txt

#combine oversampled pairs
paste oversample_cases.txt oversample_controls.txt > "$PHENO"_oversampled_random_pairs.txt

#downsample to n_pairs
module load  Python/3.8.6-GCCcore-10.2.0
python unrelated_pairs.py "$PHENO" "$n_pairs"

#collect unique included individuals
cut -f1 new_non_sibling_cases/"$PHENO"_random_pairs_stratified_by_age.txt > tmp1.txt
cut -f2 new_non_sibling_cases/"$PHENO"_random_pairs_stratified_by_age.txt > tmp2.txt
cat tmp1.txt tmp2.txt | sort | uniq | awk '{print $1, $1}' > "$PHENO"_included_individuals_age.txt


echo "Final stratified age-matched non-sibling pairs saved to: ${PHENO}_random_pairs_stratified_by_age.txt"

cat "$PHENO"_included_individuals_age.txt individuals_in_pairs_with_"$PHENO".txt ukbb_pc1_16.txt.empty > test_"$PHENO"_sibs_nonsibs.ids.no_pc
