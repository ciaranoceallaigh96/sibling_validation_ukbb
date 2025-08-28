#SCRIPT TO EXTRACT SIBLING PAIRS FROM UKBB WHERE ONE IS A CASE AND ONE IS A CONTROL (WILL REMOVE TRIOS AND HIGHER AND ONLY KEEP ONE PAIR)

import pandas as pd
import sys

phenotype = str(sys.argv[1])

cases_file = f"{phenotype}.ids"
cases = set(pd.read_csv(cases_file, header=None)[0])

siblings_file = "ukb_rel_sibs.csv"
siblings = pd.read_csv(siblings_file)

fam_file = "all_maf_geno_mind.fam"
fam_data = pd.read_csv(fam_file, delim_whitespace=True, header=None, names=["FID", "IID", "PID", "MID", "Sex", "Pheno"])
present_iids = set(fam_data["IID"])

#Filter to keep only pairs where both siblings are present in the fam file
siblings = siblings[(siblings['ID1'].isin(present_iids)) & (siblings['ID2'].isin(present_iids))]

fam_data_dict = dict(zip(fam_data["IID"], fam_data["Sex"]))
siblings['Phenotype_ID1'] = siblings['ID1'].apply(lambda x: 1 if x in cases else 0)
siblings['Phenotype_ID2'] = siblings['ID2'].apply(lambda x: 1 if x in cases else 0)

siblings['Sex_ID1'] = siblings['ID1'].map(fam_data_dict)
siblings['Sex_ID2'] = siblings['ID2'].map(fam_data_dict)

#Filter pairs where one sibling is a case and the other is not
filtered_siblings = siblings[
    ((siblings['Phenotype_ID1'] == 1) & (siblings['Phenotype_ID2'] == 0)) |
    ((siblings['Phenotype_ID1'] == 0) & (siblings['Phenotype_ID2'] == 1))
]

#filtering for sex-specific phenotypes
if phenotype in ['prostate_cancer_self_report', 'breast_cancer_self_report']:
    filtered_siblings = filtered_siblings[filtered_siblings['Sex_ID1'] == filtered_siblings['Sex_ID2']]


#Identify all individuals involved in the filtered sibling pairs
individual_pairs = pd.concat([filtered_siblings['ID1'], filtered_siblings['ID2']])
grouped_ids = individual_pairs.value_counts()

#Process each group of overlapping IDs
unique_pairs = []
used_ids = set()

for _, row in filtered_siblings.iterrows():
    id1, id2 = row['ID1'], row['ID2']
    if id1 not in used_ids and id2 not in used_ids:
        unique_pairs.append(row)
        #Mark these IDs as used
        used_ids.add(id1)
        used_ids.add(id2)

unique_filtered_siblings = pd.DataFrame(unique_pairs)

columns_to_fix = ['ID1', 'ID2', 'Phenotype_ID1', 'Phenotype_ID2']
unique_filtered_siblings[columns_to_fix] = unique_filtered_siblings[columns_to_fix].applymap(
    lambda x: str(x).replace('.0', '') if isinstance(x, float) else x
)

output_file = f"unique_filtered_sibling_pairs_with_{phenotype}.csv"
unique_filtered_siblings.to_csv(output_file, index=False)

individuals = pd.concat([
    unique_filtered_siblings[['ID1']].rename(columns={'ID1': 'IID'}),
    unique_filtered_siblings[['ID2']].rename(columns={'ID2': 'IID'})
])


individuals['IID'] = individuals['IID'].astype(str)
fam_data['IID'] = fam_data['IID'].astype(str)

individuals = individuals.merge(fam_data[['FID', 'IID']], on='IID', how='left')
individuals_output_file = f"individuals_in_pairs_with_{phenotype}.txt"
individuals.to_csv(individuals_output_file, sep=' ', index=False, header=False)

print(f"Unique sibling pairs with phenotype information saved to {output_file}")
print(f"All individuals in pairs saved to {individuals_output_file}")
