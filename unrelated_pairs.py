#Created non-sibling pairs based on sibling pair characteristics (requires unrelated_pairs.sh to be run)
import pandas as pd
import numpy as np
import sys

pheno = sys.argv[1]
n_target = int(sys.argv[2])

age_file = 'age_info_participant.csv'
sibling_file = f'unique_filtered_sibling_pairs_with_{pheno}.csv'
random_pair_file = f'{pheno}_oversampled_random_pairs.txt'

age_df = pd.read_csv(age_file)
age_df.columns = ['eid', 'age']

#Sibling Pairs 
sib_df = pd.read_csv(sibling_file)
sib_df = sib_df.merge(age_df, left_on='ID1', right_on='eid', how='left') \
               .merge(age_df, left_on='ID2', right_on='eid', how='left', suffixes=('_1', '_2'))
sib_df = sib_df.dropna(subset=['age_1', 'age_2'])
sib_df['age_diff'] = (sib_df['age_1'] - sib_df['age_2']).abs()

#Bin sibling age differences
bin_width = 2
max_diff = sib_df['age_diff'].max()
bins = np.arange(0, max_diff + bin_width, bin_width)
sib_df['bin'] = pd.cut(sib_df['age_diff'], bins)
bin_props = sib_df['bin'].value_counts(normalize=True).sort_index()

#Random Non-Sibling Pairs 
rand_df = pd.read_csv(random_pair_file, sep='\t', header=None, names=['ID1', 'ID2'])
rand_df = rand_df.merge(age_df, left_on='ID1', right_on='eid', how='left') \
                 .merge(age_df, left_on='ID2', right_on='eid', how='left', suffixes=('_1', '_2'))
rand_df = rand_df.dropna(subset=['age_1', 'age_2'])
rand_df['age_diff'] = (rand_df['age_1'] - rand_df['age_2']).abs()
rand_df['bin'] = pd.cut(rand_df['age_diff'], bins)

#Stratified sampling
stratified_sample = []
for bin_label, prop in bin_props.items():
    bin_subset = rand_df[rand_df['bin'] == bin_label]
    if bin_subset.empty:
        continue
    n_samples = int(round(n_target * prop))
    if len(bin_subset) < n_samples:
        n_samples = len(bin_subset)
    stratified_sample.append(bin_subset.sample(n=n_samples, random_state=42))

#Combine and randomly trim to exactly n_target 
final_df = pd.concat(stratified_sample)

if len(final_df) < n_target:
    deficit = n_target - len(final_df)

    # Find the bin with the most data available
    largest_bin_label = bin_props.idxmax()
    top_up_pool = rand_df[rand_df['bin'] == largest_bin_label]

    # Exclude already selected pairs
    already_selected = final_df[['ID1', 'ID2']].apply(tuple, axis=1).tolist()
    top_up_pool = top_up_pool[~top_up_pool[['ID1', 'ID2']].apply(tuple, axis=1).isin(already_selected)]

    top_up_sample = top_up_pool.sample(n=min(deficit, len(top_up_pool)), random_state=42)
    final_df = pd.concat([final_df, top_up_sample])

#If overshot, trim back down
if len(final_df) > n_target:
    final_df = final_df.sample(n=n_target, random_state=42)





output_path = f'{pheno}_random_pairs_stratified_by_age.txt'
final_df[['ID1', 'ID2']].to_csv(output_path, sep='\t', index=False, header=False)

print(f"Stratified sample saved with {len(final_df)} pairs: {output_path}")

