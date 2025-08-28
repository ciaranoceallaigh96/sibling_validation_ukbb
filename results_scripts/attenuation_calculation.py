#CALCULATION IS ON RESIDUALS (I.E. PCS INCLUDED AS COVARIATES IN PRS)

import re
import pandas as pd

phenotypes = ["CAD_UKBB", "T2D_UKBB", "breast_cancer_self_report", "prostate_cancer_self_report", "ea_quant_binary"]
snp_sizes = ["1000", "10000", "100000"]

methods = {
    "lmm": ("sibs_prs_lmm_16_pc", "non_sibs_prs_lmm_16_pc"),
    "glm": ("sibs_prs_w_16_pc", "non_sibs_prs_w_16_pc")
}

#extract classification %
pattern = re.compile(r"correctly classified the case in (\d+\.\d+)%")

data = {method: {size: {} for size in snp_sizes} for method in methods}

for method, (sibs_prefix, nonsibs_prefix) in methods.items():
    for size in snp_sizes:
        for pheno in phenotypes:
            try:
                sib_file = f"{sibs_prefix}_{pheno}_{size}.resid.log"
                with open(sib_file) as f:
                    sib_pct = float(pattern.search(f.read()).group(1))
            except:
                sib_pct = None

            try:
                nonsib_file = f"{nonsibs_prefix}_{pheno}_{size}.resid.log"
                with open(nonsib_file) as f:
                    nonsib_pct = float(pattern.search(f.read()).group(1))
            except:
                nonsib_pct = None

            data[method][size][pheno] = {"sibs": sib_pct, "nonsibs": nonsib_pct}


rows = []

for pheno in phenotypes:
    row = {"Phenotype": pheno}
    for method in methods:
        for size in snp_sizes:
            vals = data[method][size][pheno]
            sibs, nonsibs = vals["sibs"], vals["nonsibs"]

            #Attenuation calculation
            if sibs is not None and nonsibs not in (None, 50.0):
                d_sibs = sibs - 50.0
                d_nonsibs = nonsibs - 50.0
                attenuation = (1 - (d_sibs / d_nonsibs)) * 100
                row[f"{size}_{method}"] = round(attenuation, 2)
            else:
                row[f"{size}_{method}"] = None

    rows.append(row)

df_att = pd.DataFrame(rows).set_index("Phenotype")
print(df_att)


rows_sibs = []
for pheno in phenotypes:
    row = {"Phenotype": pheno}
    for method in methods:
        for size in snp_sizes:
            row[f"{size}_{method}"] = data[method][size][pheno]["sibs"]
    rows_sibs.append(row)
df_sibs = pd.DataFrame(rows_sibs).set_index("Phenotype")

#Non-sibs
rows_nonsibs = []
for pheno in phenotypes:
    row = {"Phenotype": pheno}
    for method in methods:
        for size in snp_sizes:
            row[f"{size}_{method}"] = data[method][size][pheno]["nonsibs"]
    rows_nonsibs.append(row)
df_nonsibs = pd.DataFrame(rows_nonsibs).set_index("Phenotype")



print("Raw classification % — Siblings")
print(df_sibs)

print("\nRaw classification % — Non-siblings")
print(df_nonsibs)

