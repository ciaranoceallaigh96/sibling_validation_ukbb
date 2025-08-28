#FINDS SIGNIFICANT PCS RELATIVE TO PHENOTYPE via REGRESSION
#Regress each PC individually on each phenotype 

library(dplyr)

pheno_vec <- c("CAD_UKBB", "T2D_UKBB",
               "prostate_cancer_self_report",
               "breast_cancer_self_report",
               "ea_quant_binary")

pc_file <- "ukbb_pc1_16.txt"
pcs <- read.table(
  pc_file,
  header = TRUE,
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

pval_matrix <- matrix(NA, nrow = 16, ncol = length(pheno_vec),
                      dimnames = list(paste0("PC", 1:16), pheno_vec))

for (pheno in pheno_vec) {

  cat("Phenotype:", pheno, "\n")

  pheno_path <- paste0(pheno, ".full.pheno")
  pheno_data <- read.table(pheno_path, header = FALSE)
  colnames(pheno_data) <- c("FID", "IID", "PHENO")

  merged <- merge(pheno_data, pcs, by = c("FID", "IID"))
  merged$PHENO <- as.numeric(as.character(merged$PHENO))

  merged <- merged %>% filter(PHENO %in% c(1, 2))
  merged$PHENO <- ifelse(merged$PHENO == 2, 1, 0)

  model <- glm(PHENO ~ ., data = merged[, c("PHENO", paste0("PC", 1:16))], family = binomial())
  coef_df <- as.data.frame(summary(model)$coefficients)

  pc_pvals <- coef_df[grep("PC", rownames(coef_df)), "Pr(>|z|)"]
  names(pc_pvals) <- paste0("PC", 1:16)

  pval_matrix[names(pc_pvals), pheno] <- pc_pvals
}

#LaTeX table 
cat("\n\\begin{tabular}{l", paste(rep("c", length(pheno_vec)), collapse = ""), "}\n", sep = "")
cat("PC & ", paste(pheno_vec, collapse = " & "), " \\\\\n")
cat("\\hline\n")

for (pc in rownames(pval_matrix)) {
  row_vals <- sapply(pval_matrix[pc, ], function(p) {
    if (is.na(p)) return("NA")
    formatted <- formatC(p, format = "e", digits = 1)
    if (p < 0.05) paste0("\\textbf{", formatted, "}") else formatted
  })
  cat(pc, " & ", paste(row_vals, collapse = " & "), " \\\\\n", sep = "")
}

cat("\\end{tabular}\n")
