library(ggplot2)
library(pROC)

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
pheno <- args[2]
sibs <- args[3] #True or False
sibs <- as.logical(args[3])
pc_number <- as.integer(args[4]) #number of PCs to include

log_file <- paste0(prefix, ".resid.log")
cat("Log printed to: ", log_file, "\n\n")

log_con <- file(log_file, open = "at")
sink(log_con)
sink(log_con, type = "message")


cat("Log started at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Log printed to: ", log_file, "\n\n")


#Function to calculate Nagelkerke's R2 (gives identical results to rcompanion nagelkerke function)

nagelkerke_r2 <- function(model) {
     null_model <- glm(PHENO ~ 1, data = merged_df, family = binomial) # Null model
     llf_model <- logLik(model) # Log-likelihood of the fitted model
     llf_null <- logLik(null_model) # Log-likelihood of the null model
     n <- length(merged_df$PHENO) # Number of observations
     r2 <- 1 - exp((2 / n) * (llf_null - llf_model)) # Cox & Snell's R2
     nagelkerke <- r2 / (1 - exp((2 / n) * llf_null)) # Nagelkerke's R2
     return(nagelkerke)
}


profile_file <- paste0(prefix, ".sscore") #PRS file
cat("profile_file: ", profile_file, "\n")



covar_file <- ("ukbb_pc1_16.txt")
cat("covar_file: ", covar_file, "\n\n")
covars <- paste0("PC", 1:pc_number)
cat("no. pcs for residualization: ", pc_number, "\n\n")



prs_data <- read.table(profile_file, header = TRUE)


if (length(unique(prs_data$PHENO)) == 2) {
  cat("Phenotype is BINARY therefore Nagelkerke's R2 and AUC will be calculated. \n")
  family_type <- "binomial"
  prs_data$PHENO <- as.numeric(prs_data$PHENO1) - 1
} else {
  cat("Phenotype is CONTINUOUS - Pearson correlation will be calculated .\n ")
  names(prs_data)[names(prs_data) == "PHENO1"] <- "PHENO"
  family_type <- "gaussian"
}



prs_data$SCORE <- prs_data$SCORE1_AVG


cov_data <- read.table(covar_file, header = TRUE, fill = TRUE, blank.lines.skip = TRUE)

merged_df <- merge(prs_data, cov_data, by = "FID", all.x = TRUE)



#Check
if (anyNA(merged_df[, c("SCORE", covars)])) {
  stop("Missing values detected in SCORE or covariates.")
}

#RESIDULAIZE THE PRS BY REGRESSING OUT THE COVARIATES
lm_formula <- as.formula(paste("SCORE ~", paste(covars, collapse = " + ")))
lm_model <- lm(lm_formula, data = merged_df)
merged_df$SCORE <- residuals(lm_model)
cat("\nTHE PRS HAS NOW BEEN RESIDUALIZED (REGRESSION WITH COVARIATES)\n")


if (sibs == TRUE) {
        sibling_pairs <- read.csv(paste0("unique_filtered_sibling_pairs_with_",  pheno , ".csv"))
        all_ids <- unique(c(sibling_pairs$ID1, sibling_pairs$ID2))
        missing_ids <- setdiff(all_ids, merged_df$IID.x)
        if (length(missing_ids) > 0) {
                stop(paste("CK Error: The following IDs are missing in merged_df:", paste(missing_ids, collapse = ", ")))
        }
        pair_data <- merge(sibling_pairs, merged_df[, c("IID.x", "PHENO", "SCORE", "FID")],by.x = "ID1", by.y = "IID.x", all.x = TRUE)
        names(pair_data)[names(pair_data) %in% c("PHENO", "SCORE", "FID")] <- c("PHENO1", "SCORE1", "FID1")
        pair_data <- merge(pair_data, merged_df[, c("IID.x", "PHENO", "SCORE", "FID")],by.x = "ID2", by.y = "IID.x", all.x = TRUE)
        names(pair_data)[names(pair_data) %in% c("PHENO", "SCORE", "FID")] <- c("PHENO2", "SCORE2", "FID2")
        invalid_rows <- which((pair_data$PHENO1 + pair_data$PHENO2) != 1)
        if (length(invalid_rows) > 0) {
                stop("CK Error: Some sibling pairs do not consist of one case and one control.")
        }
        tie_rows <- which(pair_data$SCORE1 == pair_data$SCORE2)
        if (length(tie_rows) > 0) {
                pair_data <- pair_data[-tie_rows, ]
        }
        correct_classifications <- sum(
        (pair_data$PHENO1 == 1 & pair_data$SCORE1 > pair_data$SCORE2) |
        (pair_data$PHENO2 == 1 & pair_data$SCORE2 > pair_data$SCORE1)
        )
        total_pairs <- nrow(pair_data)
        accuracy <- correct_classifications / total_pairs * 100
        cat(sprintf("Resid-PS correctly classified the case in %.2f%% of sibling pairs (excluding ties).\n", accuracy))
}

if (sibs == FALSE) {
        sibling_pairs <- read.table(paste0(pheno, "_random_pairs_stratified_by_age.txt"),header = FALSE,col.names = c("ID1", "ID2"),sep = "\t")
        all_ids <- unique(c(sibling_pairs$ID1, sibling_pairs$ID2))
        missing_ids <- setdiff(all_ids, merged_df$IID.x)
        if (length(missing_ids) > 0) {
                stop(paste("CK Error: The following IDs are missing in merged_df:", paste(missing_ids, collapse = ", ")))
        }
        pair_data <- merge(sibling_pairs, merged_df[, c("IID.x", "PHENO", "SCORE", "FID")],by.x = "ID1", by.y = "IID.x", all.x = TRUE)
        names(pair_data)[names(pair_data) %in% c("PHENO", "SCORE", "FID")] <- c("PHENO1", "SCORE1", "FID1")
        pair_data <- merge(pair_data, merged_df[, c("IID.x", "PHENO", "SCORE", "FID")],by.x = "ID2", by.y = "IID.x", all.x = TRUE)
        names(pair_data)[names(pair_data) %in% c("PHENO", "SCORE", "FID")] <- c("PHENO2", "SCORE2", "FID2")
        invalid_rows <- which((pair_data$PHENO1 + pair_data$PHENO2) != 1)
        if (length(invalid_rows) > 0) {
                stop("CK Error: Some NON-sibling pairs do not consist of one case and one control.")
        }
        tie_rows <- which(pair_data$SCORE1 == pair_data$SCORE2)
        if (length(tie_rows) > 0) {
                pair_data <- pair_data[-tie_rows, ]
        }
        correct_classifications <- sum(
        (pair_data$PHENO1 == 1 & pair_data$SCORE1 > pair_data$SCORE2) |
        (pair_data$PHENO2 == 1 & pair_data$SCORE2 > pair_data$SCORE1)
        )
        total_pairs <- nrow(pair_data)
        accuracy <- correct_classifications / total_pairs * 100
        cat(sprintf("Resid-PS correctly classified the case in %.2f%% of NON-sibling pairs (excluding ties).\n", accuracy))
}
