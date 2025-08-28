#mixed-effects logistic model to predict whether the case in a pair has the higher PRS (correct = 1). PC-selection method, pairType, SNP-set-size and phenotype are treated as fixed effects. Model estimates a main effect for each factor plus interaction between method and pair type.
#(i) which PRS-construction method is best overall
#(ii) how much accuracy drops when you move from unrelated population pairs to sibling pairs
#(iii) whether attenuation differs by method



pkgs <- c("data.table", "dplyr", "tidyr", "stringr",
          "purrr", "lme4", "emmeans", "car", "glue")
invisible(lapply(pkgs, require, character.only = TRUE))

prs_dir        <- "prs_results"
sib_pair_dir   <- "new_sibling_cases"
pop_pair_dir   <- "new_non_sibling_cases"
pc_file        <- "../population_data/ukbb_pc1_16.txt"

pheno_vec      <- c("CAD_UKBB", "T2D_UKBB",
                    "prostate_cancer_self_report",
                    "breast_cancer_self_report",
                    "ea_quant_binary")

method_vec     <- c("w_0", "w_16", "lmm_16")
pairType_vec   <- c("sibs", "non_sibs")
snpset_vec     <- c("1000", "10000", "100000")
covars         <- paste0("PC", 1:16)

#Read and residualise one PRS file
read_prs <- function(pairType, pheno, method, snpset_size) {

  file_path <- glue("{prs_dir}/{pairType}_prs_{method}_pc_{pheno}_{snpset_size}.sscore")
  if (!file.exists(file_path))
    stop(glue("File not found: {file_path}"))

  prs <- fread(file_path, fill = TRUE) %>%
    select(FID, IID, PHENO1, SCORE1_AVG) %>%
    transmute(
      FID,
      IID,
      SCORE   = SCORE1_AVG,
      is_case = as.integer(PHENO1) - 1L
    )

  if (method %in% c("w_16", "lmm_16")) {
    covar <- fread(pc_file, fill = TRUE, header = TRUE) %>%
      rename_with(~ sub("^#", "", .x)) %>%
      select(FID, all_of(covars))
    prs   <- left_join(prs, covar, by = "FID")
    good  <- prs %>% select(all_of(covars)) %>% complete.cases()
    if (sum(good) > 10) {
      df_lm <- data.frame(
        SCORE = prs$SCORE[good],
        prs[good, covars, with = FALSE]
      )
      form   <- as.formula(paste("SCORE ~", paste(covars, collapse = " + ")))
      lm_mod <- lm(form, data = df_lm)
      prs$SCORE[good] <- residuals(lm_mod)
      message(glue("{pairType}-{pheno}-{method}-{snpset_size}: residualised {sum(good)} rows"))
    } else {
      warning(glue("{pairType}-{pheno}-{method}-{snpset_size}: <10 complete PC rows – SCORE left raw"))
    }
  }

  prs
}

##Score one set of pairs
score_pairs <- function(pairType, pheno, method, snpset_size) {

  prs_tbl <- read_prs(pairType, pheno, method, snpset_size)

  prs_id1 <- prs_tbl %>% rename(SCORE1 = SCORE, ISCASE1 = is_case)
  prs_id2 <- prs_tbl %>% rename(SCORE2 = SCORE, ISCASE2 = is_case)

  pair_file <- if (pairType == "sibs") {
    glue("{sib_pair_dir}/unique_filtered_sibling_pairs_with_{pheno}.csv")
  } else {
    glue("{pop_pair_dir}/{pheno}_random_pairs_stratified_by_age.txt")
  }
  pairs <- fread(pair_file, header = (pairType == "sibs"))[, 1:2]
  setnames(pairs, c("ID1", "ID2"))

  df <- pairs %>%
    left_join(prs_id1, by = c("ID1" = "IID")) %>%
    left_join(prs_id2, by = c("ID2" = "IID")) %>%
    filter(!is.na(SCORE1), !is.na(SCORE2), ISCASE1 != ISCASE2) %>%
    mutate(
      case_score    = ifelse(ISCASE1 == 1, SCORE1, SCORE2),
      control_score = ifelse(ISCASE1 == 1, SCORE2, SCORE1),
      correct       = case_when(
                        case_score > control_score ~ 1,
                        case_score < control_score ~ 0,
                        TRUE                        ~ NA_real_
                      )
    ) %>%
    filter(!is.na(correct)) %>%
    mutate(
      pairType     = pairType,
      method       = method,
      pheno        = pheno,
      snpset_size  = snpset_size,
      PairID       = paste(ID1, ID2, sep = "_")
    ) %>%
    select(PairID, pairType, method, pheno, snpset_size, correct)

  df
}

message("Scoring all pairs across SNP-set sizes ...")
all_pairs <- pmap_dfr(
  expand.grid(
    pairType     = pairType_vec,
    pheno        = pheno_vec,
    method       = method_vec,
    snpset_size  = snpset_vec,
    stringsAsFactors = FALSE
  ),
  score_pairs
)

all_pairs <- all_pairs %>%
  mutate(
    method      = factor(method,      levels = method_vec),
    pairType    = factor(pairType,    levels = pairType_vec),
    pheno       = factor(pheno,       levels = pheno_vec),
    snpset_size = factor(snpset_size, levels = snpset_vec)
  )

#Fit GLMM
message("Fitting GLMM ...")
glmm_fit <- glmer(
  correct ~ method * pairType + snpset_size + pheno + (1 | PairID),
  data    = all_pairs,
  family  = binomial(link = "logit"),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl   = list(maxfun = 1e6)
  )
)

cat("\nGLMM summary\n");  print(summary(glmm_fit))

cat("\nType-III χ² tests\n")
print(car::Anova(glmm_fit, type = "III"))

cat("\nMarginal P(correct) by method & pairType\n")
emm <- emmeans(glmm_fit, ~ method | pairType, type = "response")
print(emm)

cat("\nAttenuation (sib – pop) within each method\n")
att <- contrast(emm, method = "revpairwise", by = "method")
print(summary(att, infer = c(TRUE, TRUE)))

