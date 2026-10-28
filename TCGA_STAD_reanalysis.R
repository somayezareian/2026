############################################################
# TCGA-STAD clinical validation of the frozen 15-gene signature
# Revised, auditable analysis pipeline
#
# Inputs (paths retained from the user's previous script):
#   TCGA-STAD.star_tpm.tsv.gz
#   TCGA-STAD.survival.tsv.gz
#   TCGA-STAD.clinical.tsv.gz
#
# Primary endpoint: overall survival (OS)
# Primary predictor: continuous signature score, expressed per 1 SD
# Main adjusted model: signature + age + sex + early/advanced stage
# Secondary: KM/log-rank, gene-wise Cox (BH), stage associations,
#            exploratory advanced-stage subgroup analysis.
############################################################

options(stringsAsFactors = FALSE, scipen = 999)

# 1. Packages ------------------------------------------------
required_packages <- c("data.table", "dplyr", "ggplot2", "survival",
                       "AnnotationDbi", "org.Hs.eg.db")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install these packages before running: ",
       paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# 2. File paths ----------------------------------------------
expression_file <- "C:/Users/somay/Documents/TCGA/TCGA-STAD.star_tpm.tsv.gz"
survival_file   <- "C:/Users/somay/Documents/TCGA/TCGA-STAD.survival.tsv.gz"
clinical_file   <- "C:/Users/somay/Documents/TCGA/TCGA-STAD.clinical.tsv.gz"

input_files <- c(expression = expression_file,
                 survival = survival_file,
                 clinical = clinical_file)
print(data.frame(input = names(input_files), path = unname(input_files),
                 exists = file.exists(input_files)))
if (!all(file.exists(input_files))) {
  stop("One or more input files were not found. Correct the paths in section 2.")
}
output_dir <- file.path(dirname(expression_file), "TCGA_STAD_reanalysis_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 3. Frozen signature ----------------------------------------
signature_genes <- c(
  "VIM", "ITGB1", "LGALS1", "S100A4", "PRNP", "IFITM2", "CD44",
  "LAMC2", "CNN3", "PGK1", "CLIC4", "CSRP1", "SYPL1", "CD99", "CLU"
)
minimum_genes_required <- 12L

# 4. Helper functions ----------------------------------------
normalise_name <- function(x) tolower(gsub("[^a-z0-9]", "", x))
find_col <- function(dat, candidates, required = TRUE, label = "column") {
  nms <- names(dat)
  hits <- match(normalise_name(candidates), normalise_name(nms), nomatch = 0L)
  hits <- hits[hits > 0L]
  if (!length(hits)) {
    if (required) stop("Could not identify ", label, ". Available columns: ",
                       paste(nms, collapse = ", "))
    return(NA_character_)
  }
  nms[hits[1]]
}
patient_id_from <- function(x) {
  x <- trimws(as.character(x))
  # TCGA patient barcode is the first 12 characters (TCGA-XX-XXXX).
  ifelse(grepl("^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}", x), substr(x, 1, 12), NA_character_)
}
first_nonmissing <- function(x) {
  x <- x[!is.na(x) & trimws(as.character(x)) != "" &
           !tolower(trimws(as.character(x))) %in% c("na", "n/a", "not reported", "unknown")]
  if (!length(x)) NA else x[1]
}
parse_event <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(z))
  out[z %in% c("dead", "deceased", "death", "1", "true", "yes", "event")] <- 1L
  out[z %in% c("alive", "living", "0", "false", "no", "censored", "censor")] <- 0L
  # Numeric encodings are accepted only if exactly 0/1; anything else stays missing.
  num <- suppressWarnings(as.numeric(z))
  ok <- !is.na(num) & num %in% c(0, 1)
  out[ok] <- as.integer(num[ok])
  out
}
parse_stage <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z <- gsub("^STAGE\\s*", "", z)
  z <- gsub("PATHOLOGIC", "", z)
  z <- gsub("AJCC", "", z)
  z <- gsub("STAGE", "", z)
  z <- gsub("[^IV0-9]", "", z)
  out <- rep(NA_character_, length(z))
  out[grepl("^I$|^IA$|^IB$|^IC$|^II$|^IIA$|^IIB$|^IIC$", z)] <- "Early"
  out[grepl("^III$|^IIIA$|^IIIB$|^IIIC$|^IV$|^IVA$|^IVB$", z)] <- "Advanced"
  out
}
extract_cox <- function(fit, model_name) {
  s <- summary(fit)
  cf <- as.data.frame(s$coefficients)
  ci <- as.data.frame(s$conf.int)
  data.frame(
    model = model_name,
    term = rownames(cf),
    HR = ci$`exp(coef)`,
    lower95 = ci$`lower .95`,
    upper95 = ci$`upper .95`,
    p_value = cf$`Pr(>|z|)`,
    n = fit$n,
    events = fit$nevent,
    censored = fit$n - fit$nevent,
    row.names = NULL
  )
}
write_csv <- function(x, filename) {
  data.table::fwrite(x, file.path(output_dir, filename), na = "")
}

# 5. Read and prepare STAR-TPM expression ---------------------
message("Reading STAR-TPM expression...")
expr_raw <- fread(expression_file, data.table = FALSE, check.names = FALSE)
ens_col <- find_col(expr_raw,
                    c("Ensembl_ID", "ensembl_id", "ENSEMBL", "Ensembl", "gene_id", "Gene_ID"),
                    label = "Ensembl gene-ID column")
expr_cols <- setdiff(names(expr_raw), ens_col)
if (length(expr_cols) < 10L) stop("Too few expression columns were detected.")
ens_ids <- sub("\\..*$", "", as.character(expr_raw[[ens_col]]))
expr_num <- as.matrix(expr_raw[, expr_cols, drop = FALSE])
suppressWarnings(storage.mode(expr_num) <- "numeric")
rownames(expr_num) <- ens_ids
if (any(!is.finite(expr_num) & !is.na(expr_num))) stop("Non-finite values found in expression matrix.")

# Average duplicate Ensembl rows on the TPM scale before log transformation.
if (anyDuplicated(rownames(expr_num))) {
  expr_dt <- data.frame(ENSEMBL = rownames(expr_num), expr_num, check.names = FALSE)
  expr_dt <- expr_dt %>% group_by(ENSEMBL) %>%
    summarise(across(where(is.numeric), ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)),
              .groups = "drop")
  expr_num <- as.matrix(expr_dt[, setdiff(names(expr_dt), "ENSEMBL"), drop = FALSE])
  rownames(expr_num) <- expr_dt$ENSEMBL
}
if (any(expr_num < 0, na.rm = TRUE)) stop("Negative values found; expected STAR-derived TPM values.")

# Restrict to primary tumour samples (barcode sample-type code 01) when TCGA
# sample barcodes are present. If columns are patient-level IDs, keep them.
sample_ids_all <- colnames(expr_num)
sample_type <- ifelse(grepl("^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-", sample_ids_all),
                      substr(sample_ids_all, 14, 15), NA_character_)
if (any(sample_type == "01", na.rm = TRUE)) {
  keep_sample <- is.na(sample_type) | sample_type == "01"
  expr_num <- expr_num[, keep_sample, drop = FALSE]
  message("Retained primary tumour (01) samples and any patient-level IDs; removed other TCGA sample types.")
} else {
  message("No recognizable TCGA primary-tumour sample codes detected; no sample-type filter applied.")
}

# Use log2(TPM + 1) as the expression transformation.
expr_log <- log2(expr_num + 1)

# Map Ensembl IDs to official symbols. Any multiple Ensembl rows mapping to a
# signature symbol are averaged after log2 transformation.
ens_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(rownames(expr_log)),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)
map_df <- data.frame(ENSEMBL = names(ens_map), SYMBOL = unname(ens_map)) %>%
  filter(!is.na(SYMBOL), SYMBOL != "", SYMBOL %in% signature_genes)
write_csv(map_df, "signature_ensembl_symbol_mapping.csv")
if (!nrow(map_df)) stop("No signature genes mapped. Check Ensembl IDs and annotation package version.")

matched_rows <- intersect(rownames(expr_log), map_df$ENSEMBL)
map_df <- map_df[match(matched_rows, map_df$ENSEMBL), , drop = FALSE]
expr_sig_rows <- expr_log[map_df$ENSEMBL, , drop = FALSE]
rownames(expr_sig_rows) <- map_df$SYMBOL
# Aggregate multiple Ensembl IDs per symbol.
expr_by_symbol <- do.call(rbind, lapply(signature_genes[signature_genes %in% rownames(expr_sig_rows)], function(g) {
  mat <- expr_sig_rows[rownames(expr_sig_rows) == g, , drop = FALSE]
  if (nrow(mat) == 1L) as.numeric(mat[1, ]) else colMeans(mat, na.rm = TRUE)
}))
if (is.null(dim(expr_by_symbol))) expr_by_symbol <- matrix(expr_by_symbol, nrow = 1)
rownames(expr_by_symbol) <- signature_genes[signature_genes %in% rownames(expr_sig_rows)]
colnames(expr_by_symbol) <- colnames(expr_num)

present_genes <- rownames(expr_by_symbol)
missing_genes <- setdiff(signature_genes, present_genes)
write_csv(data.frame(gene = signature_genes, present = signature_genes %in% present_genes),
          "signature_gene_coverage.csv")
message("Signature genes mapped: ", length(present_genes), "/", length(signature_genes))
if (length(present_genes) < minimum_genes_required) {
  stop("Only ", length(present_genes), " signature genes mapped; at least ",
       minimum_genes_required, " are required by the prespecified rule.")
}

# Convert expression sample columns to patient IDs. For multiple primary tumour
# samples/aliquots per patient, average the log2 expression values per patient.
expr_patient_ids <- patient_id_from(colnames(expr_by_symbol))
valid_expr <- !is.na(expr_patient_ids)
expr_by_symbol <- expr_by_symbol[, valid_expr, drop = FALSE]
expr_patient_ids <- expr_patient_ids[valid_expr]
if (!length(expr_patient_ids)) stop("Could not derive TCGA patient IDs from expression column names.")
expr_patient <- do.call(cbind, lapply(unique(expr_patient_ids), function(pid) {
  idx <- which(expr_patient_ids == pid)
  if (length(idx) == 1L) as.numeric(expr_by_symbol[, idx]) else rowMeans(expr_by_symbol[, idx, drop = FALSE], na.rm = TRUE)
}))
rownames(expr_patient) <- rownames(expr_by_symbol)
colnames(expr_patient) <- unique(expr_patient_ids)
if (anyDuplicated(colnames(expr_patient))) stop("Duplicate patient IDs remain in expression matrix.")

# 6. Read survival and define OS endpoint ---------------------
message("Reading survival data...")
surv_raw <- fread(survival_file, data.table = FALSE, check.names = FALSE)
surv_id_col <- find_col(surv_raw,
  c("sample", "sample_id", "patient", "patient_id", "submitter_id", "barcode", "case_id", "bcr_patient_barcode"),
  label = "survival patient-ID column")
surv_time_col <- find_col(surv_raw,
  c("OS.time", "OS_time", "os_time", "OS.time.days", "overall_survival_time", "time"),
  label = "overall-survival time column")
surv_event_col <- find_col(surv_raw,
  c("OS", "OS.event", "OS_event", "vital_status", "event", "death"),
  label = "overall-survival event column")

surv_df <- data.frame(
  patient_id = patient_id_from(surv_raw[[surv_id_col]]),
  OS_time_days = suppressWarnings(as.numeric(as.character(surv_raw[[surv_time_col]]))),
  OS_event_raw = as.character(surv_raw[[surv_event_col]])
) %>% mutate(OS_event = parse_event(OS_event_raw))

# Explicitly validate the event coding. Unrecognized event labels are excluded
# and counted, rather than being silently interpreted as censored.
invalid_survival <- is.na(surv_df$patient_id) | is.na(surv_df$OS_time_days) |
  !is.finite(surv_df$OS_time_days) | surv_df$OS_time_days <= 0 | is.na(surv_df$OS_event)
message("Survival rows read: ", nrow(surv_df))
message("Rows excluded for missing/invalid patient ID, time, or event: ", sum(invalid_survival))
surv_valid <- surv_df[!invalid_survival, , drop = FALSE]
if (any(!surv_valid$OS_event %in% c(0L, 1L))) stop("Event variable contains values other than 0/1 after parsing.")

# If the source has duplicate patient records, only collapse when time and event
# agree; conflicting survival records stop the script for manual review.
dup_ids <- unique(surv_valid$patient_id[duplicated(surv_valid$patient_id)])
if (length(dup_ids)) {
  conflicts <- vapply(dup_ids, function(pid) {
    x <- surv_valid[surv_valid$patient_id == pid, , drop = FALSE]
    length(unique(x$OS_time_days)) > 1L || length(unique(x$OS_event)) > 1L
  }, logical(1))
  if (any(conflicts)) stop("Conflicting survival records for patient(s): ",
                           paste(dup_ids[conflicts], collapse = ", "),
                           ". Inspect the survival input; do not choose an arbitrary record.")
  surv_valid <- surv_valid[!duplicated(surv_valid$patient_id), , drop = FALSE]
}

# 7. Read clinical data; derive stage, age and sex ----------------
message("Reading clinical data...")
clin_raw <- fread(clinical_file, data.table = FALSE, check.names = FALSE)
clin_id_col <- find_col(clin_raw,
  c("sample", "sample_id", "patient", "patient_id", "submitter_id", "case_submitter_id", "case_id", "bcr_patient_barcode"),
  label = "clinical patient-ID column")
clin_raw$patient_id <- patient_id_from(clin_raw[[clin_id_col]])
clin_raw <- clin_raw[!is.na(clin_raw$patient_id), , drop = FALSE]

stage_col <- find_col(clin_raw,
  c("ajcc_pathologic_stage.diagnoses", "ajcc_pathologic_stage", "pathologic_stage", "tumor_stage", "stage", "ajcc_stage", "clinical_stage"),
  required = FALSE)
age_col <- find_col(clin_raw,
  c("age_at_initial_pathologic_diagnosis", "age_at_diagnosis.demographic", "age_at_diagnosis", "age_at_diagnosis_years", "age"),
  required = FALSE)
sex_col <- find_col(clin_raw,
  c("gender.demographic", "gender", "sex"), required = FALSE)
message("Clinical columns selected: stage=", ifelse(is.na(stage_col), "NOT FOUND", stage_col),
        "; age=", ifelse(is.na(age_col), "NOT FOUND", age_col),
        "; sex=", ifelse(is.na(sex_col), "NOT FOUND", sex_col))

# Collapse clinical records to one row/patient, taking the first nonmissing
# value per field. This is safe for repeated copies of the same patient record;
# resulting missingness is reported below.
clin_cols <- c(stage_col, age_col, sex_col)
clin_cols <- clin_cols[!is.na(clin_cols)]
if (!length(clin_cols)) stop("No stage, age, or sex columns could be identified in clinical data.")
clin_df <- clin_raw %>% group_by(patient_id) %>%
  summarise(across(all_of(clin_cols), first_nonmissing), .groups = "drop")

if (!is.na(stage_col)) {
  clin_df$stage_raw <- as.character(clin_df[[stage_col]])
} else clin_df$stage_raw <- NA_character_
clin_df$stage_group <- parse_stage(clin_df$stage_raw)

if (!is.na(age_col)) {
  age_num <- suppressWarnings(as.numeric(as.character(clin_df[[age_col]])))
  # Prefer age-at-initial-pathologic-diagnosis (usually years). GDC age_at_diagnosis
  # is often stored in days; infer days only if the column's median is clearly >120.
  if (normalise_name(age_col) == "ageatdiagnosisdemographic" ||
      (sum(!is.na(age_num)) > 0 && median(age_num, na.rm = TRUE) > 120)) {
    age_years <- age_num / 365.25
    message("Age converted from days to years for column: ", age_col)
  } else {
    age_years <- age_num
    message("Age treated as years for column: ", age_col)
  }
  age_years[age_years < 18 | age_years > 110] <- NA_real_
  clin_df$age_years <- age_years
} else clin_df$age_years <- NA_real_

if (!is.na(sex_col)) {
  sex <- tolower(trimws(as.character(clin_df[[sex_col]])))
  sex[sex %in% c("female", "f")] <- "Female"
  sex[sex %in% c("male", "m")] <- "Male"
  sex[!sex %in% c("Female", "Male")] <- NA_character_
  clin_df$sex <- factor(sex, levels = c("Female", "Male"))
} else clin_df$sex <- factor(NA_character_, levels = c("Female", "Male"))

# 8. Join patient-level expression, survival, and clinical data ------------
common_ids <- Reduce(intersect, list(colnames(expr_patient), surv_valid$patient_id))
if (!length(common_ids)) stop("No overlap between expression and valid survival patient IDs.")
expr_patient <- expr_patient[, common_ids, drop = FALSE]
surv_valid <- surv_valid[match(common_ids, surv_valid$patient_id), , drop = FALSE]
stopifnot(identical(colnames(expr_patient), surv_valid$patient_id))

# Gene-wise z-score across the expression/survival-matched patient cohort.
expr_z <- t(scale(t(expr_patient)))
expr_z[!is.finite(expr_z)] <- NA_real_
signature_score <- colMeans(expr_z, na.rm = TRUE)
# Do not allow a patient score based on fewer genes than the declared minimum.
per_patient_gene_count <- colSums(is.finite(expr_z))
signature_score[per_patient_gene_count < minimum_genes_required] <- NA_real_

analysis <- data.frame(
  patient_id = colnames(expr_patient),
  OS_time_days = surv_valid$OS_time_days,
  OS_time_years = surv_valid$OS_time_days / 365.25,
  OS_event = surv_valid$OS_event,
  signature_score = as.numeric(signature_score),
  signature_score_SD = if (is.finite(sd(signature_score, na.rm = TRUE)) && sd(signature_score, na.rm = TRUE) > 0)
    (as.numeric(signature_score) - mean(signature_score, na.rm = TRUE)) / sd(signature_score, na.rm = TRUE) else rep(NA_real_, length(signature_score)),
  genes_in_score = as.integer(per_patient_gene_count),
  stringsAsFactors = FALSE
) %>% left_join(clin_df, by = "patient_id")
analysis$signature_group <- factor(ifelse(analysis$signature_score >= median(analysis$signature_score, na.rm = TRUE),
                                          "High", "Low"), levels = c("Low", "High"))
analysis$stage_group <- factor(analysis$stage_group, levels = c("Early", "Advanced"))
analysis$age_per_10y <- analysis$age_years / 10

# Save full patient-level analysis table and a reproducible cohort-flow table.
write_csv(analysis, "TCGA_STAD_patient_level_analysis_table.csv")
flow <- data.frame(
  step = c("Expression samples after primary-tumour filter",
           "Unique patients with signature expression",
           "Valid unique survival patients",
           "Expression-survival matched patients",
           "Primary continuous-score Cox complete cases",
           "Adjusted Cox complete cases (score + age (per 10 years) + sex + stage)"),
  n = c(ncol(expr_num), ncol(expr_patient), nrow(surv_valid), length(common_ids),
        sum(complete.cases(analysis[, c("OS_time_days", "OS_event", "signature_score")])),
        sum(complete.cases(analysis[, c("OS_time_days", "OS_event", "signature_score_SD", "age_per_10y", "sex", "stage_group")]))),
  stringsAsFactors = FALSE
)
write_csv(flow, "cohort_flow_counts.csv")
write_csv(data.frame(
  metric = c("matched_patients", "matched_deaths", "matched_censored", "missing_signature_score",
             "missing_age", "missing_sex", "missing_stage_group"),
  value = c(nrow(analysis), sum(analysis$OS_event == 1, na.rm = TRUE),
            sum(analysis$OS_event == 0, na.rm = TRUE), sum(is.na(analysis$signature_score)),
            sum(is.na(analysis$age_years)), sum(is.na(analysis$sex)), sum(is.na(analysis$stage_group)))
), "cohort_counts_and_missingness.csv")

# Helper: event and censor counts for each model-specific dataset.
cohort_count <- function(dat, label) data.frame(
  analysis = label, patients = nrow(dat), deaths = sum(dat$OS_event == 1),
  censored = sum(dat$OS_event == 0), stringsAsFactors = FALSE
)
model_cohorts <- list()
cox_tables <- list()
ph_tables <- list()

# 9. Primary univariable continuous-score Cox ------------------
primary_df <- analysis %>% filter(complete.cases(OS_time_days, OS_event, signature_score_SD))
if (nrow(primary_df) < 20 || sum(primary_df$OS_event) < 5) stop("Too few complete cases/events for the primary Cox model.")
model_cohorts[["primary_univariable"]] <- cohort_count(primary_df, "Primary univariable Cox")
fit_primary <- coxph(Surv(OS_time_years, OS_event) ~ signature_score_SD,
                     data = primary_df, x = TRUE, y = TRUE)
cox_tables[["primary_univariable"]] <- extract_cox(fit_primary, "Primary univariable Cox")
ph_primary <- cox.zph(fit_primary)
ph_tables[["primary_univariable"]] <- data.frame(model = "Primary univariable Cox",
  term = rownames(ph_primary$table), chisq = ph_primary$table[, "chisq"],
  df = ph_primary$table[, "df"], p_value = ph_primary$table[, "p"])

# 10. Median-split Kaplan-Meier (secondary/descriptive) --------
km_df <- primary_df %>% filter(!is.na(signature_group))
km_fit <- survfit(Surv(OS_time_years, OS_event) ~ signature_group, data = km_df)
logrank <- survdiff(Surv(OS_time_years, OS_event) ~ signature_group, data = km_df)
logrank_p <- pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
km_counts <- km_df %>% group_by(signature_group) %>% summarise(
  patients = n(), deaths = sum(OS_event == 1), censored = sum(OS_event == 0), .groups = "drop")
km_counts$logrank_p <- logrank_p
write_csv(km_counts, "KM_group_counts_and_logrank.csv")
# Build a simple step plot without requiring survminer.
km_sum <- summary(km_fit)
km_plot_df <- data.frame(time_years = km_sum$time, survival = km_sum$surv,
                         strata = as.character(km_sum$strata))
p_km <- ggplot(km_plot_df, aes(x = time_years, y = survival, colour = strata)) +
  geom_step(linewidth = 0.9) + coord_cartesian(ylim = c(0, 1)) +
  labs(title = "TCGA-STAD overall survival by median signature score",
       x = "Overall survival (years)", y = "Overall survival probability", colour = "Signature group") +
  theme_classic(base_size = 12)
ggsave(file.path(output_dir, "TCGA_STAD_signature_KM.pdf"), p_km, width = 6.2, height = 5.2)

# 11. Multivariable Cox: age (per 10 years) + sex + stage, then add signature --------
adj_df <- analysis %>% filter(complete.cases(OS_time_years, OS_event, signature_score_SD,
                                              age_per_10y, sex, stage_group))
if (nrow(adj_df) >= 30 && sum(adj_df$OS_event) >= 10 && nlevels(droplevels(adj_df$sex)) >= 2 &&
    nlevels(droplevels(adj_df$stage_group)) >= 2) {
  adj_df$sex <- droplevels(adj_df$sex)
  adj_df$stage_group <- droplevels(adj_df$stage_group)
  model_cohorts[["adjusted_complete_case"]] <- cohort_count(adj_df, "Adjusted Cox complete cases")
  fit_clinical <- coxph(Surv(OS_time_years, OS_event) ~ age_per_10y + sex + stage_group,
                        data = adj_df, x = TRUE, y = TRUE)
  fit_adjusted <- coxph(Surv(OS_time_years, OS_event) ~ signature_score_SD + age_per_10y + sex + stage_group,
                        data = adj_df, x = TRUE, y = TRUE)
  cox_tables[["clinical_only"]] <- extract_cox(fit_clinical, "Clinical covariates only: age (per 10 years) + sex + stage")
  cox_tables[["adjusted"]] <- extract_cox(fit_adjusted, "Adjusted: signature + age (per 10 years) + sex + stage")
  lrt <- anova(fit_clinical, fit_adjusted, test = "LRT")
  lrt_out <- data.frame(comparison = "Clinical covariates vs clinical covariates + signature",
                        n = fit_adjusted$n, events = fit_adjusted$nevent,
                        chisq = lrt$Chisq[2], df = lrt$Df[2], p_value = lrt$`Pr(>|Chi|)`[2])
  write_csv(lrt_out, "signature_incremental_value_LRT.csv")
  ph_adj <- cox.zph(fit_adjusted)
  ph_tables[["adjusted"]] <- data.frame(model = "Adjusted: signature + age (per 10 years) + sex + stage",
    term = rownames(ph_adj$table), chisq = ph_adj$table[, "chisq"],
    df = ph_adj$table[, "df"], p_value = ph_adj$table[, "p"])
  saveRDS(list(primary = fit_primary, clinical_only = fit_clinical, adjusted = fit_adjusted,
               ph_primary = ph_primary, ph_adjusted = ph_adj),
          file.path(output_dir, "TCGA_STAD_cox_models.rds"))
} else {
  warning("Adjusted model not fitted: insufficient complete cases/events or a covariate has fewer than two observed levels. Check missingness and output cohort tables.")
}

# 12. Gene-wise univariable Cox (15-gene secondary analysis) --------
gene_cox <- lapply(intersect(signature_genes, rownames(expr_z)), function(g) {
  dat <- analysis
  dat$gene_z <- as.numeric(expr_z[g, dat$patient_id])
  dat <- dat %>% filter(complete.cases(OS_time_years, OS_event, gene_z))
  if (nrow(dat) < 20 || sum(dat$OS_event) < 5 || sd(dat$gene_z) == 0) return(NULL)
  fit <- coxph(Surv(OS_time_years, OS_event) ~ gene_z, data = dat)
  s <- summary(fit)
  data.frame(gene = g, HR_per_SD = s$conf.int[1, "exp(coef)"],
    lower95 = s$conf.int[1, "lower .95"], upper95 = s$conf.int[1, "upper .95"],
    p_value = s$coefficients[1, "Pr(>|z|)"], n = fit$n, deaths = fit$nevent,
    censored = fit$n - fit$nevent)
})
gene_cox_df <- bind_rows(gene_cox)
if (nrow(gene_cox_df)) gene_cox_df$FDR_BH <- p.adjust(gene_cox_df$p_value, method = "BH")
write_csv(gene_cox_df, "gene_wise_univariable_cox_BH.csv")

# 13. Signature score by stage; ordered stage correlation --------
stage_df <- analysis %>% filter(!is.na(stage_group), !is.na(signature_score))
if (nrow(stage_df) >= 10 && nlevels(droplevels(stage_df$stage_group)) == 2) {
  stage_wilcox <- wilcox.test(signature_score ~ stage_group, data = stage_df, exact = FALSE)
  stage_summary <- stage_df %>% group_by(stage_group) %>% summarise(
    n = n(), median_score = median(signature_score), IQR_score = IQR(signature_score), .groups = "drop")
  stage_summary$p_wilcoxon <- stage_wilcox$p.value
  write_csv(stage_summary, "signature_score_early_advanced_stage.csv")
  p_stage <- ggplot(stage_df, aes(x = stage_group, y = signature_score)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.12, alpha = 0.45, size = 1.4) +
    labs(title = "TCGA-STAD signature score by pathological stage",
         x = "Pathological stage", y = "Signature score (mean gene-wise z-score)") +
    theme_classic(base_size = 12)
  ggsave(file.path(output_dir, "TCGA_STAD_signature_stage_boxplot.pdf"), p_stage, width = 5.8, height = 5.0)
}

# Map stage to ordinal I-IV for gene-wise Spearman correlations.
stage_ordinal <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z <- gsub("^STAGE\\s*", "", z)
  z <- gsub("^PATHOLOGIC\\s*", "", z)
  z <- gsub("^AJCC\\s*", "", z)
  out <- rep(NA_real_, length(z))
  out[grepl("^I($|A$|B$|C$)", z)] <- 1
  out[grepl("^II($|A$|B$|C$)", z)] <- 2
  out[grepl("^III($|A$|B$|C$)", z)] <- 3
  out[grepl("^IV($|A$|B$|C$)", z)] <- 4
  out
}
analysis$stage_ordinal <- stage_ordinal(analysis$stage_raw)
stage_cor <- lapply(intersect(signature_genes, rownames(expr_z)), function(g) {
  gz <- as.numeric(expr_z[g, analysis$patient_id])
  ok <- is.finite(gz) & is.finite(analysis$stage_ordinal)
  if (sum(ok) < 10 || length(unique(analysis$stage_ordinal[ok])) < 2) return(NULL)
  ct <- suppressWarnings(cor.test(gz[ok], analysis$stage_ordinal[ok], method = "spearman", exact = FALSE))
  data.frame(gene = g, spearman_rho = unname(ct$estimate), p_value = ct$p.value, n = sum(ok))
})
stage_cor_df <- bind_rows(stage_cor)
if (nrow(stage_cor_df)) stage_cor_df$FDR_BH <- p.adjust(stage_cor_df$p_value, method = "BH")
write_csv(stage_cor_df, "gene_wise_spearman_stage_BH.csv")

# 14. Exploratory advanced-stage survival subgroup ----------------
advanced_df <- analysis %>% filter(stage_group == "Advanced",
  complete.cases(OS_time_years, OS_event, signature_score_SD))
if (nrow(advanced_df) >= 20 && sum(advanced_df$OS_event) >= 5) {
  fit_adv <- coxph(Surv(OS_time_years, OS_event) ~ signature_score_SD, data = advanced_df)
  adv_res <- extract_cox(fit_adv, "Exploratory advanced-stage subgroup: continuous signature")
  adv_res$interpretation <- "Exploratory; not adjusted for other covariates; account for subgroup/multiple testing."
  write_csv(adv_res, "exploratory_advanced_stage_cox.csv")
  adv_median <- median(advanced_df$signature_score, na.rm = TRUE)
  advanced_df$signature_group_adv <- factor(ifelse(advanced_df$signature_score >= adv_median, "High", "Low"),
                                             levels = c("Low", "High"))
  adv_lr <- survdiff(Surv(OS_time_years, OS_event) ~ signature_group_adv, data = advanced_df)
  adv_p <- pchisq(adv_lr$chisq, df = length(adv_lr$n) - 1, lower.tail = FALSE)
  write_csv(data.frame(analysis = "Exploratory advanced-stage median split", n = nrow(advanced_df),
    deaths = sum(advanced_df$OS_event == 1), censored = sum(advanced_df$OS_event == 0),
    logrank_p = adv_p), "exploratory_advanced_stage_KM_logrank.csv")
} else {
  write_csv(data.frame(message = "Advanced-stage subgroup not analysed because of insufficient complete cases/events.",
                        n = nrow(advanced_df), deaths = sum(advanced_df$OS_event == 1)),
            "exploratory_advanced_stage_not_run.csv")
}

# 15. Save all model-level counts and results ------------------
model_cohort_df <- bind_rows(model_cohorts)
cox_result_df <- bind_rows(cox_tables)
ph_result_df <- bind_rows(ph_tables)
write_csv(model_cohort_df, "model_specific_patient_death_counts.csv")
write_csv(cox_result_df, "cox_model_results.csv")
write_csv(ph_result_df, "proportional_hazards_tests.csv")

# Include explicit survival-event definition and run metadata.
metadata <- data.frame(
  item = c("OS endpoint", "Event definition", "Censoring definition", "Survival-time units in source",
           "Survival-time units used in models", "Expression transformation", "Signature score",
           "Score scaling in Cox models", "Minimum genes required", "Primary tumour filter",
           "Primary model", "Adjusted model", "Advanced-stage analysis"),
  value = c("Overall survival", "OS_event = 1 for death/deceased/dead or numeric 1",
            "OS_event = 0 for alive/living or numeric 0",
            "Assumed days from the UCSC Xena OS.time field; verify against source file documentation",
            "Years (source time divided by 365.25)", "log2(TPM + 1)",
            "Mean of gene-wise z-scores across matched patients for mapped signature genes",
            "Composite score standardized to 1 SD; HR is per 1-SD increase",
            as.character(minimum_genes_required), "TCGA sample-type code 01 when identifiable",
            "Univariable Cox: signature score per SD", "Cox: signature + age (per 10 years) + sex + early/advanced stage",
            "Exploratory; report separately and account for multiplicity"),
  stringsAsFactors = FALSE
)
write_csv(metadata, "analysis_definitions_and_metadata.csv")

# Print key outputs to the console for the manuscript and supervisor response.
cat("\n================ ANALYSIS COMPLETE ================\n")
cat("Output directory: ", output_dir, "\n", sep = "")
cat("Signature genes mapped: ", length(present_genes), "/15\n", sep = "")
cat("Matched patients: ", nrow(analysis), "\n", sep = "")
cat("Matched deaths: ", sum(analysis$OS_event == 1, na.rm = TRUE), "\n", sep = "")
cat("Matched censored: ", sum(analysis$OS_event == 0, na.rm = TRUE), "\n", sep = "")
cat("\nCohort flow:\n"); print(flow)
cat("\nModel-specific patient/death counts:\n"); print(model_cohort_df)
cat("\nCox results:\n"); print(cox_result_df)
cat("\nPH tests:\n"); print(ph_result_df)
cat("\nImportant: inspect the input survival columns and confirm OS.time is in days and OS coding is 1=death, 0=censored before finalizing manuscript wording.\n")
