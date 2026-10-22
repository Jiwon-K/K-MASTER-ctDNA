
# Figure 4: ctDNA-derived molecular prognostic features ------------------
# Code accompanying the K-MASTER pan-cancer ctDNA study.
# The script reproduces principal summaries and plots from prepared inputs.
# Cohort curation, detailed OncoKB matching, resampling and cross-validation
# are performed upstream; those patient-level inputs are controlled-access.
# Input filenames below describe expected local files, not public downloads.

# Packages ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)
library(survival)
library(survminer)
library(ggrepel)

# Input files ------------------------------------------------------------
# One record per selected patient; see required columns below.
patient_file <- "data/Fig4/Figure4_patient_level.tsv"
# OncoKB v5.0 alteration-specific Cox results (tumor and evidence matching).
oncokb_file <- "data/Fig4/Figure4_OncoKB_Cox.tsv"
# Iteration-level Cox results from 100 repeated random subsamples.
mpi_file <- "data/Fig4/Figure4_MPI_iterations.tsv"
# Selected (non-zero) coefficients from 100 elastic-net iterations.
enet_file <- "data/Fig4/Figure4_elastic_net_iterations.tsv"
# Long-format presence/absence of selected LC features, per patient.
lc_features_file <- "data/Fig4/Figure4_LC_feature_matrix.tsv"
# Independent test-set outcomes and risk scores from the selected CV model.
test_file <- "data/Fig4/Figure4_test_risk_scores.tsv"
# Feature-level SHAP summaries from the selected prognostic model.
shap_file <- "data/Fig4/Figure4_SHAP_summary.tsv"


# Input files ------------------------------------------------------------
# Filtered AXEN ctDNA mutation calls.
maf_file <- "KM00_AXEN_FILTERED_MAF_CONCAT_final_250131.maf.gz"

# K-MASTER clinical annotations, survival, age, stage and subtype.
clinical_file <- "KM00_MASTER_DATASET_240726/02.3.KM00_CLINICAL_onlyONEvisit_OS_MONTH+panel+oncotree.xlsx"

# Original 87-gene panel list (original absolute path is retained for local verification).
panel_file <- "/Users/jiwonkim/Desktop/BLCA:UTUC/MAF/K-MASTER_panels_genelist.xlsx"

# Original curated OncoKB annotation workbook, sheet 'v5.0'.
oncokb_file <- "oncokb_biomarker_drug_associations_KM00_JSH_final.xlsx"

# Intermediate files used explicitly in the original R script.
# axen_final3.rds was generated from the filtered mutation and clinical data.
processed_maf_file <- "axen_final3.rds"

# The 44,191-row result table from 100 repeated MPI analyses.
mpi_file <- "score_axen_minmut2+mut3+onco+alt+znor_1pat1sam_44191.txt"

# Additional intermediate file used in later treatment-related model experiments.
platinum_file <- "axen_1l_cleaned.Rds"

# External validation input files used near the end of the original script.
external_maf_file <- "nsclc_ctdx_msk_2022/data_mutations.txt"
external_clinical_file <- "nsclc_ctdx_msk_2022/data_clinical_sample_merge.txt"


# Data source checks -----------------------------------------------------
source_files <- c(MAF = maf_file, Clinical = clinical_file, Panel = panel_file,
                  OncoKB = oncokb_file, Processed_MAF = processed_maf_file,
                  MPI_100_iterations = mpi_file, Platinum = platinum_file,
                  External_MAF = external_maf_file, External_Clinical = external_clinical_file)

source_check <- data.frame(source = names(source_files), file = unname(source_files),
                           exists = file.exists(source_files), row.names = NULL)
print(source_check, row.names = FALSE)

# Load verified source data -----------------------------------------------
# These input formats and filenames occur in the original script.
# This file checks the inputs; it does not substitute for the modeling and
# plotting steps that generated the original Figure 4D and 4F-I.

if (!file.exists(maf_file)) stop("Missing original MAF: ", maf_file)
if (!file.exists(clinical_file)) stop("Missing original clinical workbook: ", clinical_file)

axen_final2 <- read.delim(maf_file, sep = "\t", check.names = FALSE) %>%
  filter(!Tumor %in% c("Double cancer", "Other"))

selected_tumors <- axen_final2 %>% distinct(Tumor, `SUBJECT ID`) %>%
  dplyr::count(Tumor, name = "n_subjects") %>% filter(n_subjects >= 10) %>% pull(Tumor)
axen_filtered <- axen_final2 %>% filter(Tumor %in% selected_tumors)

cli <- readxl::read_xlsx(clinical_file, skip = 1)
cli$Diagnosis_date_edit <- as.Date(suppressWarnings(as.numeric(cli$Diagnosis_date_edit)),
                                   origin = "1899-12-30")

message("K-MASTER cohort: ", n_distinct(axen_filtered$Tumor), " tumor types, ",
        n_distinct(axen_filtered$`SUBJECT ID`), " patients, ",
        n_distinct(axen_filtered$Tumor_Sample_Barcode), " samples")

if (file.exists(processed_maf_file)) {
  axen_final3 <- readRDS(processed_maf_file)
  message("Loaded original preprocessed MAF: ", processed_maf_file)
} else {
  message("axen_final3.rds is not present; original cohort-level preprocessing is required.")
}

if (file.exists(mpi_file)) {
  random100_all2 <- read.delim(mpi_file, sep = "\t", check.names = FALSE)
  message("Loaded MPI iteration table: ", nrow(random100_all2), " rows")
  print(names(random100_all2))
} else {
  message("The original 100-iteration MPI table is not present: ", mpi_file)
}

if (file.exists(oncokb_file)) {
  onco <- readxl::read_xlsx(oncokb_file, sheet = "v5.0")
  message("Loaded original OncoKB v5.0 sheet: ", nrow(onco), " rows")
}



outdir <- "results/Fig4"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_input <- function(path, cols) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  dat <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(cols, names(dat))
  if (length(missing)) stop(basename(path), ": missing columns: ", paste(missing, collapse = ", "))
  dat
}



# Figure 4A-B: cohort survival and clinical characteristics ---------------
pt <- read_input(patient_file, c("SUBJECT ID", "Tumor", "OS_MONTH_250417", "Survival_250417",
                                 "NGS_stage_edit", "Gender", "Age_edit", "count_edit", "count")) %>%
  distinct(`SUBJECT ID`, .keep_all = TRUE) %>%
  mutate(OS_MONTH_250417 = as.numeric(OS_MONTH_250417), count = as.numeric(count), Survival_250417 = as.integer(Survival_250417))

plot_pt <- pt %>% filter(!is.na(OS_MONTH_250417)) %>% arrange(desc(OS_MONTH_250417)) %>%
  mutate(patient_order = row_number(), Survival = if_else(Survival_250417 == 1, "Death", "Survival"))

p4a_os <- ggplot(plot_pt, aes(patient_order, OS_MONTH_250417, fill = Tumor)) +
  geom_col(width = 1) + theme_classic() +
  labs(x = NULL, y = "Overall survival (months)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")

anno <- plot_pt %>% transmute(patient_order, Stage = as.character(NGS_stage_edit), Survival,
                              Sex = Gender, Age = Age_edit, `ctDNA count` = count_edit) %>%
  pivot_longer(-patient_order, names_to = "Track", values_to = "Value") %>%
  mutate(Track = factor(Track, levels = rev(c("Stage", "Survival", "Sex", "Age", "ctDNA count"))))
p4a_anno <- ggplot(anno, aes(patient_order, Track, fill = Value)) +
  geom_tile(width = 1, height = 0.92) + theme_void() +
  theme(axis.text.y = element_text(), legend.position = "right")

p4a <- p4a_os / p4a_anno + plot_layout(heights = c(3, 1))
ggsave(file.path(outdir, "Figure4A_patient_survival.pdf"), p4a, width = 13, height = 5)


p4b_df <- pt %>% transmute(Stage = as.character(NGS_stage_edit),
                           Survival = if_else(Survival_250417 == 1, "Death", "Survival"), 
                           Sex = Gender, Age = Age_edit) %>%
  pivot_longer(everything(), names_to = "Characteristic", values_to = "Category") %>%
  filter(!is.na(Category)) %>% count(Characteristic, Category) %>%
  group_by(Characteristic) %>% mutate(proportion = n / sum(n)) %>% ungroup()
p4b <- ggplot(p4b_df, aes(x = 1, y = proportion, fill = Category)) +
  geom_col(width = 1, color = "white") + coord_polar(theta = "y") +
  facet_wrap(~Characteristic, ncol = 2) + theme_void() +
  labs(fill = "Category")
ggsave(file.path(outdir, "Figure4B_clinical_pies.pdf"), p4b, width = 8, height = 6)



# Figure 4C: Cox association with continuous mutation count --------------
# HR represents one additional ctDNA mutation, as in the original analysis.
cox_count <- pt %>% filter(!is.na(OS_MONTH_250417), OS_MONTH_250417 > 0,
                           !is.na(Survival_250417), !is.na(count))

cox_by_tumor <- lapply(split(cox_count, cox_count$Tumor), function(dat) {
  if (n_distinct(dat$count) < 2 || sum(dat$Survival_250417 == 1) < 2) return(NULL)
  fit <- tryCatch(coxph(Surv(OS_MONTH_250417, Survival_250417) ~ count, data = dat),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(Tumor = first(dat$Tumor), HR = s$conf.int[1, "exp(coef)"],
             lower = s$conf.int[1, "lower .95"], upper = s$conf.int[1, "upper .95"],
             pval = s$coefficients[1, "Pr(>|z|)"])
}) %>% bind_rows() %>% filter(is.finite(HR), is.finite(lower), is.finite(upper), lower > 0)

p4c <- ggplot(cox_by_tumor, aes(HR, fct_reorder(Tumor, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
  geom_segment(aes(x = lower, xend = upper, yend = Tumor)) + geom_point(size = 2) +
  scale_x_log10() + theme_classic() + labs(x = "HR per mutation", y = NULL)
ggsave(file.path(outdir, "Figure4C_mutation_count_cox.pdf"), p4c, width = 6, height = 5)



# Figure 4D: OncoKB v5.0 alterations and OS ------------------------------
# Input is the outcome of the original tumor/alteration-level annotation and Cox models.
# Expected: Tumor, feature, Level (1/2/3), HR, pval, shape_flag (match/other).
onco <- read_input(oncokb_file, c("Tumor", "feature", "Level", "HR", "pval", "shape_flag")) %>%
  mutate(Level = factor(paste0("L", sub("^L", "", Level)), levels = c("L1", "L2", "L3")),
         neglogp = -log10(pmax(as.numeric(pval), .Machine$double.xmin)),
         log2HR = log2(as.numeric(HR))) %>%
  filter(!is.na(Level), is.finite(log2HR), is.finite(neglogp))

p4d <- ggplot(onco, aes(feature, Tumor)) +
  geom_point(aes(size = neglogp, fill = log2HR, shape = shape_flag), color = "gray25") +
  scale_shape_manual(values = c(match = 21, other = 22)) +
  scale_fill_gradient2(low = "coral", mid = "white", high = "blue", midpoint = 0) +
  facet_grid(~Level, scales = "free_x", space = "free_x") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = NULL, size = "-log10(P)", fill = "log2(HR)", shape = "Annotation")
ggsave(file.path(outdir, "Figure4D_OncoKB_Cox.pdf"), p4d, width = 13, height = 5)



# Figure 4E: Molecular Prognostic Index (100 iterations) -----------------
# The upstream analyses calculated tumor-specific univariable Cox estimates.
# Prevalence is the fraction of feature-positive patients (range: 0-1).
mpi <- read_input(mpi_file, c("Tumor", "iter", "label", "HR", "cox_p", "prevalence")) %>%
  mutate(HR = as.numeric(HR), cox_p = as.numeric(cox_p), prevalence = as.numeric(prevalence)) %>%
  filter(is.finite(HR), HR > 0, is.finite(cox_p), cox_p > 0) %>%
  mutate(MPI = -log2(HR) * -log10(cox_p) * prevalence) %>%
  group_by(Tumor, iter) %>%
  mutate(score_z_iter = if (sd(MPI, na.rm = TRUE) > 0) as.numeric(scale(MPI)) else NA_real_) %>% ungroup()

rank_mpi <- mpi %>% group_by(label) %>% summarise(median_MPI = median(score_z_iter, na.rm = TRUE), .groups = "drop")
selected_mpi <- bind_rows(rank_mpi %>% slice_max(median_MPI, n = 30),
                          rank_mpi %>% slice_min(median_MPI, n = 30)) %>% distinct(label)
p4e_df <- mpi %>% semi_join(selected_mpi, by = "label") %>%
  left_join(rank_mpi, by = "label") %>% mutate(Direction = if_else(median_MPI >= 0, "Favorable", "Unfavorable"))
p4e <- ggplot(p4e_df, aes(fct_reorder(label, median_MPI), score_z_iter, fill = Tumor)) +
  geom_violin(scale = "width", trim = FALSE) + geom_boxplot(width = 0.15, outlier.size = 0.3) +
  facet_wrap(~Direction, scales = "free_x", ncol = 1) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(x = NULL, y = "MPI (z score)")
ggsave(file.path(outdir, "Figure4E_MPI.pdf"), p4e, width = 13, height = 8)



# Figure 4F-G: repeated elastic-net selection ----------------------------
# Expected: Tumor, iter, feature, coefficient (non-zero Cox model coefficient).
# For display, reversed Cox coefficient sign follows the original figure code.
enet <- read_input(enet_file, c("Tumor", "iter", "feature", "coefficient")) %>%
  mutate(coefficient = as.numeric(coefficient)) %>% filter(is.finite(coefficient), coefficient != 0)
enet_summary <- enet %>% group_by(Tumor, feature) %>%
  summarise(times = n_distinct(iter), mean_weight = -mean(coefficient), .groups = "drop")

p4f <- ggplot(enet_summary, aes(mean_weight, times)) +
  geom_point(aes(fill = Tumor), shape = 21, color = "gray30", size = 2.5) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text_repel(data = subset(enet_summary, times >= 50), aes(label = feature), size = 2.7, max.overlaps = 50) + 
  theme_classic() +
  labs(x = "Averaged weights", y = "Selection frequency")
ggsave(file.path(outdir, "Figure4F_elastic_net.pdf"), p4f, width = 8, height = 6)

lc_coef <- enet_summary %>% filter(Tumor == "Lung cancer", times >= 50) %>%
  arrange(mean_weight) %>% mutate(feature = factor(feature, levels = feature))
p4g_coef <- ggplot(lc_coef, aes(mean_weight, feature, fill = mean_weight >= 0)) +
  geom_col() + scale_fill_manual(values = c("TRUE" = "darkorange", "FALSE" = "saddlebrown"), guide = "none") +
  geom_text(aes(label = paste0("(", times, ")")), hjust = if_else(lc_coef$mean_weight > 0, -0.15, 1.1), size = 3) +
  theme_classic() + labs(x = "Mean elastic-net coefficient (sign reversed)", y = NULL)
ggsave(file.path(outdir, "Figure4G_LC_coefficients.pdf"), p4g_coef, width = 6, height = 5)


# The LC patient-by-feature heatmap and OS track use the selected-feature table.
# Expected: SUBJECT ID, OS_MONTH_250417, feature, present (0/1).
lc <- read_input(lc_features_file, c("SUBJECT ID", "OS_MONTH_250417", "feature", "present")) %>%
  filter(feature %in% as.character(lc_coef$feature))
lc_order <- lc %>% distinct(`SUBJECT ID`, OS_MONTH_250417) %>% arrange(OS_MONTH_250417) %>% mutate(patient = row_number())
lc_plot <- lc %>% left_join(lc_order, by = c("SUBJECT ID", "OS_MONTH_250417"))
p4g_heat <- ggplot(lc_plot, aes(patient, feature, fill = factor(present))) +
  geom_tile() + 
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +
  theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, y = NULL, fill = "Altered")

p4g_surv <- ggplot(lc_order, aes(patient, OS_MONTH_250417)) + geom_line() +
  theme_classic() + 
  labs(x = "Patient", y = "OS (months)")
ggsave(file.path(outdir, "Figure4G_LC_patient_features.pdf"), 
       p4g_heat / p4g_surv + plot_layout(heights = c(3, 1)), width = 8, height = 5)



# Figure 4H: independent test-set overall survival -----------------------
# Risk scores are obtained upstream from the best cross-validated model.
test <- read_input(test_file, c("SUBJECT ID", "OS_MONTH_250417", "Survival_250417", "risk_score")) %>%
  distinct(`SUBJECT ID`, .keep_all = TRUE) %>%
  mutate(OS_MONTH_250417 = as.numeric(OS_MONTH_250417),
         Survival_250417 = as.integer(Survival_250417), risk_score = as.numeric(risk_score)) %>%
  mutate(risk_group = factor(if_else(risk_score > median(risk_score, na.rm = TRUE), "High", "Low"), levels = c("Low", "High")))
fit_test <- survfit(Surv(OS_MONTH_250417, Survival_250417) ~ risk_group, data = test)
fit_hr <- coxph(Surv(OS_MONTH_250417, Survival_250417) ~ risk_group, data = test)
p4h <- ggsurvplot(fit_test, data = test, pval = TRUE, risk.table = TRUE,
                  conf.int = FALSE, xlab = "Months", ylab = "Overall survival",
                  legend.title = "Risk group", legend.labs = c("Low", "High"))
pdf(file.path(outdir, "Figure4H_test_set_survival.pdf"), width = 6, height = 6)
print(p4h)
dev.off()



# Figure 4I: SHAP importance --------------------------------------------
# The detailed model fit and SHAP estimation were performed upstream.
# Expected: feature, mean_abs_shap; optional Direction and logrank_p.
shap <- read_input(shap_file, c("feature", "mean_abs_shap")) %>%
  mutate(mean_abs_shap = as.numeric(mean_abs_shap)) %>% arrange(desc(mean_abs_shap)) %>% slice_head(n = 20)
if (!"Direction" %in% names(shap)) shap$Direction <- "Not specified"
p4i <- ggplot(shap, aes(mean_abs_shap, fct_reorder(feature, mean_abs_shap), fill = Direction)) +
  geom_col() + theme_classic() + labs(x = "Mean |SHAP value|", y = NULL)
ggsave(file.path(outdir, "Figure4I_SHAP_importance.pdf"), p4i, width = 6, height = 5)



# Figure 4: ctDNA-derived prognostic features ----------------------------
# K-MASTER pan-cancer ctDNA study. Analysis/plotting code for Figure 4A-I.
# The input files were used in the original Figure4(2).R analysis; they are
# controlled-access or study-specific and are not provided in this repository.
# The saved MPI table is an upstream output of 100 repeated Cox analyses.
# Model and figure output must be compared with the submitted figure before release.

# Packages ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(forcats)
library(purrr)
library(stringr)
library(survival)
library(survminer)
library(glmnet)
library(ggrepel)
library(patchwork)
library(kernelshap)

# Input files: original analysis filenames -------------------------------
# The MAF and clinical workbook are upstream sources used to create
# axen_final3.rds. This script reads the saved RDS directly.
maf_file <- "KM00_AXEN_FILTERED_MAF_CONCAT_final_250131.maf.gz"
clinical_file <- "KM00_MASTER_DATASET_240726/02.3.KM00_CLINICAL_onlyONEvisit_OS_MONTH+panel+oncotree.xlsx"
selected_sample_file <- "KM00_MASTER_DATASET_240726/03.3.KM00_NGS_FINAL_METADATA_FINALannotation_axenserial.xlsx"
panel_file <- "/Users/jiwonkim/Desktop/BLCA:UTUC/MAF/K-MASTER_panels_genelist.xlsx"
oncokb_file <- "oncokb_biomarker_drug_associations_KM00_JSH_final.xlsx"
processed_maf_file <- "axen_final3.rds"
mpi_file <- "score_axen_minmut2+mut3+onco+alt+znor_1pat1sam_44191.txt"
outdir <- "results/Figure4"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# The preprocessed RDS is loaded in the original Figure 4 workflow and
# retains patient-level clinical annotations used by later analyses.
required <- c(processed_maf_file, selected_sample_file, panel_file, oncokb_file, mpi_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Original input file(s) not found:\n", paste(missing, collapse = "\n"))

axen_final3 <- readRDS(processed_maf_file)
needed <- c("SUBJECT ID", "SAMPLE ID", "Tumor", "Tumor_Sample_Barcode", "Hugo_Symbol",
            "HGVSp_Short", "OS_MONTH_250417", "Survival_250417", "count_edit", "count",
            "NGS_stage_edit", "Gender", "Age_edit")
if (!all(needed %in% names(axen_final3))) {
  stop("axen_final3.rds is missing required fields: ", paste(setdiff(needed, names(axen_final3)), collapse = ", "))
}
axen_final3 <- axen_final3 %>% mutate(OS_MONTH_250417 = as.numeric(OS_MONTH_250417),
                                      Survival_250417 = as.numeric(Survival_250417), count = as.numeric(count),
                                      HGVSp_Short_edit = sub("^p\\.", "", as.character(HGVSp_Short)))
message("Figure 4 data: ", n_distinct(axen_final3$`SUBJECT ID`), " patients; ",
        n_distinct(axen_final3$Tumor_Sample_Barcode), " ctDNA samples")

# Select one AXEN sample per patient where required for modeling.
axen_serial <- read_xlsx(selected_sample_file) %>% filter(PANEL == "AXEN", FINAL == "O", is.na(Except))
patient_samples <- axen_final3 %>% distinct(`SUBJECT ID`, `SAMPLE ID`, Tumor_Sample_Barcode)
patient_samples <- patient_samples %>% mutate(preferred = `SAMPLE ID` %in% axen_serial$`SAMPLE ID`) %>%
  arrange(`SUBJECT ID`, desc(preferred)) %>% group_by(`SUBJECT ID`) %>% slice_head(n = 1) %>% ungroup()
axen_one <- axen_final3 %>% semi_join(patient_samples, by = c("SUBJECT ID", "SAMPLE ID", "Tumor_Sample_Barcode"))
stopifnot(n_distinct(axen_one$`SUBJECT ID`) == n_distinct(axen_one$Tumor_Sample_Barcode))

first_valid <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA }
patient <- axen_one %>% group_by(`SUBJECT ID`, Tumor) %>%
  summarise(OS_MONTH_250417 = first_valid(OS_MONTH_250417), Survival_250417 = first_valid(Survival_250417),
            Age_edit = first_valid(Age_edit), Gender = first_valid(Gender),
            NGS_stage_edit = first_valid(NGS_stage_edit), count = first_valid(count),
            count_edit = first_valid(count_edit), .groups = "drop")

cancer_colors <- c("Lung cancer" = "#D55E00", "Breast cancer" = "#0072B2", "Intestinal cancer" = "#009E73",
                   "Gallbladder cancer" = "#F0E442", "Pancreatic cancer" = "#56B4E9", "Gastric cancer" = "#E69F00",
                   "Prostate cancer" = "#B97B60", "Head and neck cancer" = "#D4A10A", "Sarcoma" = "#6C8E3E",
                   "Bladder cancer" = "#E1C16E", "Liver cancer" = "#9C3B3B", "Renal cell cancer" = "#5582B0",
                   "Esophageal cancer" = "#F1B1B1", "Neuroendocrine tumor" = "#AD6F6F", "Ovarian cancer" = "#AABF82",
                   "Thyroid cancer" = "#7F5A6D", "Cervical cancer" = "#4D8AB1", "Melanoma" = "#6D4583", "Thymoma" = "#5F4B8B")

# Figure 4A: overall survival and clinical annotation --------------------
# Original Figure 4A uses one row per ctDNA sample, including 1,243 samples.
fig4a <- axen_final3 %>% distinct(Tumor_Sample_Barcode, .keep_all = TRUE) %>%
  filter(is.finite(OS_MONTH_250417)) %>% arrange(desc(OS_MONTH_250417)) %>% 
  mutate(sample_order = row_number(), Survival = if_else(Survival_250417 == 1, "Death", "Survival"))

p4a_surv <- ggplot(fig4a, aes(sample_order, OS_MONTH_250417, fill = Tumor)) +
  geom_col(width = 1) + scale_fill_manual(values = cancer_colors) +
  theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "Overall survival (months)")

p4a_anno <- fig4a %>% transmute(sample_order, Stage = as.character(NGS_stage_edit), Survival,
                                Sex = Gender, Age = Age_edit, `ctDNA count` = count_edit) %>%
  pivot_longer(-sample_order, names_to = "Track", values_to = "Status") %>%
  mutate(Track = factor(Track, levels = rev(c("Stage", "Survival", "Sex", "Age", "ctDNA count")))) %>%
  ggplot(aes(sample_order, Track, fill = Status)) + geom_tile(height = 0.96) +
  theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, y = NULL)

ggsave(file.path(outdir, "Figure4A.pdf"), p4a_surv / p4a_anno + plot_layout(heights = c(3, 1)), width = 14, height = 6)



# Figure 4B: clinical characteristics ------------------------------------
fig4b <- axen_final3 %>% distinct(Tumor_Sample_Barcode, .keep_all = TRUE) %>%
  transmute(Stage = as.character(NGS_stage_edit), Survival = if_else(Survival_250417 == 1, "Death", "Survival"),
            Sex = Gender, Age = Age_edit) %>% pivot_longer(everything(), names_to = "Characteristic", values_to = "Category") %>%
  filter(!is.na(Category)) %>% dplyr::count(Characteristic, Category) %>% group_by(Characteristic) %>%
  mutate(proportion = n / sum(n)) %>% ungroup()
p4b <- ggplot(fig4b, aes(x = 1, y = proportion, fill = Category)) + geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") + facet_wrap(~Characteristic, ncol = 2) + theme_void()
ggsave(file.path(outdir, "Figure4B.pdf"), p4b, width = 7, height = 6)



# Figure 4C: mutation count and tumor-specific OS ------------------------
hr_count <- axen_final3 %>% filter(is.finite(OS_MONTH_250417), OS_MONTH_250417 > 0,
                                   !is.na(Survival_250417), !is.na(count)) %>% arrange(`SUBJECT ID`, desc(count)) %>%
  distinct(Tumor, `SUBJECT ID`, .keep_all = TRUE) %>% group_split(Tumor) %>% map_dfr(function(dat) {
    if (n_distinct(dat$count) < 2 || sum(dat$Survival_250417 == 1) < 2) return(NULL)
    fit <- tryCatch(coxph(Surv(OS_MONTH_250417, Survival_250417) ~ count, data = dat), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    sm <- summary(fit)
    tibble(Tumor = first(dat$Tumor), HR = sm$conf.int[1, "exp(coef)"],
           lower = sm$conf.int[1, "lower .95"], upper = sm$conf.int[1, "upper .95"],
           p = sm$coef[1, "Pr(>|z|)"], count_max = max(dat$count))
  }) %>% filter(is.finite(HR), is.finite(lower), is.finite(upper), lower > 0, count_max > 3)
p4c <- ggplot(hr_count, aes(HR, fct_reorder(Tumor, HR))) +
  geom_vline(xintercept = 1, color = "gray70", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.1) + geom_point(size = 2.5) +
  scale_x_log10() + theme_classic() + labs(x = "HR per mutation (95% CI)", y = NULL)
ggsave(file.path(outdir, "Figure4C.pdf"), p4c, width = 6, height = 5)



# Figure 4D: OncoKB v5.0 alteration / OS association ---------------------
# The annotation source and matching categories follow the original script.
panel_genes <- read_xlsx(panel_file, sheet = "Sheet1")$`K&C&A_87` %>% as.character() %>% na.omit()
onco <- read_xlsx(oncokb_file, sheet = "v5.0") %>%
  distinct(Level, Gene, Alterations, `Cancer Types`, .keep_all = TRUE) %>%
  filter(Gene %in% panel_genes, Tumor != "Liquid", !is.na(Tumor), is.na(annotation),
         !is.na(Alterations), category %in% c("all_mutation", "specific_mutation") | Alterations == "Wildtype") %>%
  mutate(Level = sub("A$", "", sub("^L", "", as.character(Level))),
         Alterations = str_replace_all(Alterations, "\u2013", "-"))

fig4d_maf <- axen_final3 %>% filter((Tumor == "Lung cancer" & Tumor_subtype_JW %in% c("LUAD", "LUSC", "LCLC", "NSCLCPD")) |
                                      (Tumor == "Intestinal cancer" & Tumor_subtype_JW == "COADREAD") | 
                                      !(Tumor %in% c("Lung cancer", "Intestinal cancer"))) %>%
  mutate(HGVSp_Short_edit = str_replace_all(HGVSp_Short_edit, "\u2013", "-"))
pt4d <- fig4d_maf %>% distinct(`SUBJECT ID`, Tumor, OS_MONTH_250417, Survival_250417) %>%
  filter(is.finite(OS_MONTH_250417), OS_MONTH_250417 > 0, !is.na(Survival_250417))

exon19 <- c("E746_A750del", "L747_S752del", "L747_T751del", "L747_P753delinsS")
feature_label <- function(gene, alt) {
  if (gene == "EGFR" && alt %in% exon19) return("EGFR Exon 19 in-frame deletions")
  paste(gene, ifelse(alt == "all", "", alt)) %>% trimws()
}

match_annotation <- function(rule) {
  tumor <- as.character(rule$Tumor); gene <- as.character(rule$Gene); alt <- as.character(rule$Alterations)
  tumor_subset <- if (tumor == "Pan-cancer") unique(fig4d_maf$Tumor) else tumor
  out <- lapply(tumor_subset, function(tm) {
    dat <- fig4d_maf %>% filter(Tumor == tm)
    if (alt == "Wildtype") {
      if (!(tm == "Intestinal cancer" && gene %in% c("KRAS", "NRAS"))) return(NULL)
      sub_ids <- pt4d %>% filter(Tumor == tm) %>% pull(`SUBJECT ID`)
      hit_ids <- setdiff(sub_ids, dat %>% filter(Hugo_Symbol == gene) %>% pull(`SUBJECT ID`))
    } else if (as.character(rule$category) == "all_mutation") {
      hit_ids <- dat %>% filter(Hugo_Symbol == gene) %>% pull(`SUBJECT ID`)
    } else {
      hit_ids <- dat %>% filter(Hugo_Symbol == gene, HGVSp_Short_edit == alt) %>% pull(`SUBJECT ID`)
    }
    tibble(`SUBJECT ID` = unique(hit_ids), Tumor = tm, feature = feature_label(gene, alt), Level = as.character(rule$Level))
  })
  bind_rows(out)
}

annotated <- map_dfr(seq_len(nrow(onco)), function(i) match_annotation(onco[i, ])) %>%
  filter(Level %in% c("1", "2", "3", "4")) %>% distinct()
feature_level <- annotated %>% group_by(Tumor, feature) %>%
  summarise(Level = as.character(min(as.numeric(Level))), .groups = "drop")
global_level <- feature_level %>% group_by(feature) %>%
  summarise(best_level = as.character(min(as.numeric(Level))), .groups = "drop")

fit_onco <- function(dat, ft) {
  positive <- annotated %>% filter(Tumor == first(dat$Tumor), feature == ft) %>% pull(`SUBJECT ID`) %>% unique()
  dat <- dat %>% mutate(present = as.integer(`SUBJECT ID` %in% positive))
  if (sum(dat$present) < 2 || n_distinct(dat$present) < 2) return(NULL)
  events <- tapply(dat$Survival_250417, dat$present, sum)
  if (length(events) < 2 || any(events < 1)) return(NULL)
  fit <- tryCatch(coxph(Surv(OS_MONTH_250417, Survival_250417) ~ present, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)
  if (any(!is.finite(sm$conf.int[1, c("lower .95", "upper .95")]))) return(NULL)
  tibble(Tumor = first(dat$Tumor), feature = ft,
         HR = sm$conf.int[1, "exp(coef)"], p = sm$coef[1, "Pr(>|z|)"])
}

onco_cox <- map_dfr(split(pt4d, pt4d$Tumor), function(dat) {
  map_dfr(unique(annotated$feature), ~fit_onco(dat, .x))
}) %>% left_join(feature_level, by = c("Tumor", "feature")) %>%
  left_join(global_level, by = "feature") %>%
  filter(best_level %in% c("1", "2", "3"), is.finite(HR), HR > 0, is.finite(p), p > 0) %>%
  mutate(level_band = factor(paste0("L", best_level), levels = c("L1", "L2", "L3")),
         match = if_else(Level == best_level, "match", "other"),
         neglogp = -log10(p), log2HR = -log2(HR))
p4d <- ggplot(onco_cox, aes(feature, Tumor)) +
  geom_point(aes(fill = log2HR, size = neglogp, shape = match, alpha = p <= 0.05), color = "gray30") +
  scale_shape_manual(values = c(match = 21, other = 22)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
  scale_fill_gradient2(low = "coral", mid = "white", high = "blue", midpoint = 0) +
  #facet_grid(~level_band, scales = "free_x", space = "free_x") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "-log2(HR)", size = "-log10(P)")
ggsave(file.path(outdir, "Figure4D.pdf"), p4d)



# Figure 4E: MPI from the original 100-iteration output ------------------
# Do not replace the actual repeated-analysis table with a fabricated TSV.
mpi <- read.delim(mpi_file, sep = "\t", check.names = FALSE)
stopifnot(all(c("label", "score_z_iter", "Tumor", "cox_pval") %in% names(mpi)))
mpi_rank <- mpi %>% group_by(label) %>% summarise(median_score = median(score_z_iter, na.rm = TRUE), .groups = "drop")
show_n <- 25L
selected_mpi <- bind_rows(mpi_rank %>% slice_max(median_score, n = show_n),
                          mpi_rank %>% slice_min(median_score, n = show_n)) %>% distinct(label)
mpi_plot <- mpi %>% semi_join(selected_mpi, by = "label") %>% left_join(mpi_rank, by = "label") %>%
  mutate(Direction = if_else(median_score >= 0, "Favorable", "Unfavorable"))
p4e <- ggplot(mpi_plot, aes(fct_reorder(label, median_score), score_z_iter, fill = Tumor)) +
  geom_violin(scale = "width", trim = FALSE) + geom_boxplot(width = 0.12, outlier.size = 0.3, fill = "white") +
  facet_wrap(~Direction, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = cancer_colors) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) + labs(x = NULL, y = "MPI (z score)")
ggsave(file.path(outdir, "Figure4E.pdf"), p4e, width = 13, height = 9)



# Figure 4F-G: repeated Cox elastic-net ----------------------------------
# Original design: alpha=0.5, lambda.min, 100 repeated 80/90% subsamples.
# This is a new run of the original analysis, not a pre-filled result table.
exclude_genes <- c("ETV1", "ETV4", "ETV5", "FLT3", "GNAQ", "GNAS", "JAK1", "JAK2", "JAK3",
                   "MAP2K1", "MAP2K2", "MAP2K4", "MPL", "NPM1", "RUNX1", "SMARCB1")
level_groups <- list(L1_vs_L4N = c("1"), L12_vs_L4N = c("1", "2"), L123_vs_L4N = c("1", "2", "3"))

make_tumor_matrix <- function(tm) {
  maf_tm <- axen_one %>% filter(Tumor == tm, is.finite(OS_MONTH_250417), !is.na(Survival_250417))
  core <- maf_tm %>% distinct(`SUBJECT ID`, OS_MONTH_250417, Survival_250417, count_edit)
  genes <- maf_tm %>% filter(!Hugo_Symbol %in% exclude_genes) %>%
    distinct(`SUBJECT ID`, Hugo_Symbol) %>% mutate(present = 1L) %>%
    pivot_wider(names_from = Hugo_Symbol, values_from = present, values_fill = 0)
  mat <- core %>% left_join(genes, by = "SUBJECT ID") %>%
    mutate(mutHL = as.integer(count_edit == "mut-H")) %>% dplyr::select(-count_edit)
  has_level <- annotated %>% filter(Tumor == tm) %>% group_by(`SUBJECT ID`) %>%
    summarise(Level = min(as.numeric(Level)), .groups = "drop")
  mat <- mat %>% left_join(has_level, by = "SUBJECT ID")
  for (nm in names(level_groups)) {
    mat[[nm]] <- ifelse(mat$Level %in% as.numeric(level_groups[[nm]]), 1L, ifelse(is.na(mat$Level) | mat$Level >= 4, 0L, NA_integer_))
  }
  mat <- mat %>% dplyr::select(-Level)
  alt <- annotated %>% filter(Tumor == tm, str_detect(feature, " ")) %>%
    distinct(`SUBJECT ID`, feature, Level) %>%
    mutate(feature_id = make.names(paste0("L", Level, "_", feature)), present = 1L) %>%
    dplyr::select(`SUBJECT ID`, feature_id, present) %>% distinct() %>%
    pivot_wider(names_from = feature_id, values_from = present, values_fill = 0)
  if (nrow(alt)) mat <- mat %>% left_join(alt, by = "SUBJECT ID")
  mat %>% mutate(across(-c(`SUBJECT ID`, OS_MONTH_250417, Survival_250417), ~replace_na(.x, 0))) %>%
    distinct(`SUBJECT ID`, .keep_all = TRUE)
}

boot_coefs <- list()
full_mats <- list()
for (tm in unique(axen_one$Tumor)) {
  dat <- make_tumor_matrix(tm)
  if (nrow(dat) < 10 || sum(dat$Survival_250417 == 1, na.rm = TRUE) < 3) next
  full_mats[[tm]] <- dat
  ratio <- if (tm %in% c("Lung cancer", "Breast cancer", "Intestinal cancer")) 0.8 else 0.9
  set.seed(42)
  boot_coefs[[tm]] <- map_dfr(seq_len(100), function(iter) {
    sub <- dat %>% slice_sample(n = floor(nrow(dat) * ratio), replace = FALSE)
    feat <- setdiff(names(sub), c("SUBJECT ID", "OS_MONTH_250417", "Survival_250417"))
    keep <- feat[vapply(sub[feat], function(x) sum(x == 1, na.rm = TRUE) >= 2 && sum(x == 0, na.rm = TRUE) >= 2, logical(1))]
    if (!length(keep)) return(NULL)
    x <- data.matrix(sub[keep]); y <- Surv(sub$OS_MONTH_250417, sub$Survival_250417)
    cv <- tryCatch(glmnet::cv.glmnet(x, y, family = "cox", alpha = 0.5), error = function(e) NULL)
    if (is.null(cv)) return(NULL)
    cf <- as.matrix(coef(cv, s = "lambda.min")); cf <- cf[cf[, 1] != 0, , drop = FALSE]
    if (!nrow(cf)) return(NULL)
    tibble(Tumor = tm, feature = rownames(cf), coefficient = as.numeric(cf[, 1]), iter = iter)
  })
  message("Elastic-net subsampling finished: ", tm)
}
enet <- bind_rows(boot_coefs)
if (!nrow(enet)) stop("No elastic-net features were selected; inspect input and Cox events")
enet_summary <- enet %>% group_by(Tumor, feature) %>%
  summarise(times = n_distinct(iter), coef_mean = mean(coefficient), mean_weight = -mean(coefficient), .groups = "drop")

feature_p <- map_dfr(names(full_mats), function(tm) {
  dat <- full_mats[[tm]]
  fs <- enet_summary %>% filter(Tumor == tm) %>% pull(feature)
  map_dfr(fs, function(f) {
    if (!f %in% names(dat) || n_distinct(dat[[f]]) < 2) return(NULL)
    z <- dat %>% transmute(OS_MONTH_250417, Survival_250417, marker = .data[[f]])
    fit <- tryCatch(coxph(Surv(OS_MONTH_250417, Survival_250417) ~ marker, data = z), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    tibble(Tumor = tm, feature = f, p = summary(fit)$coef[1, "Pr(>|z|)"])
  })
})
plot_enet <- enet_summary %>% left_join(feature_p, by = c("Tumor", "feature")) %>%
  mutate(neglogp = -log10(pmax(p, .Machine$double.xmin)))
p4f <- ggplot(plot_enet, aes(mean_weight, times)) +
  geom_point(aes(size = neglogp, fill = if_else(times >= 50, Tumor, "Other")), shape = 21, alpha = 0.85) +
  geom_text_repel(data = subset(plot_enet, times >= 50), aes(label = paste(Tumor, feature, sep = ": ")),
                  size = 2.5, max.overlaps = 70) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = 50, linetype = "dashed") +
  scale_fill_manual(values = c(cancer_colors, Other = "gray85")) + theme_classic() +
  labs(x = "Averaged weights", y = "Bootstrap frequency", size = "-log10(P)")
ggsave(file.path(outdir, "Figure4F.pdf"), p4f, width = 9, height = 7)


lc <- enet_summary %>% filter(Tumor == "Lung cancer", times >= 50) %>% arrange(mean_weight)
if (nrow(lc)) {
  lc <- lc %>% mutate(feature = factor(feature, levels = feature))
  p4g_coef <- ggplot(lc, aes(mean_weight, feature, fill = coef_mean < 0)) +
    geom_col() + 
    scale_fill_manual(values = c(`TRUE` = "#FF9933", `FALSE` = "#8A3700"), guide = "none") +
    theme_classic() + 
    labs(x = "Mean elastic-net coefficient (sign reversed)", y = NULL)
  lc_mat <- full_mats[["Lung cancer"]] %>% arrange(OS_MONTH_250417) %>% mutate(patient_order = row_number())
  p4g_heat <- lc_mat %>% dplyr::select(patient_order, all_of(as.character(lc$feature))) %>%
    pivot_longer(-patient_order, names_to = "feature", values_to = "present") %>%
    ggplot(aes(patient_order, feature, fill = factor(present))) + geom_tile() +
    scale_fill_manual(values = c(`0` = "white", `1` = "black")) +
    theme_classic() + theme(axis.text.x = element_blank(), legend.position = "none")
  p4g_os <- ggplot(lc_mat, aes(patient_order, OS_MONTH_250417)) + geom_line() +
    theme_classic() + labs(x = "Patient", y = "OS (months)")
  ggsave(file.path(outdir, "Figure4G.pdf"), p4g_coef | (p4g_heat / p4g_os), width = 13, height = 6)
}



# Figure 4H-I: tumor-stratified split, Cox models, and SHAP ---------------
# Original final code uses an 80:20 train/test split and five folds within
# the training cohort (seed 2026031922). H: Model 1, fold 1; I: Model 3.
model_pt <- axen_final3 %>% group_by(`SUBJECT ID`, Tumor) %>%
  summarise(OS_MONTH_250417 = first_valid(OS_MONTH_250417), Survival_250417 = first_valid(Survival_250417),
            Age_edit = first_valid(Age_edit), Gender = first_valid(Gender),
            NGS_stage_edit = first_valid(NGS_stage_edit), count = max(count, na.rm = TRUE),
            count_edit = if_else(any(count_edit == "mut-H", na.rm = TRUE), "mut-H", "mut-L"), .groups = "drop") %>%
  mutate(stage2 = case_when(NGS_stage_edit %in% c("1", "2") ~ "Low",
                            NGS_stage_edit %in% c("3", "4") ~ "High", TRUE ~ NA_character_),
         stage2 = factor(stage2, levels = c("High", "Low")), Gender = factor(Gender, levels = c("M", "F")),
         Age_edit = factor(Age_edit, levels = c("Older", "Younger")),
         count_edit = factor(count_edit, levels = c("mut-H", "mut-L"))) %>%
  filter(is.finite(OS_MONTH_250417), OS_MONTH_250417 > 0, !is.na(Survival_250417),
         !is.na(stage2), !is.na(Age_edit), !is.na(Gender), !is.na(count_edit)) %>% arrange(`SUBJECT ID`)

split_by_tumor <- function(dat, test_prop = 0.2, seed = 2026031922) {
  set.seed(seed)
  splits <- split(dat, dat$Tumor)
  test <- bind_rows(lapply(splits, function(d) {
    n <- nrow(d); if (n <= 1) return(d[FALSE, ])
    n_test <- min(max(1L, round(n * test_prop)), n - 1L)
    idx <- sample.int(n, size = n, replace = FALSE)
    d[idx[seq_len(n_test)], , drop = FALSE]
  }))
  list(train = anti_join(dat, test, by = "SUBJECT ID"), test = test)
}
make_folds <- function(dat, k = 5, seed = 2026031922) {
  set.seed(seed)
  bind_rows(lapply(split(dat, dat$Tumor), function(d) {
    d$fold_id <- sample(rep(seq_len(k), length.out = nrow(d)))
    d
  }))
}

split <- split_by_tumor(model_pt)
train_folds <- make_folds(split$train)
train1 <- train_folds %>% filter(fold_id != 1) %>% dplyr::select(-fold_id)
test_pt <- split$test

# Model 1: clinical covariates, tumor-specific baseline hazards.
fit1 <- coxph(Surv(OS_MONTH_250417, Survival_250417) ~ Age_edit + Gender + stage2 + count_edit + strata(Tumor),
              data = train1, x = TRUE, y = TRUE)
test_pt$risk_score <- as.numeric(predict(fit1, newdata = test_pt, type = "lp"))
test_pt <- test_pt %>% mutate(risk_group = factor(if_else(risk_score > median(risk_score, na.rm = TRUE), "High", "Low"),
                                                  levels = c("Low", "High")))
fit4h <- survfit(Surv(OS_MONTH_250417, Survival_250417) ~ risk_group, data = test_pt)
p4h <- ggsurvplot(fit4h, data = test_pt, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                  surv.median.line = "hv", legend.title = "Risk group", legend.labs = c("Low", "High"), ggtheme = theme_bw())
pdf(file.path(outdir, "Figure4H.pdf"), width = 7, height = 6)
print(p4h)
dev.off()

# Model 3: clinical covariates and patient-level gene mutation indicators.
pt_gene <- axen_final3 %>% distinct(`SUBJECT ID`, Hugo_Symbol) %>% mutate(present = 1L)
gene_freq <- pt_gene %>% semi_join(model_pt, by = "SUBJECT ID") %>%
  dplyr::count(Hugo_Symbol, name = "n") %>% mutate(prevalence = n / nrow(model_pt))
selected_genes <- gene_freq %>% filter(prevalence >= 0.01) %>% pull(Hugo_Symbol)
gene_mat <- pt_gene %>% filter(Hugo_Symbol %in% selected_genes) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = present, values_fill = 0)
names(gene_mat)[-1] <- make.names(names(gene_mat)[-1], unique = TRUE)
model3_pt <- model_pt %>% left_join(gene_mat, by = "SUBJECT ID") %>%
  mutate(across(all_of(make.names(selected_genes)), ~replace_na(.x, 0L)))
train3 <- model3_pt %>% semi_join(train1, by = "SUBJECT ID")
valid3 <- model3_pt %>% semi_join(train_folds %>% filter(fold_id == 1), by = "SUBJECT ID")
formula3 <- reformulate(c("Age_edit", "Gender", "stage2", "count_edit", make.names(selected_genes),
                          "strata(Tumor)"), response = "Surv(OS_MONTH_250417, Survival_250417)")
fit3 <- coxph(formula3, data = train3, x = TRUE, y = TRUE)

# Kernel SHAP as used in the original analysis (not a summary TSV).
x_explain <- as.data.frame(model.matrix(fit3, data = valid3))
coef_vec <- coef(fit3)
feat <- intersect(names(coef_vec), names(x_explain))
x_explain <- x_explain[, feat, drop = FALSE]
coef_vec <- coef_vec[feat]
predict_lp <- function(object, newdata) as.numeric(as.matrix(newdata[, feat, drop = FALSE]) %*% coef_vec)
set.seed(2026031922)
ks <- kernelshap(fit3, X = x_explain, bg_X = x_explain, pred_fun = predict_lp)
shap_summary <- as.data.frame(ks$S) %>% pivot_longer(everything(), names_to = "feature", values_to = "shap") %>%
  group_by(feature) %>% summarise(mean_abs_shap = mean(abs(shap), na.rm = TRUE), .groups = "drop") %>%
  slice_max(mean_abs_shap, n = 20) %>% arrange(mean_abs_shap)
p4i <- ggplot(shap_summary, aes(mean_abs_shap, factor(feature, levels = feature))) +
  geom_col(fill = "gray45") + theme_bw() + labs(x = "Mean |SHAP value|", y = NULL)
ggsave(file.path(outdir, "Figure4I.pdf"), p4i, width = 7, height = 6)



