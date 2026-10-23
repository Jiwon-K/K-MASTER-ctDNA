
# Figure 4: K-MASTER pan-cancer ctDNA analysis ------------------------
# Figure-level analysis and plots using controlled-access K-MASTER input data.

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


# Input files ------------------------------------------------------------
selected_sample_file <- "KM00_CLINICAL_serial.xlsx"
panel_file <- path.expand("K-MASTER_panels_genelist.xlsx")
oncokb_file <- "oncokb_biomarker.xlsx"
processed_maf_file <- read.delim("KM00_AXEN.maf.gz", sep = "\t", check.names = F) 
outdir <- "results/Figure4"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required <- c(processed_maf_file, selected_sample_file, panel_file, oncokb_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Original input file(s) not found:\n", paste(missing, collapse = "\n"))

ctdna_data <- readRDS(processed_maf_file)
needed <- c("SUBJECT ID", "SAMPLE ID", "Tumor", "Tumor_Sample_Barcode", "Hugo_Symbol",
            "HGVSp_Short", "OS_MONTH_250417", "Survival_250417", "count_edit", "count",
            "NGS_stage_edit", "Gender", "Age_edit", "Tumor_subtype_JW")
if (!all(needed %in% names(ctdna_data))) {
  stop("Input ctDNA file is missing required fields: ", paste(setdiff(needed, names(ctdna_data)), collapse = ", "))
}
# Rename source-specific fields once; all subsequent steps use descriptive names.
ctdna_data <- ctdna_data %>%
  dplyr::rename(patient_id = `SUBJECT ID`, sample_id = `SAMPLE ID`,
                sample_barcode = Tumor_Sample_Barcode, gene_symbol = Hugo_Symbol,
                protein_annotation = HGVSp_Short, os_months = OS_MONTH_250417,
                os_event = Survival_250417, mutation_count = count,
                mutation_group = count_edit, clinical_stage = NGS_stage_edit,
                sex = Gender, age_group = Age_edit, tumor_subtype = Tumor_subtype_JW) %>%
  mutate(os_months = as.numeric(os_months), os_event = as.numeric(os_event),
         mutation_count = as.numeric(mutation_count),
         protein_change = sub("^p\\.", "", as.character(protein_annotation)))

message("Figure 4 data: ", n_distinct(ctdna_data$patient_id), " patients; ",
        n_distinct(ctdna_data$sample_barcode), " ctDNA samples")

sample_metadata <- read_xlsx(selected_sample_file) %>%
  dplyr::rename(sample_id = `SAMPLE ID`, panel = PANEL, selected = FINAL, exclusion_flag = Except) %>%
  filter(panel == "AXEN", selected == "O", is.na(exclusion_flag))
patient_samples <- ctdna_data %>% distinct(patient_id, sample_id, sample_barcode)
patient_samples <- patient_samples %>% mutate(preferred = sample_id %in% sample_metadata$sample_id) %>%
  arrange(patient_id, desc(preferred)) %>% group_by(patient_id) %>% slice_head(n = 1) %>% ungroup()
selected_ctdna <- ctdna_data %>% semi_join(patient_samples, by = c("patient_id", "sample_id", "sample_barcode"))
stopifnot(n_distinct(selected_ctdna$patient_id) == n_distinct(selected_ctdna$sample_barcode))

first_valid <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA }

cancer_colors <- c("Lung cancer" = "#D55E00", "Breast cancer" = "#0072B2", "Intestinal cancer" = "#009E73",
                   "Gallbladder cancer" = "#F0E442", "Pancreatic cancer" = "#56B4E9", "Gastric cancer" = "#E69F00",
                   "Prostate cancer" = "#B97B60", "Head and neck cancer" = "#D4A10A", "Sarcoma" = "#6C8E3E",
                   "Bladder cancer" = "#E1C16E", "Liver cancer" = "#9C3B3B", "Renal cell cancer" = "#5582B0",
                   "Esophageal cancer" = "#F1B1B1", "Neuroendocrine tumor" = "#AD6F6F", "Ovarian cancer" = "#AABF82",
                   "Thyroid cancer" = "#7F5A6D", "Cervical cancer" = "#4D8AB1", "Melanoma" = "#6D4583", "Thymoma" = "#5F4B8B")



# Figure 4A: overall survival and clinical annotation --------------------
survival_samples <- ctdna_data %>% distinct(sample_barcode, .keep_all = TRUE) %>%
  filter(is.finite(os_months)) %>% arrange(desc(os_months)) %>% mutate(sample_order = row_number(),
                                                                       Survival = if_else(os_event == 1, "Death", "Survival"))

p4a_surv <- ggplot(survival_samples, aes(sample_order, os_months, fill = Tumor)) +
  geom_col(width = 1) + scale_fill_manual(values = cancer_colors) +
  theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "Overall survival (months)")

p4a_anno <- survival_samples %>% transmute(sample_order, Stage = as.character(clinical_stage), Survival,
                                           Sex = sex, Age = age_group, `ctDNA count` = mutation_group) %>%
  pivot_longer(-sample_order, names_to = "Track", values_to = "Status") %>%
  mutate(Track = factor(Track, levels = rev(c("Stage", "Survival", "Sex", "Age", "ctDNA count")))) %>%
  ggplot(aes(sample_order, Track, fill = Status)) + geom_tile(height = 0.96) +
  theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + labs(x = NULL, y = NULL)

ggsave(file.path(outdir, "Figure4A.pdf"), p4a_surv / p4a_anno + plot_layout(heights = c(3, 1)), width = 14, height = 6)



# Figure 4B: clinical characteristics ------------------------------------
clinical_summary <- ctdna_data %>% distinct(sample_barcode, .keep_all = TRUE) %>%
  transmute(Stage = as.character(clinical_stage), Survival = if_else(os_event == 1, "Death", "Survival"),
            Sex = sex, Age = age_group) %>% pivot_longer(everything(), names_to = "Characteristic", values_to = "Category") %>%
  filter(!is.na(Category)) %>% dplyr::count(Characteristic, Category) %>% group_by(Characteristic) %>%
  mutate(proportion = n / sum(n)) %>% ungroup()
p4b <- ggplot(clinical_summary, aes(x = 1, y = proportion, fill = Category)) + geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") + facet_wrap(~Characteristic, ncol = 2) + theme_void()
ggsave(file.path(outdir, "Figure4B.pdf"), p4b, width = 7, height = 6)



# Figure 4C: mutation count and tumor-specific OS ------------------------
mutation_cox_results <- ctdna_data %>% filter(is.finite(os_months), os_months > 0,
                                              !is.na(os_event), !is.na(mutation_count)) %>% arrange(patient_id, desc(mutation_count)) %>%
  distinct(Tumor, patient_id, .keep_all = TRUE) %>% group_split(Tumor) %>% map_dfr(function(dat) {
    if (n_distinct(dat$mutation_count) < 2 || sum(dat$os_event == 1) < 2) return(NULL)
    fit <- tryCatch(coxph(Surv(os_months, os_event) ~ mutation_count, data = dat), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    sm <- summary(fit)
    tibble(Tumor = first(dat$Tumor), HR = sm$conf.int[1, "exp(coef)"],
           lower = sm$conf.int[1, "lower .95"], upper = sm$conf.int[1, "upper .95"],
           p = sm$coef[1, "Pr(>|z|)"], count_max = max(dat$mutation_count))
  }) %>% filter(is.finite(HR), is.finite(lower), is.finite(upper), lower > 0, count_max > 3)
p4c <- ggplot(mutation_cox_results, aes(HR, fct_reorder(Tumor, HR))) +
  geom_vline(xintercept = 1, color = "gray70", linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.1) + geom_point(size = 2.5) +
  scale_x_log10() + theme_classic() + labs(x = "HR per mutation (95% CI)", y = NULL)
ggsave(file.path(outdir, "Figure4C.pdf"), p4c, width = 6, height = 5)



# Figure 4D: OncoKB v5.0 alteration / OS association ---------------------
panel_genes <- read_xlsx(panel_file, sheet = "Sheet1")$`K&C&A_87` %>% as.character() %>% na.omit()
oncokb_annotations <- read_xlsx(oncokb_file, sheet = "v5.0") %>%
  distinct(Level, Gene, Alterations, `Cancer Types`, .keep_all = TRUE) %>%
  filter(Gene %in% panel_genes, Tumor != "Liquid", !is.na(Tumor), is.na(annotation),
         !is.na(Alterations), category %in% c("all_mutation", "specific_mutation") | Alterations == "Wildtype") %>%
  mutate(Level = sub("A$", "", sub("^L", "", as.character(Level))),
         Alterations = str_replace_all(Alterations, "\u2013", "-"))

oncokb_ctdna <- ctdna_data %>% filter((Tumor == "Lung cancer" & tumor_subtype %in% c("LUAD", "LUSC", "LCLC", "NSCLCPD")) |
                                        (Tumor == "Intestinal cancer" & tumor_subtype == "COADREAD") | !(Tumor %in% c("Lung cancer", "Intestinal cancer"))) %>%
  mutate(protein_change = str_replace_all(protein_change, "\u2013", "-"))
oncokb_patients <- oncokb_ctdna %>% distinct(patient_id, Tumor, os_months, os_event) %>%
  filter(is.finite(os_months), os_months > 0, !is.na(os_event))

exon19 <- c("E746_A750del", "L747_S752del", "L747_T751del", "L747_P753delinsS")
feature_label <- function(gene, alt) {
  if (gene == "EGFR" && alt %in% exon19) return("EGFR Exon 19 in-frame deletions")
  paste(gene, ifelse(alt == "all", "", alt)) %>% trimws()
}

match_annotation <- function(rule) {
  tumor <- as.character(rule$Tumor); gene <- as.character(rule$Gene); alt <- as.character(rule$Alterations)
  tumor_subset <- if (tumor == "Pan-cancer") unique(oncokb_ctdna$Tumor) else tumor
  out <- lapply(tumor_subset, function(tm) {
    dat <- oncokb_ctdna %>% filter(Tumor == tm)
    if (alt == "Wildtype") {
      if (!(tm == "Intestinal cancer" && gene %in% c("KRAS", "NRAS"))) return(NULL)
      sub_ids <- oncokb_patients %>% filter(Tumor == tm) %>% pull(patient_id)
      hit_ids <- setdiff(sub_ids, dat %>% filter(gene_symbol == gene) %>% pull(patient_id))
    } else if (as.character(rule$category) == "all_mutation") {
      hit_ids <- dat %>% filter(gene_symbol == gene) %>% pull(patient_id)
    } else {
      hit_ids <- dat %>% filter(gene_symbol == gene, protein_change == alt) %>% pull(patient_id)
    }
    tibble(patient_id = unique(hit_ids), Tumor = tm, feature = feature_label(gene, alt), Level = as.character(rule$Level))
  })
  bind_rows(out)
}

annotated <- map_dfr(seq_len(nrow(oncokb_annotations)), function(i) match_annotation(oncokb_annotations[i, ])) %>%
  filter(Level %in% c("1", "2", "3", "4")) %>% distinct()
feature_level <- annotated %>% group_by(Tumor, feature) %>%
  summarise(Level = as.character(min(as.numeric(Level))), .groups = "drop")
global_level <- feature_level %>% group_by(feature) %>%
  summarise(best_level = as.character(min(as.numeric(Level))), .groups = "drop")

fit_onco <- function(dat, ft) {
  positive <- annotated %>% filter(Tumor == first(dat$Tumor), feature == ft) %>% pull(patient_id) %>% unique()
  dat <- dat %>% mutate(present = as.integer(patient_id %in% positive))
  if (sum(dat$present) < 2 || n_distinct(dat$present) < 2) return(NULL)
  events <- tapply(dat$os_event, dat$present, sum)
  if (length(events) < 2 || any(events < 1)) return(NULL)
  fit <- tryCatch(coxph(Surv(os_months, os_event) ~ present, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)
  if (any(!is.finite(sm$conf.int[1, c("lower .95", "upper .95")]))) return(NULL)
  tibble(Tumor = first(dat$Tumor), feature = ft, HR = sm$conf.int[1, "exp(coef)"], p = sm$coef[1, "Pr(>|z|)"])
}

onco_cox <- map_dfr(split(oncokb_patients, oncokb_patients$Tumor), function(dat) {
  map_dfr(unique(annotated$feature), ~fit_onco(dat, .x))
}) %>% left_join(feature_level, by = c("Tumor", "feature")) %>% left_join(global_level, by = "feature") %>%
  filter(best_level %in% c("1", "2", "3"), is.finite(HR), HR > 0, is.finite(p), p > 0) %>%
  mutate(level_band = factor(paste0("L", best_level), levels = c("L1", "L2", "L3")),
         match = if_else(Level == best_level, "match", "other"), neglogp = -log10(p), log2HR = -log2(HR))
p4d <- ggplot(onco_cox, aes(feature, Tumor)) +
  geom_point(aes(fill = log2HR, size = neglogp, shape = match, alpha = p <= 0.05), color = "gray30") +
  scale_shape_manual(values = c(match = 21, other = 22)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
  scale_fill_gradient2(low = "coral", mid = "white", high = "blue", midpoint = 0) +
  facet_grid(~level_band, scales = "free_x", space = "free_x") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "-log2(HR)", size = "-log10(P)")
ggsave(file.path(outdir, "Figure4D.pdf"), p4d, width = 13, height = 5)



# Figure 4E: feature-level Molecular Prognostic Index --------------------
# MPI = -log2(Cox HR) * -log10(Cox Wald P) * feature prevalence.
# Repeated subsampling is performed independently for four feature classes.
# Gene mutation features
set.seed(42)

gene_mpi_iterations <- list()
target_cancers_80 <- c("Lung cancer", "Breast cancer", "Intestinal cancer")

for (tumor_type in unique(selected_ctdna$Tumor)) {
  patient_df_master <- selected_ctdna %>% filter(Tumor == tumor_type) %>%
    group_by(patient_id) %>% summarise(os_months = first(os_months), os_event = first(os_event),
                                       mutated_genes = list(unique(gene_symbol)), .groups = "drop")
  
  patient_ids <- patient_df_master$patient_id
  n_total <- length(patient_ids)
  
  genes_to_test <- unique(unlist(patient_df_master$mutated_genes))
  
  for (i in 1:100) {
    sample_ratio <- ifelse(tumor_type %in% target_cancers_80, 0.8, 0.9)
    sampled_patient_ids <- sample(patient_ids, size = floor(n_total * sample_ratio))
    df_patient_sampled <- patient_df_master %>% filter(patient_id %in% sampled_patient_ids)
    
    for (gene in genes_to_test) {
      
      df_analysis <- df_patient_sampled %>%
        mutate(mutated = as.integer(purrr::map_lgl(mutated_genes, ~ gene %in% .x)))
      
      if (sum(df_analysis$mutated) < 2 || length(unique(df_analysis$mutated)) < 2) next
      
      cox <- tryCatch(coxph(Surv(os_months, os_event) ~ mutated, data = df_analysis), error = function(e) NULL)
      if (is.null(cox)) next
      
      hr <- summary(cox)$coefficients["mutated", "exp(coef)"]
      pval <- summary(cox)$coefficients["mutated", "Pr(>|z|)"]
      prop <- mean(df_analysis$mutated == 1)
      mut_count <- sum(df_analysis$mutated == 1)
      sample_count <- nrow(df_analysis)
      score <- (-log2(hr)) * (-log10(pval)) * prop
      
      gene_mpi_iterations[[length(gene_mpi_iterations) + 1]] <- data.frame(
        Tumor = tumor_type, Gene = gene, HR = hr, pval = pval,
        prop = prop, score = score, mut_count = mut_count, sample_count = sample_count,
        label = paste(tumor_type, gene, sep = ":"), iter = i)}}}

gene_mpi_results <- bind_rows(gene_mpi_iterations)

# Mutation count groups
set.seed(42)
target_cancers_80 <- c("Lung cancer", "Breast cancer", "Intestinal cancer")
mutation_group_iterations <- list()

for (tumor_type in unique(selected_ctdna$Tumor)) {
  df_all <- selected_ctdna %>% filter(Tumor == tumor_type) %>% filter(!is.na(mutation_group)) %>% distinct(sample_barcode, .keep_all = TRUE)
  
  if (!("mutation_group" %in% colnames(df_all))) next
  if (length(unique(df_all$mutation_group)) < 2) next
  
  sample_ids_all <- unique(df_all$sample_barcode)
  n_total <- length(sample_ids_all)
  
  for (i in 1:100) {
    sample_ratio <- ifelse(tumor_type %in% target_cancers_80, 0.8, 0.9)
    sampled_ids <- sample(sample_ids_all, size = floor(n_total * sample_ratio))
    
    df <- df_all %>% filter(sample_barcode %in% sampled_ids)
    
    for (target_group in c("mut-H", "mut-L")) {
      df$mutation_group <- factor(df$mutation_group, levels = if (target_group == "mut-H") c("mut-L", "mut-H") else c("mut-H", "mut-L"))
      
      cox <- tryCatch(coxph(Surv(os_months, os_event) ~ mutation_group, data = df), error = function(e) NULL)
      
      if (!is.null(cox)) {
        coef_name <- paste0("mutation_group", target_group)
        if (coef_name %in% rownames(summary(cox)$coefficients)) {
          hr <- summary(cox)$coefficients[coef_name, "exp(coef)"]
          pval <- summary(cox)$coefficients[coef_name, "Pr(>|z|)"]
          prop <- mean(df$mutation_group == target_group)
          mut_count <- sum(df$mutation_group == target_group)
          sample_count <- nrow(df)
          score <- (-log2(hr)) * (-log10(pval)) * prop
          
          mutation_group_iterations[[length(mutation_group_iterations) + 1]] <- data.frame(
            Tumor = tumor_type, Gene = "mutation_group", group_comp = target_group,
            HR = hr, pval = pval, prop = prop, score = score,
            mut_count = mut_count, sample_count = sample_count, label = paste(tumor_type, target_group, sep = ":"), iter = i)}}}}}

mutation_group_mpi_results <- bind_rows(mutation_group_iterations)

# OncoKB evidence-level contrasts
set.seed(42)

oncokb_level_iterations <- list()
level_priority <- c("1", "2", "3", "4", "None")
level_groups <- list(L1_vs_L4N = c("1"), L12_vs_L4N = c("1", "2"), L123_vs_L4N = c("1", "2", "3"))

mpi_onco_cohort <- selected_ctdna %>%
  filter((Tumor == "Lung cancer" & tumor_subtype %in% c("LCLC", "LUAD", "LUSC", "NSCLCPD")) |
           (Tumor == "Intestinal cancer" & tumor_subtype == "COADREAD") |
           Tumor %in% c("Breast cancer") |
           !(Tumor %in% c("Lung cancer", "Intestinal cancer", "Breast cancer")))

for (tumor in unique(mpi_onco_cohort$Tumor)) {
  
  tumor_display_name <- case_when(tumor == "Lung cancer" ~ "NSCLC", tumor == "Intestinal cancer" ~ "CRC", TRUE ~ tumor)
  
  df_tumor_all <- mpi_onco_cohort %>% filter(Tumor == tumor)
  
  onco_sub <- oncokb_annotations %>% filter(category %in% c("all_mutation", "specific_mutation")) %>%
    filter(Tumor == tumor | Tumor == "Pan-cancer") %>% mutate(Level = as.character(Level))
  
  ann_all <- onco_sub %>% filter(category == "all_mutation") %>%
    inner_join(df_tumor_all, by = c("Gene" = "gene_symbol", "Tumor")) %>% dplyr::select(sample_barcode, Level)
  
  ann_spec <- onco_sub %>% filter(category == "specific_mutation") %>%
    inner_join(df_tumor_all, by = c("Gene" = "gene_symbol", "Alterations" = "protein_change", "Tumor")) %>%
    dplyr::select(sample_barcode, Level)
  
  ann_combined <- bind_rows(ann_all, ann_spec) %>% group_by(sample_barcode) %>%
    summarise(Levels = list(Level), .groups = "drop") %>%
    mutate(Level = sapply(Levels, function(x) {
      x <- unique(x)
      x[which.min(match(x, level_priority))]}))
  
  df_master_annotated <- df_tumor_all %>% group_by(sample_barcode) %>%
    summarise(os_months = first(os_months), os_event = first(os_event), .groups = "drop") %>%
    left_join(ann_combined %>% dplyr::select(-Levels), by = "sample_barcode") %>%
    mutate(Level = ifelse(is.na(Level), "None", Level))
  
  sample_ids <- df_master_annotated$sample_barcode # 1 patient 1 sample
  n_total <- length(sample_ids)
  
  for (i in 1:100) {
    sample_ratio <- ifelse((tumor == "Lung cancer") | (tumor == "Intestinal cancer") | (tumor == "Breast cancer"), 0.8, 0.9)
    sampled_ids <- sample(sample_ids, size = floor(n_total * sample_ratio))
    
    df_sampled <- df_master_annotated %>% filter(sample_barcode %in% sampled_ids)
    
    for (group_name in names(level_groups)) {
      level_in <- level_groups[[group_name]]
      
      df_merged <- df_sampled %>%
        mutate(group = case_when(Level %in% level_in ~ "High", Level %in% c("4", "None") ~ "Low", TRUE ~ NA_character_)) %>%
        filter(!is.na(group))
      
      n_high <- sum(df_merged$group == "High")
      n_low <- sum(df_merged$group == "Low")
      n_total_group <- nrow(df_merged)
      
      if (n_high < 2 | n_low < 2) next
      
      cox <- tryCatch(coxph(Surv(os_months, os_event) ~ group, data = df_merged), error = function(e) NULL)
      if (is.null(cox)) next
      
      hr <- summary(cox)$coefficients["groupLow", "exp(coef)"]
      pval <- summary(cox)$coefficients["groupLow", "Pr(>|z|)"]
      prop <- n_high / n_total_group
      score <- (-log2(hr)) * (-log10(pval)) * prop
      
      oncokb_level_iterations[[length(oncokb_level_iterations) + 1]] <- data.frame(
        Tumor = tumor_display_name, Gene = group_name, HR = hr, pval = pval, prop = prop, score = score,
        mut_count = n_high, sample_count = n_total_group, label = paste(tumor_display_name, group_name, sep = ":"), iter = i)}}}

oncokb_level_mpi_results <- bind_rows(oncokb_level_iterations)

# OncoKB specific alterations
set.seed(42)
alteration_mpi_iterations <- list()
level_priority <- c("1", "2", "3", "4", "None")

for (tumor_type in unique(selected_ctdna$Tumor)) {
  
  tumor_display_name <- case_when(tumor_type == "Lung cancer" ~ "NSCLC", tumor_type == "Intestinal cancer" ~ "CRC", TRUE ~ tumor_type)
  
  df_tumor_all <- selected_ctdna %>% filter(Tumor == tumor_type) %>%
    (function(d) {
      if (tumor_type == "Lung cancer") {
        filter(d, tumor_subtype %in% c("LUAD", "NSCLCPD", "LUSC", "LCLC"))
      } else if (tumor_type == "Intestinal cancer") {
        filter(d, tumor_subtype == "COADREAD")
      } else {
        d
      }}) %>%
    group_by(sample_barcode) %>%
    summarise(os_months = first(os_months), os_event = first(os_event),
              genes = list(unique(gene_symbol)), alterations = list(unique(protein_change)), .groups = "drop")
  
  sample_ids_all <- unique(df_tumor_all$sample_barcode)
  n_total <- length(sample_ids_all)
  
  onco_sub <- oncokb_annotations %>% filter(category == "specific_mutation") %>% filter(Tumor == tumor_type | Tumor == "Pan-cancer")
  
  mutation_list <- onco_sub %>% dplyr::select(Gene, Alterations, Level) %>% distinct()
  
  for (iter in 1:100) {
    sample_ratio <- ifelse(tumor_type %in% c("Lung cancer", "Intestinal cancer", "Breast cancer"), 0.8, 0.9)
    sampled_ids <- sample(sample_ids_all, size = floor(n_total * sample_ratio))
    df_tumor <- df_tumor_all %>% filter(sample_barcode %in% sampled_ids)
    
    for (i in seq_len(nrow(mutation_list))) {
      gene <- mutation_list$Gene[i]
      alteration <- mutation_list$Alterations[i]
      level <- mutation_list$Level[i]
      
      if (gene == "EGFR") {
        exon19_dels <- c("E746_A750del", "L747_S752del", "L747_T751del", "L747_P753delinsS")
        exon20_ins <- "A767_V769dup"
        
        if (alteration %in% exon19_dels) {
          alt_label <- "Exon19_In_frame_Del"
          mutated_samples <- df_tumor %>%
            filter(purrr::map2_lgl(genes, alterations, ~ "EGFR" %in% .x && any(exon19_dels %in% .y))) %>%
            pull(sample_barcode)
        } else if (alteration %in% exon20_ins) {
          alt_label <- "Exon20_In_frame_Ins"
          mutated_samples <- df_tumor %>%
            filter(purrr::map2_lgl(genes, alterations, ~ "EGFR" %in% .x && exon20_ins %in% .y)) %>%
            pull(sample_barcode)
        } else {
          alt_label <- alteration
          mutated_samples <- df_tumor %>% filter(purrr::map2_lgl(genes, alterations, ~ gene %in% .x && alteration %in% .y)) %>%
            pull(sample_barcode)
        }}
      
      else if ((gene == "KRAS" || gene == "NRAS") && alteration == "Wildtype" && tumor_type == "Intestinal cancer") {
        alt_label <- "Wildtype"
        gene_mutated <- df_tumor %>% filter(purrr::map_lgl(genes, ~ gene %in% .x)) %>% pull(sample_barcode)
        mutated_samples <- setdiff(df_tumor$sample_barcode, gene_mutated)
      }
      
      else {
        alt_label <- alteration
        mutated_samples <- df_tumor %>% filter(purrr::map2_lgl(genes, alterations, ~ gene %in% .x && alteration %in% .y)) %>%
          pull(sample_barcode)
      }
      
      df_gene <- df_tumor %>% distinct(sample_barcode, os_months, os_event) %>%
        mutate(mutated = ifelse(sample_barcode %in% mutated_samples, 1, 0))
      
      if (sum(df_gene$mutated) < 2 || length(unique(df_gene$mutated)) < 2) next
      
      cox <- tryCatch(coxph(Surv(os_months, os_event) ~ mutated, data = df_gene), error = function(e) NULL)
      if (is.null(cox)) next
      
      hr <- summary(cox)$coefficients["mutated", "exp(coef)"]
      pval <- summary(cox)$coefficients["mutated", "Pr(>|z|)"]
      prop <- mean(df_gene$mutated == 1)
      mut_count <- sum(df_gene$mutated == 1)
      sample_count <- nrow(df_gene)
      score <- (-log2(hr)) * (-log10(pval)) * prop
      label <- paste0(tumor_display_name, ":L", level, "_", gene, "_", alt_label)
      
      alteration_mpi_iterations[[length(alteration_mpi_iterations) + 1]] <- data.frame(
        Tumor = tumor_display_name, Gene = paste0("L", level, "_", gene, "_", alt_label), HR = hr, pval = pval, prop = prop, score = score,
        mut_count = mut_count, sample_count = sample_count, label = label, iter = iter)
    }}}

alteration_mpi_results <- bind_rows(alteration_mpi_iterations)

# Standardize scores within each cancer type and iteration.
mpi <- bind_rows(gene_mpi_results, dplyr::select(mutation_group_mpi_results, -group_comp),
                 oncokb_level_mpi_results, alteration_mpi_results) %>%
  group_by(iter, Tumor) %>%
  mutate(score_z_iter = as.numeric(scale(score)), score_z_iter_abs = abs(score_z_iter)) %>% ungroup()
message("MPI score rows calculated: ", nrow(mpi))


mpi_rank <- mpi %>% group_by(label) %>% summarise(median_score = median(score_z_iter, na.rm = TRUE), .groups = "drop")
show_n <- 25L
selected_mpi <- bind_rows(mpi_rank %>% slice_max(median_score, n = show_n),
                          mpi_rank %>% slice_min(median_score, n = show_n)) %>% distinct(label)
mpi_plot <- mpi %>% semi_join(selected_mpi, by = "label") %>% left_join(mpi_rank, by = "label") %>%
  mutate(Direction = if_else(median_score >= 0, "Favorable", "Unfavorable"))
mpi_colors <- c(cancer_colors, NSCLC = unname(cancer_colors["Lung cancer"]),
                CRC = unname(cancer_colors["Intestinal cancer"]))
p4e <- ggplot(mpi_plot, aes(fct_reorder(label, median_score), score_z_iter, fill = Tumor)) +
  geom_violin(scale = "width", trim = FALSE) + geom_boxplot(width = 0.12, outlier.size = 0.3, fill = "white") +
  facet_wrap(~Direction, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = mpi_colors) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7)) + labs(x = NULL, y = "MPI (z score)")
ggsave(file.path(outdir, "Figure4E.pdf"), p4e, width = 13, height = 9)



# Figure 4F-G: repeated Cox elastic-net ----------------------------------
exclude_genes <- c("ETV1", "ETV4", "ETV5", "FLT3", "GNAQ", "GNAS", "JAK1", "JAK2", "JAK3",
                   "MAP2K1", "MAP2K2", "MAP2K4", "MPL", "NPM1", "RUNX1", "SMARCB1")
level_groups <- list(L1_vs_L4N = c("1"), L12_vs_L4N = c("1", "2"), L123_vs_L4N = c("1", "2", "3"))

make_tumor_matrix <- function(tm) {
  maf_tm <- selected_ctdna %>% filter(Tumor == tm, is.finite(os_months), !is.na(os_event))
  core <- maf_tm %>% distinct(patient_id, os_months, os_event, mutation_group)
  genes <- maf_tm %>% filter(!gene_symbol %in% exclude_genes) %>%
    distinct(patient_id, gene_symbol) %>% mutate(present = 1L) %>%
    pivot_wider(names_from = gene_symbol, values_from = present, values_fill = 0)
  mat <- core %>% left_join(genes, by = "patient_id") %>%
    mutate(high_mutation_count = as.integer(mutation_group == "mut-H")) %>% dplyr::select(-mutation_group)
  has_level <- annotated %>% filter(Tumor == tm) %>% group_by(patient_id) %>%
    summarise(Level = min(as.numeric(Level)), .groups = "drop")
  mat <- mat %>% left_join(has_level, by = "patient_id")
  for (nm in names(level_groups)) {
    mat[[nm]] <- ifelse(mat$Level %in% as.numeric(level_groups[[nm]]), 1L,
                        ifelse(is.na(mat$Level) | mat$Level >= 4, 0L, NA_integer_))
  }
  mat <- mat %>% dplyr::select(-Level)
  alt <- annotated %>% filter(Tumor == tm, str_detect(feature, " ")) %>% distinct(patient_id, feature, Level) %>%
    mutate(feature_id = make.names(paste0("L", Level, "_", feature)), present = 1L) %>%
    dplyr::select(patient_id, feature_id, present) %>% distinct() %>%
    pivot_wider(names_from = feature_id, values_from = present, values_fill = 0)
  if (nrow(alt)) mat <- mat %>% left_join(alt, by = "patient_id")
  mat %>% mutate(across(-c(patient_id, os_months, os_event), ~replace_na(.x, 0))) %>%
    distinct(patient_id, .keep_all = TRUE)
}

boot_coefs <- list()
full_mats <- list()
for (tm in unique(selected_ctdna$Tumor)) {
  dat <- make_tumor_matrix(tm)
  if (nrow(dat) < 10 || sum(dat$os_event == 1, na.rm = TRUE) < 3) next
  full_mats[[tm]] <- dat
  ratio <- if (tm %in% c("Lung cancer", "Breast cancer", "Intestinal cancer")) 0.8 else 0.9
  set.seed(42)
  boot_coefs[[tm]] <- map_dfr(seq_len(100), function(iter) {
    sub <- dat %>% slice_sample(n = floor(nrow(dat) * ratio), replace = FALSE)
    feat <- setdiff(names(sub), c("patient_id", "os_months", "os_event"))
    keep <- feat[vapply(sub[feat], function(x) sum(x == 1, na.rm = TRUE) >= 2 &&
                          sum(x == 0, na.rm = TRUE) >= 2, logical(1))]
    if (!length(keep)) return(NULL)
    x <- data.matrix(sub[keep]); y <- Surv(sub$os_months, sub$os_event)
    cv <- tryCatch(glmnet::cv.glmnet(x, y, family = "cox", alpha = 0.5), error = function(e) NULL)
    if (is.null(cv)) return(NULL)
    cf <- as.matrix(coef(cv, s = "lambda.min")); cf <- cf[cf[, 1] != 0, , drop = FALSE]
    if (!nrow(cf)) return(NULL)
    tibble(Tumor = tm, feature = rownames(cf), coefficient = as.numeric(cf[, 1]), iter = iter)
  })
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
    z <- dat %>% transmute(os_months, os_event, marker = .data[[f]])
    fit <- tryCatch(coxph(Surv(os_months, os_event) ~ marker, data = z), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    tibble(Tumor = tm, feature = f, p = summary(fit)$coef[1, "Pr(>|z|)"])
  })
})
plot_enet <- enet_summary %>% left_join(feature_p, by = c("Tumor", "feature")) %>%
  mutate(neglogp = -log10(pmax(p, .Machine$double.xmin)))
p4f <- ggplot(plot_enet, aes(mean_weight, times)) +
  geom_point(aes(size = neglogp, fill = if_else(times >= 50, Tumor, "Other")), shape = 21, alpha = 0.85) +
  geom_text_repel(data = subset(plot_enet, times >= 50), aes(label = paste(Tumor, feature, sep = ": ")),
                  size = 2.5, max.overlaps = 70) + geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = 50, linetype = "dashed") +
  scale_fill_manual(values = c(cancer_colors, Other = "gray85")) + theme_classic() +
  labs(x = "Averaged weights", y = "Bootstrap frequency", size = "-log10(P)")
ggsave(file.path(outdir, "Figure4F.pdf"), p4f, width = 9, height = 7)


lc <- enet_summary %>% filter(Tumor == "Lung cancer", times >= 50) %>% arrange(mean_weight)
if (nrow(lc)) {
  lc <- lc %>% mutate(feature = factor(feature, levels = feature))
  p4g_coef <- ggplot(lc, aes(mean_weight, feature, fill = coef_mean < 0)) +
    geom_col() + scale_fill_manual(values = c(`TRUE` = "#FF9933", `FALSE` = "#8A3700"), guide = "none") +
    theme_classic() + labs(x = "Mean elastic-net coefficient (sign reversed)", y = NULL)
  lc_mat <- full_mats[["Lung cancer"]] %>% arrange(os_months) %>% mutate(patient_order = row_number())
  p4g_heat <- lc_mat %>% dplyr::select(patient_order, all_of(as.character(lc$feature))) %>%
    pivot_longer(-patient_order, names_to = "feature", values_to = "present") %>%
    ggplot(aes(patient_order, feature, fill = factor(present))) + geom_tile() +
    scale_fill_manual(values = c(`0` = "white", `1` = "black")) +
    theme_classic() + theme(axis.text.x = element_blank(), legend.position = "none")
  p4g_os <- ggplot(lc_mat, aes(patient_order, os_months)) + geom_line() +
    theme_classic() + labs(x = "Patient", y = "OS (months)")
  ggsave(file.path(outdir, "Figure4G.pdf"), p4g_coef | (p4g_heat / p4g_os), width = 13, height = 6)
}



# Figure 4H-I: tumor-stratified split, Cox models, and SHAP ---------------
model_patients <- ctdna_data %>% group_by(patient_id, Tumor) %>%
  summarise(os_months = first_valid(os_months), os_event = first_valid(os_event),
            age_group = first_valid(age_group), sex = first_valid(sex),
            clinical_stage = first_valid(clinical_stage), mutation_count = max(mutation_count, na.rm = TRUE),
            mutation_group = if_else(any(mutation_group == "mut-H", na.rm = TRUE), "mut-H", "mut-L"), .groups = "drop") %>%
  mutate(stage_group = case_when(clinical_stage %in% c("1", "2") ~ "Low",
                                 clinical_stage %in% c("3", "4") ~ "High", TRUE ~ NA_character_),
         stage_group = factor(stage_group, levels = c("High", "Low")), sex = factor(sex, levels = c("M", "F")),
         age_group = factor(age_group, levels = c("Older", "Younger")),
         mutation_group = factor(mutation_group, levels = c("mut-H", "mut-L"))) %>%
  filter(is.finite(os_months), os_months > 0, !is.na(os_event),
         !is.na(stage_group), !is.na(age_group), !is.na(sex), !is.na(mutation_group)) %>% arrange(patient_id)

split_by_tumor <- function(dat, test_prop = 0.2, seed = 2026031922) {
  set.seed(seed)
  splits <- split(dat, dat$Tumor)
  test <- bind_rows(lapply(splits, function(d) {
    n <- nrow(d); if (n <= 1) return(d[FALSE, ])
    n_test <- min(max(1L, round(n * test_prop)), n - 1L)
    idx <- sample.int(n, size = n, replace = FALSE)
    d[idx[seq_len(n_test)], , drop = FALSE]
  }))
  list(train = anti_join(dat, test, by = "patient_id"), test = test)
}
make_folds <- function(dat, k = 5, seed = 2026031922) {
  set.seed(seed)
  bind_rows(lapply(split(dat, dat$Tumor), function(d) {
    d$fold_id <- sample(rep(seq_len(k), length.out = nrow(d)))
    d
  }))
}

split <- split_by_tumor(model_patients)
train_folds <- make_folds(split$train)
fold1_training <- train_folds %>% filter(fold_id != 1) %>% dplyr::select(-fold_id)
test_patients <- split$test


# Save the survival curve and risk table together in a single PDF.
save_survival_plot <- function(plot_obj, filename) {
  arranged <- survminer::arrange_ggsurvplots(list(plot_obj), print = FALSE,
                                             ncol = 1, nrow = 1, risk.table.height = 0.25)
  ggplot2::ggsave(file.path(outdir, filename), plot = arranged,
                  width = 7, height = 6, units = "in", device = "pdf")
}

# Supplementary figure: Model 1 (clinical variables only).
fit1 <- coxph(Surv(os_months, os_event) ~ age_group + sex + stage_group + mutation_group + strata(Tumor),
              data = fold1_training, x = TRUE, y = TRUE)
test1 <- test_patients %>% mutate(risk_score = as.numeric(predict(fit1, newdata = test_patients, type = "lp"))) %>%
  mutate(risk_group = factor(if_else(risk_score > median(risk_score, na.rm = TRUE), "High", "Low"),
                             levels = c("Low", "High")))
fit_supp_model1 <- survfit(Surv(os_months, os_event) ~ risk_group, data = test1)
p_supp_model1 <- ggsurvplot(fit_supp_model1, data = test1, pval = TRUE, risk.table = TRUE,
                            conf.int = FALSE, surv.median.line = "hv", legend.title = "Risk group",
                            legend.labs = c("Low", "High"), ggtheme = theme_bw())
save_survival_plot(p_supp_model1, "Supplementary_Model1_test_survival.pdf")

# Model 3: clinical variables + patient-level gene mutation indicators.
patient_gene_features <- ctdna_data %>% distinct(patient_id, gene_symbol) %>% mutate(present = 1L)
gene_freq <- patient_gene_features %>% semi_join(model_patients, by = "patient_id") %>%
  dplyr::count(gene_symbol, name = "n") %>% mutate(prevalence = n / nrow(model_patients))
selected_genes <- gene_freq %>% filter(prevalence >= 0.01) %>% pull(gene_symbol)
gene_feature_matrix <- patient_gene_features %>% filter(gene_symbol %in% selected_genes) %>%
  pivot_wider(names_from = gene_symbol, values_from = present, values_fill = 0)
names(gene_feature_matrix)[-1] <- make.names(names(gene_feature_matrix)[-1], unique = TRUE)
model3_data <- model_patients %>% left_join(gene_feature_matrix, by = "patient_id") %>%
  mutate(across(all_of(make.names(selected_genes)), ~replace_na(.x, 0L)))
model3_training <- model3_data %>% semi_join(fold1_training, by = "patient_id")
model3_validation <- model3_data %>% semi_join(train_folds %>% filter(fold_id == 1), by = "patient_id")
model3_test <- model3_data %>% semi_join(test_patients, by = "patient_id")
formula3 <- reformulate(c("age_group", "sex", "stage_group", "mutation_group", make.names(selected_genes),
                          "strata(Tumor)"), response = "Surv(os_months, os_event)")
fit3 <- coxph(formula3, data = model3_training, x = TRUE, y = TRUE)

# Main Figure 4H: survival curves using Model 3 risk scores in the test set.
model3_test <- model3_test %>% mutate(risk_score = as.numeric(predict(fit3, newdata = model3_test, type = "lp"))) %>%
  mutate(risk_group = factor(if_else(risk_score > median(risk_score, na.rm = TRUE), "High", "Low"),
                             levels = c("Low", "High")))
fit4h <- survfit(Surv(os_months, os_event) ~ risk_group, data = model3_test)
p4h <- ggsurvplot(fit4h, data = model3_test, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                  surv.median.line = "hv", legend.title = "Risk group",
                  legend.labs = c("Low", "High"), ggtheme = theme_bw())
save_survival_plot(p4h, "Figure4H.pdf")

# Figure 4I: Kernel SHAP for Model 3 in the fold-1 validation cohort.
x_explain <- as.data.frame(model.matrix(fit3, data = model3_validation))
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

