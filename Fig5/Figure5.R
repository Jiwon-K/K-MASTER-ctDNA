
# Figure 5: ctDNA features and outcomes after first-line platinum chemotherapy ----
# Figure-level analysis using controlled-access data and the curated treatment cohort.
# The treatment cohort contains the first-line platinum treatment date, follow-up
# date or next-treatment date, and the outcome fields used in the original analysis.

# Packages ----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(forcats)
library(purrr)
library(survival)
library(survminer)
library(ggrepel)
library(patchwork)


# Input files --------------------------------------------------------------
mutation_file <- read.delim("KM00_AXEN.maf.gz", sep = "\t", check.names = F) 
metadata_file <- "KM00_clinical.xlsx"
treatment_file <- "KM00_TTNT.txt"
output_dir <- "results/Figure5"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(mutation_file, metadata_file, treatment_file)
if (any(!file.exists(input_files))) stop("Missing input files: ", paste(input_files[!file.exists(input_files)], collapse = ", "))


# Data preparation ---------------------------------------------------------
mutations <- readRDS(mutation_file)
required_mutation_cols <- c("SUBJECT ID", "SAMPLE ID", "Tumor_Sample_Barcode", "Tumor",
                            "Hugo_Symbol", "count", "Survival_250417")
if (!all(required_mutation_cols %in% names(mutations))) {
  stop("Required mutation columns are missing: ", paste(setdiff(required_mutation_cols, names(mutations)), collapse = ", "))
}
mutations <- mutations %>%
  rename(patient_id = `SUBJECT ID`, sample_id = `SAMPLE ID`, sample_barcode = Tumor_Sample_Barcode,
         tumor_type = Tumor, gene = Hugo_Symbol, mutation_count = count, os_event = Survival_250417) %>%
  mutate(os_event = as.numeric(os_event), mutation_count = as.numeric(mutation_count),
         mutation_group = case_when(mutation_count >= 3 ~ "mut-H", mutation_count < 3 ~ "mut-L", TRUE ~ NA_character_))

sample_metadata <- read_xlsx(metadata_file, sheet = "KM00_NGS_FINAL_METADATA_FINALan")
stopifnot(all(c("SAMPLE ID", "FINAL_new") %in% names(sample_metadata)))
selected_ids <- sample_metadata %>% filter(FINAL_new == "O") %>% pull(`SAMPLE ID`) %>% unique()
mutations <- mutations %>% filter(sample_id %in% selected_ids)


treatment <- read.delim(treatment_file, sep = "\t", fileEncoding = "euc-kr", check.names = FALSE)
required_treatment_cols <- c("SUBJECT ID", "Tumor", "TTNT_months", "status_edit")
if (!all(required_treatment_cols %in% names(treatment))) {
  stop("Required treatment columns are missing: ", paste(setdiff(required_treatment_cols, names(treatment)), collapse = ", "))
}
treatment <- treatment %>%
  rename(patient_id = `SUBJECT ID`, tumor_type = Tumor, ttnt_months = TTNT_months,
         followup_event = status_edit) %>%
  mutate(ttnt_months = pmax(as.numeric(ttnt_months), 0),
         followup_event = suppressWarnings(as.numeric(followup_event))) %>%
  distinct(patient_id, .keep_all = TRUE)

if (nrow(treatment) != 398L) warning("Figure 5A cohort differs from the original 398 patients: ", nrow(treatment))

major_tumors <- c("Bladder cancer", "Breast cancer", "Gallbladder cancer",
                  "Gastric cancer", "Head and neck cancer", "Lung cancer")
tumor_labels <- c("Bladder cancer" = "BLCA", "Breast cancer" = "BC",
                  "Gallbladder cancer" = "BTC", "Gastric cancer" = "GC",
                  "Head and neck cancer" = "HNSC", "Lung cancer" = "LC")
tumor_colors <- c("Lung cancer" = "#D55E00", "Breast cancer" = "#0072B2",
                  "Gallbladder cancer" = "#F0E442", "Gastric cancer" = "#E69F00",
                  "Head and neck cancer" = "#D4A10A", "Bladder cancer" = "#E1C16E",
                  "Esophageal cancer" = "#F1B1B1", "Neuroendocrine tumor" = "#AD6F6F",
                  "Ovarian cancer" = "#AABF82", "Prostate cancer" = "#B97B60",
                  "Liver cancer" = "#9C3B3B", "Intestinal cancer" = "#009E73",
                  "Renal cell cancer" = "#5582B0", "Sarcoma" = "#6C8E3E",
                  "Pancreatic cancer" = "#56B4E9", "Thyroid cancer" = "#7F5A6D")

cohort_mutations <- mutations %>% filter(tumor_type %in% major_tumors) %>%
  semi_join(treatment, by = c("patient_id", "tumor_type")) %>%
  left_join(treatment %>% dplyr::select(patient_id, ttnt_months, followup_event), by = "patient_id")
cohort_patients <- cohort_mutations %>% distinct(patient_id, .keep_all = TRUE)
if (nrow(cohort_patients) != n_distinct(cohort_patients$sample_barcode)) {
  stop("The selected cohort includes more than one sample per patient. Check the sample selection metadata.")
}
message("First-line platinum patients: ", nrow(treatment), "; six-tumor ctDNA cohort: ", nrow(cohort_patients))
if (nrow(cohort_patients) != 359L) warning("Original six-tumor ctDNA cohort had 359 patients; observed ", nrow(cohort_patients))



# Figure 5A: tumor distribution -------------------------------------------
tumor_dist <- treatment %>% dplyr::count(tumor_type, name = "patients") %>%
  mutate(plot_group = if_else(patients < 10, "Other", tumor_type)) %>%
  group_by(plot_group) %>% summarise(patients = sum(patients), .groups = "drop") %>%
  arrange(desc(patients)) %>%
  mutate(fraction = patients / sum(patients), plot_group = factor(plot_group, levels = plot_group),
         label = paste0(round(100 * fraction, 1), "%\n(n=", patients, ")"))

p5a <- ggplot(tumor_dist, aes(x = 2, y = fraction, fill = plot_group)) +
  geom_col(width = 1, color = "white") + coord_polar(theta = "y", direction = -1) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  xlim(0.5, 2.5) + theme_void() +
  scale_fill_manual(values = c(tumor_colors, "Other" = "grey70")) +
  annotate("text", x = 0.5, y = 0, label = paste0("First-line platinum\n(n = ", nrow(treatment), ")"),
           fontface = "bold", size = 3.5) + labs(fill = "Tumor type")
ggsave(file.path(output_dir, "Figure5A.pdf"), p5a, width = 8, height = 6)



# Figure 5B: TTNT across six tumor types ---------------------------------
comparison <- treatment %>% filter(tumor_type %in% major_tumors, is.finite(ttnt_months)) %>%
  mutate(tumor_type = factor(tumor_type, levels = major_tumors))
kw_p <- kruskal.test(ttnt_months ~ tumor_type, data = comparison)$p.value
message("Figure 5B Kruskal-Wallis P: ", signif(kw_p, 4))
p5b <- ggplot(comparison, aes(tumor_type, ttnt_months, fill = tumor_type)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.14, fill = "white", outlier.size = 0.6) +
  scale_x_discrete(labels = tumor_labels) + scale_fill_manual(values = tumor_colors) +
  annotate("text", x = 4.4, y = max(comparison$ttnt_months, na.rm = TRUE),
           label = paste0("P = ", signif(kw_p, 3), " (Kruskal-Wallis)"), size = 3) +
  theme_classic() + theme(legend.position = "none") + labs(x = NULL, y = "Time-to-next-treatment (months)")
ggsave(file.path(output_dir, "Figure5B.pdf"), p5b, width = 8, height = 5)

# Survival plot helper: retain survival curves and risk tables in valid PDFs.
save_survival <- function(data, formula, filename, palette = NULL, legend_title = NULL, size = c(7, 6)) {
  fit <- survival::survfit(formula, data = data)
  lr <- survival::survdiff(formula, data = data)
  logrank_p <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  message(filename, ": log-rank P = ", signif(logrank_p, 4))
  chart <- survminer::ggsurvplot(fit, data = data, pval = TRUE, risk.table = TRUE,
                                 conf.int = FALSE, palette = palette, legend.title = legend_title,
                                 xlab = "Months", ylab = "Time-to-next-treatment probability",
                                 ggtheme = theme_classic())
  arranged <- survminer::arrange_ggsurvplots(list(chart), print = FALSE, ncol = 1, nrow = 1)
  ggplot2::ggsave(file.path(output_dir, filename), plot = arranged,
                  width = size[1], height = size[2], units = "in", device = "pdf")
  invisible(fit)
}



# Figure 5C: mutation count group and TTNT -------------------------------
count_survival <- cohort_patients %>%
  filter(is.finite(ttnt_months), !is.na(followup_event), mutation_group %in% c("mut-H", "mut-L")) %>%
  mutate(mutation_group = factor(mutation_group, levels = c("mut-H", "mut-L")))
save_survival(count_survival, Surv(ttnt_months, followup_event) ~ mutation_group,
              "Figure5C.pdf", palette = c("#F5B800", "#0072B2"), legend_title = "Mutation group")


# Common Cox helper ------------------------------------------------------
cox_binary <- function(data, feature, time_col = "ttnt_months", event_col = "followup_event") {
  present <- as.integer(feature)
  dat <- data.frame(time = as.numeric(data[[time_col]]), event = as.numeric(data[[event_col]]), present = present)
  dat <- dat[is.finite(dat$time) & !is.na(dat$event) & !is.na(dat$present), , drop = FALSE]
  if (nrow(dat) < 4 || length(unique(dat$present)) < 2 || min(table(dat$present)) < 2 ||
      any(tapply(dat$event, dat$present, sum) < 1)) return(NULL)
  fit <- tryCatch(survival::coxph(Surv(time, event) ~ present, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)
  hr <- unname(sm$conf.int[1, "exp(coef)"])
  p <- unname(sm$coefficients[1, "Pr(>|z|)"])
  if (!is.finite(hr) || !is.finite(p) || hr <= 0 || p <= 0) return(NULL)
  tibble(hr = hr, p = p, n_altered = sum(dat$present == 1), n_total = nrow(dat))
}



# Figure 5D: gene-level TTNT association ---------------------------------
# Original gene-level panel used the available OS-status field as the event.
patient_outcomes <- cohort_patients %>% dplyr::select(patient_id, sample_barcode, tumor_type, ttnt_months, os_event)
mutation_calls <- cohort_mutations %>% filter(!is.na(gene)) %>% distinct(sample_barcode, tumor_type, gene)
gene_results <- map_dfr(major_tumors, function(tm) {
  patients_tm <- patient_outcomes %>% filter(tumor_type == tm)
  mutations_tm <- mutation_calls %>% filter(tumor_type == tm)
  map_dfr(unique(mutations_tm$gene), function(g) {
    positive <- mutations_tm %>% filter(gene == g) %>% pull(sample_barcode)
    test <- cox_binary(patients_tm, patients_tm$sample_barcode %in% positive, event_col = "os_event")
    if (is.null(test)) return(NULL)
    mutate(test, tumor_type = tm, gene = g)
  })
})
if (!nrow(gene_results)) stop("Figure 5D: no gene-level Cox results. Verify the cohort and event fields.")
gene_order <- gene_results %>% group_by(gene) %>% summarise(n = sum(n_altered), .groups = "drop") %>%
  arrange(desc(n)) %>% pull(gene)
heatmap_data <- gene_results %>%
  tidyr::complete(tumor_type = major_tumors, gene = gene_order) %>%
  mutate(tumor_type = factor(tumor_type, levels = rev(major_tumors)),
         gene = factor(gene, levels = gene_order),
         log2hr = pmax(-2, pmin(2, log2(hr))),
         mark = case_when(p < 0.005 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""))
p5d <- ggplot(heatmap_data, aes(gene, tumor_type, fill = log2hr)) +
  geom_tile(color = "grey92", linewidth = 0.2) + geom_text(aes(label = mark), size = 2.3) +
  scale_y_discrete(labels = tumor_labels) +
  scale_fill_gradient2(low = "#33664D", mid = "white", high = "#F29130", midpoint = 0,
                       limits = c(-2, 2), na.value = "grey96", name = "log2(HR)") +
  theme_minimal() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, size = 7)) +
  labs(x = NULL, y = NULL)
ggsave(file.path(output_dir, "Figure5D.pdf"), p5d, width = 15, height = 5)



# Figure 5E: ERBB2 status in GC and BTC ---------------------------------
erbb2_patients <- cohort_patients %>% filter(tumor_type %in% c("Gastric cancer", "Gallbladder cancer"),
                                             is.finite(ttnt_months), !is.na(followup_event))
erbb2_calls <- mutation_calls %>% filter(gene == "ERBB2") %>% pull(sample_barcode)
erbb2_patients <- erbb2_patients %>% mutate(erbb2_status = factor(if_else(sample_barcode %in% erbb2_calls,
                                                                          "ERBB2 Mut", "ERBB2 WT"),
                                                                  levels = c("ERBB2 Mut", "ERBB2 WT")))
erbb2_panels <- list()
for (tm in c("Gastric cancer", "Gallbladder cancer")) {
  dat <- erbb2_patients %>% filter(tumor_type == tm)
  if (n_distinct(dat$erbb2_status) < 2) stop("Figure 5E: ERBB2 comparison unavailable in ", tm)
  suffix <- if (tm == "Gastric cancer") "GC" else "BTC"
  save_survival(dat, Surv(ttnt_months, followup_event) ~ erbb2_status,
                paste0("Figure5E_", suffix, ".pdf"), palette = c("red", "grey65"), legend_title = NULL)
  erbb2_panels[[suffix]] <- ggsurvplot(survfit(Surv(ttnt_months, followup_event) ~ erbb2_status, data = dat),
                                       data = dat, pval = TRUE, risk.table = FALSE,
                                       palette = c("red", "grey65"), ggtheme = theme_classic())$plot +
    labs(title = suffix, x = "Months", y = "TTNT probability")
}
ggsave(file.path(output_dir, "Figure5E.pdf"), erbb2_panels[["GC"]] | erbb2_panels[["BTC"]],
       width = 11, height = 5)



# Figure 5F-H: oncogenic pathway analysis --------------------------------
# Pathway definitions are taken from the original figure analysis.
pathway_genes <- list(
  "Cell Cycle" = c("ATM", "TP53", "CDKN2A", "MDM2", "CCND1", "CDKN1A", "RB1", "CCNE1", "FBXW7", "E2F3"),
  "DDR" = c("ATM", "MLH1", "BRCA1", "BRCA2", "FANCA", "FANCC", "FANCF", "FANCG"),
  "FA" = c("FANCA", "FANCC", "BRCA1", "BRCA2", "FANCF", "FANCG"),
  "HR" = c("BRCA1", "BRCA2", "RAD51", "RAD51B", "RAD51C", "RAD51D", "PALB2", 
           "BRIP1", "BARD1", "ATM", "ATR", "CHEK1", "CHEK2", "FANCD2", "FANCI"),
  "JAK-STAT" = c("JAK1", "JAK2", "JAK3"),
  "MYC" = c("MYC", "MYCN", "MYCL", "MAX", "MGA", "MXD1", "MXD3", "MXD4", "MXI1", "MNT", "MLX", "MLXIP", "MLXIPL"),
  "NOTCH" = c("NOV", "CNTN6", "JAG2", "ARRDC1", "FBXW7", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", 
              "DNER", "PSEN2", "CUL1", "NCOR1", "NCOR2", "SPEN", "KDM5A", "MALM3", "KAT2B", "CREBBP", "EP300", "NRARP"),
  "PI3K" = c("PTEN", "INPP4B", "PIK3CA", "PIK3R1", "PIK3R2", "PIK3R3", 
             "AKT1", "AKT2", "AKT3", "STK11", "PPP2R1A", "TSC1", "TSC2", "RHEB", "RICTOR", "MTOR", "RPTOR"),
  "RTK-RAS" = c("EGFR", "ERBB2", "ERBB3", "ERBB4", "MET", "PDGFRA", "FGFR1", "FGFR2", "FGFR3", "FGFR4", 
                "KIT", "IGF1R", "RET", "ROS1", "ALK", "FLT3", "NTRK1", "NTRK2", "CBL", "ERRFI1", "PTPN11", "SOS1", "NF1", 
                "RASA1", "KRAS", "HRAS", "NRAS", "RIT1", "ARAF", "BRAF", "RAF1", "RAC1", "MAPK1", "MAP2K1", "MAP2K2"),
  "TP53" = c("MDM2", "MDM4", "CDKN2A", "TP53", "ATM", "CHEK2", "PRS6KA3"),
  "WNT/b-catenin" = c("CTNNB1", "APC", "FBXW7", "RUNX1"))

pathway_status <- patient_outcomes %>% distinct(sample_barcode, .keep_all = TRUE) %>%
  left_join(cohort_patients %>% dplyr::select(sample_barcode, followup_event), by = "sample_barcode")
for (pathway in names(pathway_genes)) {
  carriers <- mutation_calls %>% filter(gene %in% pathway_genes[[pathway]]) %>% pull(sample_barcode)
  pathway_status[[pathway]] <- as.integer(pathway_status$sample_barcode %in% carriers)
}



# Figure 5F: pathway-specific Cox models by tumor type.
pathway_results <- map_dfr(major_tumors, function(tm) {
  dat <- pathway_status %>% filter(tumor_type == tm)
  map_dfr(names(pathway_genes), function(pathway) {
    fit <- cox_binary(dat, dat[[pathway]])
    if (is.null(fit)) return(NULL)
    mutate(fit, tumor_type = tm, pathway = pathway)
  })
})
if (!nrow(pathway_results)) stop("Figure 5F: no pathway-level Cox results")
p5f <- pathway_results %>%
  mutate(log2hr = pmax(-2, pmin(2, log2(hr))), neglogp = pmin(5, -log10(p)),
         significance = if_else(p < 0.05, "*", ""),
         tumor_type = factor(tumor_type, levels = major_tumors),
         pathway = factor(pathway, levels = names(pathway_genes))) %>%
  ggplot(aes(pathway, tumor_type)) +
  geom_point(aes(fill = log2hr, size = neglogp), shape = 21, color = "grey75") +
  geom_text(aes(label = significance), size = 3) +
  scale_y_discrete(labels = tumor_labels) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", limits = c(-2, 2),
                       name = "log2(HR)") + scale_size(range = c(1, 7), name = "-log10(P)") +
  theme_minimal() + theme(panel.grid = element_line(color = "grey92"), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL)
ggsave(file.path(output_dir, "Figure5F.pdf"), p5f, width = 11, height = 5)



# Figure 5G: SELECT co-alteration scores and Cox associations.
# Redundant pathway definitions were excluded from this comparison in the original analysis.
select_pathways <- c("Cell Cycle", "DDR", "JAK-STAT", "MYC", "NOTCH", "PI3K", "RTK-RAS", "WNT/b-catenin")
select_matrix <- pathway_status %>% dplyr::select(sample_barcode, all_of(select_pathways)) %>% distinct() %>%
  tibble::column_to_rownames("sample_barcode") %>% as.matrix()
mode(select_matrix) <- "numeric"
sample_classes <- pathway_status$tumor_type[match(rownames(select_matrix), pathway_status$sample_barcode)]
names(sample_classes) <- rownames(select_matrix)
alteration_classes <- setNames(rep("pathway", length(select_pathways)), select_pathways)

if (!requireNamespace("select", quietly = TRUE)) stop("The SELECT R package used in the original analysis is required for Figure 5G")
dir.create(file.path(output_dir, "SELECT_pathway"), recursive = TRUE, showWarnings = FALSE)
select_results <- select::select(M = select_matrix, sample.class = sample_classes,
                                 alteration.class = alteration_classes,
                                 folder = file.path(output_dir, "SELECT_pathway"), n.cores = 4,
                                 min.feature.support = 2, min.feature.freq = 0.001,
                                 calculate_APC_threshold = FALSE, calculate_FDR = TRUE, verbose = TRUE)
select_scores <- select_results %>%
  mutate(select_score = case_when(direction == "CO" ~ abs(MI_diff), direction == "ME" ~ -abs(MI_diff), TRUE ~ 0),
         pathway_pair = paste(pmin(SFE_1, SFE_2), pmax(SFE_1, SFE_2), sep = " & "))

pair_names <- combn(select_pathways, 2, simplify = FALSE)
pair_cox <- map_dfr(pair_names, function(pair) {
  dat <- pathway_status %>% filter(is.finite(ttnt_months), !is.na(followup_event))
  both <- dat[[pair[1]]] == 1 & dat[[pair[2]]] == 1
  test <- cox_binary(dat, both)
  if (is.null(test)) return(NULL)
  test %>% mutate(pathway_pair = paste(sort(pair), collapse = " & "))
})

comparison_scores <- select_scores %>% left_join(pair_cox, by = "pathway_pair") %>%
  filter(is.finite(select_score), is.finite(hr), is.finite(p), p > 0) %>%
  mutate(x_score = -sign(select_score) * log2(abs(select_score) + 1e-4),
         log2hr = log2(hr), neglogp = -log10(p),
         select_significant = if_else(is.na(FDR), FALSE, as.logical(FDR)))

p5g <- ggplot(comparison_scores, aes(x_score, log2hr)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(aes(size = neglogp, color = select_significant), alpha = 0.85) +
  geom_text_repel(aes(label = pathway_pair), size = 2.7, max.overlaps = 30) +
  scale_color_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey65"),
                     labels = c(`FALSE` = "n.s.", `TRUE` = "Significant"), name = "SELECT") +
  scale_size(range = c(1.5, 5), name = "-log10(P)") +
  theme_classic() + labs(x = "log2(SELECT score)", y = "log2(HR)")
ggsave(file.path(output_dir, "Figure5G.pdf"), p5g, width = 9, height = 7)



# Figure 5H: NOTCH and WNT/b-catenin combined alterations.
pair_survival <- pathway_status %>% filter(is.finite(ttnt_months), !is.na(followup_event)) %>%
  mutate(pathway_group = case_when(`NOTCH` == 1 & `WNT/b-catenin` == 1 ~ "Both altered",
                                   `NOTCH` == 1 | `WNT/b-catenin` == 1 ~ "Either altered",
                                   TRUE ~ "Neither altered"),
         pathway_group = factor(pathway_group, levels = c("Both altered", "Either altered", "Neither altered")))
if (nrow(pair_survival) != 353L) warning("Original Figure 5H analyzed 353 patients; observed ", nrow(pair_survival))
message("Figure 5H group sizes: ", paste(names(table(pair_survival$pathway_group)),
                                         as.integer(table(pair_survival$pathway_group)), collapse = "; "))
save_survival(pair_survival, Surv(ttnt_months, followup_event) ~ pathway_group,
              "Figure5H.pdf", palette = c("#D73027", "#377EB8", "grey50"), legend_title = NULL)

message("Figure 5A-H PDF files saved to: ", output_dir)

