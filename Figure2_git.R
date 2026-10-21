
# Figure 2: K-MASTER versus MSK-ACCESS ctDNA profiles
# Code accompanying: Pan-cancer circulating tumor DNA atlas establishes a
# precision oncology framework for prognostic and therapeutic stratification
# Controlled-access K-MASTER data and study-specific input tables are not included.
# Run from the repository root after placing the required input files locally.


# Packages ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(ComplexHeatmap)
library(grid)


# Input files ------------------------------------------------------------
# Filtered, annotated K-MASTER AXEN ctDNA MAF.
# Relevant columns: SUBJECT ID, Tumor, Tumor_Sample_Barcode, Hugo_Symbol.
maf_file <- "KM00_AXEN.maf.gz"

# K-MASTER clinical annotations including Tumor_subtype_JW.
clinical_file <- "KM00_CLINICAL.xlsx"

# MSK-ACCESS mutation data after the source filtering steps (normal-read count,
# tumor depth, coding-variant class) and merged clinical annotations.
# Relevant columns include Cancer Type, Cancer Type Detailed, Race and barcode.
msk_file <- "msk_ctdna_vte_2024/data_mutations_filtered.txt"

# Study-curated comparison table for the final OncoKB v5.0 frequency matrix.
# Derived upstream from OncoKB annotation matching in the original analysis.
# Expected columns: Gene, Tumor, AXEN, MSK, Level_edit, p_value, label.
oncokb_table <- "onco1-4_axen+msk_freq_nsclc_crc_v5.0_white.txt"

outdir <- "results/Figure2"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# Cohort preparation -----------------------------------------------------
axen <- read.delim(maf_file, sep = "\t", check.names = FALSE) %>% filter(!Tumor %in% c("Double cancer", "Other"))

selected_tumors <- axen %>% distinct(Tumor, `SUBJECT ID`) %>%
  dplyr::count(Tumor, name = "n_patients") %>% filter(n_patients >= 10) %>% pull(Tumor)
axen <- axen %>% filter(Tumor %in% selected_tumors)
stopifnot(n_distinct(axen$Tumor) == 17)

cli <- readxl::read_xlsx(clinical_file, skip = 1)
if (!"SUBJECT ID" %in% names(cli) && "Subject ID" %in% names(cli)) {
  cli <- cli %>% rename(`SUBJECT ID` = `Subject ID`)
}
axen <- axen %>% dplyr::select(-any_of("Tumor_subtype_JW")) %>%
  left_join(cli %>% dplyr::select(`SUBJECT ID`, Tumor_subtype_JW) %>% distinct(`SUBJECT ID`, .keep_all = TRUE), by = "SUBJECT ID") %>%
  mutate(Tumor = sub("Neuroendocine tumor", "Neuroendocrine tumor", Tumor))

# The final comparisons restrict lung and intestinal cancers to NSCLC and CRC.
axen <- axen %>% filter((Tumor == "Lung cancer" & Tumor_subtype_JW %in% c("LCLC", "LUAD", "LUSC", "NSCLCPD")) |
                          (Tumor == "Intestinal cancer" & Tumor_subtype_JW == "COADREAD") |
                          !Tumor %in% c("Lung cancer", "Intestinal cancer"))

msk <- read.delim(msk_file, sep = "\t", check.names = FALSE)
stopifnot(all(c("Race", "Cancer Type", "Cancer Type Detailed", "Tumor_Sample_Barcode") %in% names(msk)))
msk <- msk %>% filter(toupper(trimws(as.character(Race))) == "WHITE") %>%
  mutate(Tumor = case_when(
    `Cancer Type` %in% c("Anal Cancer", "Appendiceal Cancer", "Colorectal Cancer", "Small Bowel Cancer") ~ "Intestinal cancer",
    `Cancer Type` == "Bladder Cancer" ~ "Bladder cancer",
    `Cancer Type` == "Breast Cancer" ~ "Breast cancer",
    `Cancer Type` == "Gastrointestinal Neuroendocrine Tumor" ~ "Neuroendocrine tumor",
    `Cancer Type` %in% c("Gastrointestinal Stromal Tumor", "Soft Tissue Sarcoma", "Uterine Sarcoma") ~ "Sarcoma",
    `Cancer Type` %in% c("Head and Neck Cancer", "Salivary Gland Cancer") ~ "Head and neck cancer",
    `Cancer Type` == "Hepatobiliary Cancer" ~ "Gallbladder cancer",
    `Cancer Type` == "Melanoma" ~ "Melanoma",
    `Cancer Type` %in% c("Non-Small Cell Lung Cancer", "Small Cell Lung Cancer") ~ "Lung cancer",
    `Cancer Type` == "Ovarian Cancer" ~ "Ovarian cancer",
    `Cancer Type` == "Pancreatic Cancer" ~ "Pancreatic cancer",
    `Cancer Type` == "Prostate Cancer" ~ "Prostate cancer",
    `Cancer Type` == "Renal Cell Carcinoma" ~ "Renal cell cancer",
    `Cancer Type` == "Thyroid Cancer" ~ "Thyroid cancer",
    `Cancer Type Detailed` %in% c("Adenocarcinoma of the Gastroesophageal Junction", "Esophageal Adenocarcinoma",
                                  "Esophageal Squamous Cell Carcinoma", "Esophagogastric Adenocarcinoma") ~ "Esophageal cancer",
    `Cancer Type Detailed` %in% c("Adenosquamous Carcinoma of the Stomach", "Diffuse Type Stomach Adenocarcinoma",
                                  "Stomach Adenocarcinoma", "Poorly Differentiated Carcinoma of the Stomach") ~ "Gastric cancer",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(Tumor)) %>%
  filter((Tumor == "Lung cancer" & `Cancer Type` == "Non-Small Cell Lung Cancer") |
           (Tumor == "Intestinal cancer" & `Cancer Type` == "Colorectal Cancer") |
           !Tumor %in% c("Lung cancer", "Intestinal cancer"))

# Compare cancer types with at least 50 samples in both cohorts.
common_genes <- intersect(unique(axen$Hugo_Symbol), unique(msk$Hugo_Symbol))

axen_n <- axen %>% distinct(Tumor, Tumor_Sample_Barcode) %>%dplyr::count(Tumor, name = "n_axen")
msk_n <- msk %>% distinct(Tumor, Tumor_Sample_Barcode) %>% dplyr::count(Tumor, name = "n_msk")

cohort_n <- inner_join(axen_n, msk_n, by = "Tumor") %>% filter(n_msk >= 30)

axen <- axen %>% semi_join(cohort_n, by = "Tumor")
msk <- msk %>% semi_join(cohort_n, by = "Tumor")
message("Tumor types for Figure 2A: ", paste(cohort_n$Tumor, collapse = ", "))



# Figure 2A: tumor type-gene mutation frequencies -----------------------
# Each sample contributes at most once per gene. Mutation-frequency difference
# is K-MASTER minus MSK-ACCESS (percentage points); P is from Fisher's exact test.
# Keep genes observed in both datasets, as in the original comparison script.
common_genes <- intersect(unique(axen$Hugo_Symbol), unique(msk$Hugo_Symbol))

get_gene_counts <- function(x, sample_col, cohort) {
  n_by_tumor <- x %>% transmute(Tumor, sample_id = .data[[sample_col]]) %>%
    distinct() %>% dplyr::count(Tumor, name = "total")
  
  x %>% filter(Hugo_Symbol %in% common_genes) %>%
    distinct(Tumor, sample_id = .data[[sample_col]], Hugo_Symbol) %>%
    dplyr::count(Tumor, Hugo_Symbol, name = "mutated") %>%
    tidyr::complete(Tumor = n_by_tumor$Tumor, Hugo_Symbol = common_genes, fill = list(mutated = 0)) %>%
    left_join(n_by_tumor, by = "Tumor") %>% mutate(Cohort = cohort)
}

axen_freq <- get_gene_counts(axen, "Tumor_Sample_Barcode", "KM")
msk_freq <- get_gene_counts(msk, "Tumor_Sample_Barcode", "MSK")

fig2a <- inner_join(axen_freq, msk_freq, by = c("Tumor", "Hugo_Symbol"), suffix = c("_km", "_msk")) %>% rowwise() %>%
  mutate(p_value = fisher.test(matrix(c(mutated_km, total_km - mutated_km,
                                        mutated_msk, total_msk - mutated_msk), nrow = 2, byrow = TRUE))$p.value) %>% ungroup() %>%
  mutate(delta_per = 100 * (mutated_km / total_km - mutated_msk / total_msk),
         log10_p = -log10(pmax(p_value, .Machine$double.xmin)), label = paste0(Tumor, "_", Hugo_Symbol))

color_map <- c("Lung cancer" = "#D55E00", "Breast cancer" = "#0072B2", "Intestinal cancer" = "#009E73",
               "Gallbladder cancer" = "#F0E442", "Pancreatic cancer" = "#56B4E9", "Gastric cancer" = "#E69F00",
               "Prostate cancer" = "#B97B60", "Head and neck cancer" = "#D4A10A", "Sarcoma" = "#6C8E3E",
               "Bladder cancer" = "#E1C16E", "Esophageal cancer" = "#F1B1B1", "Neuroendocrine tumor" = "#AD6F6F",
               "Ovarian cancer" = "#AABF82", "Renal cell cancer" = "#5582B0", "Melanoma" = "#6D4583")

# Main volcano plot. The enlarged callouts in the manuscript were assembled
# separately; the labeled version below can be used for that layout.
p2a <- ggplot(fig2a, aes(x = delta_per, y = log10_p)) +
  geom_point(aes(fill = if_else(p_value < 0.05, Tumor, "n.s."), size = total_km),
             shape = 21, color = "black", stroke = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray70") +
  scale_fill_manual(values = c(color_map, "n.s." = "gray85")) +
  scale_size_continuous(range = c(1, 6)) + theme_classic() +
  labs(x = "Mutation frequency difference (KM - MSK, %)", y = expression(-log[10](italic(P))), size = "KM samples", fill = "Tumor type")
ggsave(file.path(outdir, "Figure2A_volcano.pdf"), p2a, width = 9, height = 6)

# Labeled view for selecting pairs in the final manuscript layout.
label_data <- fig2a %>% filter(p_value < 0.05, abs(delta_per) >= 3)
p2a_labels <- p2a + geom_text_repel(data = label_data, aes(label = label), size = 2.5, max.overlaps = 25)
ggsave(file.path(outdir, "Figure2A_labeled.pdf"), p2a_labels, width = 10, height = 7)



# Figure 2B: OncoKB v5.0 alteration frequency matrix --------------------
# This panel uses a curated, precomputed OncoKB comparison table. Its upstream
# alteration-matching and annotation work is NOT regenerated by this script.
# AXEN and MSK frequencies, evidence levels, and Fisher P values are sourced
# from that table; no patient-level data are written to the repository.
stopifnot(file.exists(oncokb_table))
onco <- read.delim(oncokb_table, sep = "\t", check.names = FALSE)
stopifnot(all(c("Gene", "Tumor", "AXEN", "MSK", "Level_edit", "p_value", "label") %in% names(onco)))

six_tumors <- c("Breast cancer", "Lung cancer", "Prostate cancer", "Intestinal cancer", "Gallbladder cancer", "Pancreatic cancer")
tumor_labels <- c("BC", "NSCLC", "PRAD", "CRC", "BTC", "PAAD")
names(tumor_labels) <- six_tumors
label_priority <- c("pancan_oncogenic", "specific_oncogenic", "pancan_specific", "specific_specific")

onco <- onco %>% filter(Tumor %in% six_tumors, (!is.na(AXEN) & AXEN > 0) | (!is.na(MSK) & MSK > 0)) %>%
  mutate(AXEN = as.numeric(AXEN), MSK = as.numeric(MSK), p_value = as.numeric(p_value),
         Level_edit = as.character(Level_edit), label = factor(label, levels = label_priority)) %>%
  arrange(Tumor, Gene, label) %>% distinct(Tumor, Gene, .keep_all = TRUE)

# Paired columns for KM and MSK within each tumor type.
row_order <- sort(unique(onco$Gene))
column_order <- as.vector(t(outer(six_tumors, c("KM", "MSK"), paste, sep = "_")))
plot_df <- onco %>% pivot_longer(c(AXEN, MSK), names_to = "Cohort", values_to = "frequency") %>%
  mutate(Cohort = recode(Cohort, AXEN = "KM"), Tumor_Cohort = paste(Tumor, Cohort, sep = "_")) %>%
  dplyr::select(Gene, Tumor_Cohort, frequency, Level_edit, p_value) %>%
  tidyr::complete(Gene = row_order, Tumor_Cohort = column_order)

make_matrix <- function(x, field) {
  m <- x %>% dplyr::select(Gene, Tumor_Cohort, all_of(field)) %>%
    pivot_wider(names_from = Tumor_Cohort, values_from = all_of(field)) %>% as.data.frame()
  rownames(m) <- m$Gene
  m$Gene <- NULL
  as.matrix(m[row_order, column_order, drop = FALSE])
}

freq_mat <- make_matrix(plot_df, "frequency")
level_mat <- make_matrix(plot_df, "Level_edit")
p_mat <- make_matrix(plot_df, "p_value")
delta_mat <- make_matrix(plot_df, "frequency")
for (tumor in six_tumors) {
  km_col <- paste0(tumor, "_KM"); msk_col <- paste0(tumor, "_MSK")
  delta_mat[, km_col] <- freq_mat[, km_col] - freq_mat[, msk_col]
  delta_mat[, msk_col] <- delta_mat[, km_col]
}

level_colors <- c("1" = "#44B549", "2" = "#3768B2", "3" = "#6E2E91", "4" = "#626164")
shade_color <- function(freq, lvl, max_freq = 30) {
  if (is.na(freq) || freq <= 0 || is.na(lvl) || !lvl %in% names(level_colors)) return("white")
  base <- grDevices::col2rgb(level_colors[[lvl]])[, 1] / 255
  fraction <- min(1, log1p(freq) / log1p(max_freq))
  grDevices::rgb(1 - fraction * (1 - base[1]), 1 - fraction * (1 - base[2]), 1 - fraction * (1 - base[3]))
}

ht <- Heatmap(freq_mat, name = "Mutation frequency (%)", cluster_rows = FALSE, cluster_columns = FALSE,
              show_heatmap_legend = FALSE, column_split = factor(rep(tumor_labels, each = 2), levels = tumor_labels),
              column_labels = rep(c("KM", "MSK"), length(six_tumors)), column_names_rot = 0,
              row_names_side = "right", column_gap = unit(2, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                f <- freq_mat[i, j]; lvl <- level_mat[i, j]; p <- p_mat[i, j]; delta <- delta_mat[i, j]
                border <- if (!is.na(p) && p < 0.05 && !is.na(delta)) {
                  if (delta > 0) "#E51813" else if (delta < 0) "#18499E" else "white"
                } else "gray90"
                grid.rect(x, y, width, height, gp = gpar(fill = shade_color(f, lvl), col = border,
                                                         lwd = ifelse(border %in% c("#E51813", "#18499E"), 1.6, 0.3)))
                if (!is.na(f) && f > 0) grid.text(format(round(f, 1), trim = TRUE), x, y, gp = gpar(fontsize = 7))
              })

pdf(file.path(outdir, "Figure2B_OncoKB_matrix.pdf"), width = 12, height = 12)
draw(ht, heatmap_legend_list = list(Legend(title = "OncoKB level", at = names(level_colors),
                                           legend_gp = gpar(fill = unname(level_colors)))))
dev.off()
