
# Figure 3: Comparison of tissue and ctDNA genomic profiles
# Analysis code accompanying the K-MASTER pan-cancer ctDNA manuscript.
# Controlled-access genomic and clinical input files are not included in this repository.
# Starts from project-filtered variant tables, not raw sequencing reads.

# Packages ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(readxl)
library(ggrepel)
library(ggforce)
library(ggnewscale)
library(RColorBrewer)

# Input files ------------------------------------------------------------
# Final filtered ctDNA MAF and tissue variant table (mutation-positive samples).
axen_file <- "KM00_AXEN.maf.gz"
tissue_file <- "tissue_mutation_maf.txt"

# Less stringently filtered files used to define mutation-frequency denominators:
axen_denominator_file <- "KM00_AXEN_filterpass.maf"
tissue_denominator_file <- "tissue_filterpass.txt"

# K-MASTER clinical annotation, shared 87-gene panel list, pathway and OncoKB tables.
clinical_file <- "KM00_CLINICAL.xlsx"
clinical_subtype_file <- "KKM00_CLINICAL_subtupe.xlsx"
panel_file <- "K-MASTER_genelist.xlsx"
pathway_file <- "oncogenic_sig_pathways.xlsx"
oncokb_file <- "oncokb_biomarker.xlsx"

# Sample-selection metadata used for the paired donut plots (Fig. 3F).
sample_info_file <- "KM00_sample.xlsx"
outdir <- "results/Figure3"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Data preparation -------------------------------------------------------
clinical <- read_xlsx(clinical_file, skip = 1) %>%
  rename(`SUBJECT ID` = `Subject ID`)
clinical_subtype <- read_xlsx(clinical_subtype_file, skip = 1) %>%
  rename(`SUBJECT ID` = `Subject ID`)

# The original input files are project-specific. Attach annotations only when absent.
add_annotations <- function(x, metadata, fields) {
  fields <- intersect(fields, names(metadata))
  needed <- setdiff(fields, names(x))
  if (length(needed) > 0) {
    x <- x %>% left_join(metadata %>% dplyr::select(`SUBJECT ID`, all_of(needed)) %>%
                           distinct(`SUBJECT ID`, .keep_all = TRUE), by = "SUBJECT ID")
  }
  x
}

axen <- read.delim(axen_file, sep = "\t", check.names = FALSE) %>%
  add_annotations(clinical, c("Tumor")) %>%
  filter(!Tumor %in% c("Double cancer", "Other"))
tissue <- read.delim(tissue_file, sep = "\t", check.names = FALSE) %>%
  add_annotations(clinical, c("Tumor", "PANEL")) %>%
  filter(!Tumor %in% c("Double cancer", "Other"))

# FIRST was excluded in the original tissue-sequencing comparison.
panel_col <- if ("PANEL.x" %in% names(tissue)) "PANEL.x" else "PANEL"
stopifnot(panel_col %in% names(tissue))
tissue <- tissue %>% filter(.data[[panel_col]] != "FIRST")

# Tumor types with at least 15 patients in each filtered cohort.
axen_n <- axen %>% distinct(Tumor, `SUBJECT ID`) %>% dplyr::count(Tumor, name = "axen_n")
tissue_n <- tissue %>% distinct(Tumor, `SUBJECT ID`) %>% dplyr::count(Tumor, name = "tissue_n")
tumors <- inner_join(axen_n, tissue_n, by = "Tumor") %>%
  filter(axen_n >= 15, tissue_n >= 15) %>% pull(Tumor)
stopifnot(length(tumors) == 13)

panel_genes <- read_xlsx(panel_file, sheet = "Sheet1") %>% pull(`K&C&A_87`) %>%
  as.character() %>% na.omit() %>% unique()
axen <- axen %>% filter(Tumor %in% tumors, Hugo_Symbol %in% panel_genes)
tissue <- tissue %>% filter(Tumor %in% tumors, Hugo_Symbol %in% panel_genes)

# Unlike the filtered variant tables, these files retain samples for denominators.
axen_base <- read.delim(axen_denominator_file, sep = "\t", check.names = FALSE) %>%
  add_annotations(clinical, c("Tumor")) %>% filter(Tumor %in% tumors)
tissue_base <- read.delim(tissue_denominator_file, sep = "\t", check.names = FALSE) %>%
  add_annotations(clinical, c("Tumor")) %>% filter(Tumor %in% tumors)

axen_totals <- axen_base %>% distinct(Tumor, Tumor_Sample_Barcode) %>%
  dplyr::count(Tumor, name = "liquid_n")
tissue_totals <- tissue_base %>% distinct(Tumor, Tumor_Sample_Barcode) %>%
  dplyr::count(Tumor, name = "tissue_n")

# Common plotting colors
tumor_colors <- c("Lung cancer" = "#D55E00", "Breast cancer" = "#0072B2",
                  "Intestinal cancer" = "#009E73", "Gallbladder cancer" = "#F0E442",
                  "Pancreatic cancer" = "#56B4E9", "Gastric cancer" = "#E69F00",
                  "Prostate cancer" = "#B97B60", "Head and neck cancer" = "#D4A10A",
                  "Sarcoma" = "#6C8E3E", "Bladder cancer" = "#E1C16E",
                  "Liver cancer" = "#9C3B3B", "Renal cell cancer" = "#5582B0",
                  "Esophageal cancer" = "#F1B1B1")



# Figure 3A-B: gene-level mutation frequency ----------------------------
# Count a gene only once per sample. Denominators use the broader QC-filtered files.
count_gene <- function(maf, sample_n, count_name, total_name) {
  maf %>% distinct(Tumor, `SAMPLE ID`, Hugo_Symbol) %>%
    dplyr::count(Tumor, Hugo_Symbol, name = count_name) %>%
    tidyr::complete(Tumor = tumors, Hugo_Symbol = panel_genes,
                    fill = setNames(list(0L), count_name)) %>%
    left_join(sample_n, by = "Tumor") %>%
    mutate(frequency = 100 * .data[[count_name]] / .data[[total_name]])
}

liquid_freq <- count_gene(axen, axen_totals, "liquid_mut", "liquid_n") %>%
  rename(liquid_freq = frequency)
tissue_freq <- count_gene(tissue, tissue_totals, "tissue_mut", "tissue_n") %>%
  rename(tissue_freq = frequency)

freq <- inner_join(liquid_freq, tissue_freq,
                   by = c("Tumor", "Hugo_Symbol")) %>%
  mutate(delta_pp = tissue_freq - liquid_freq)
freq$p_value <- vapply(seq_len(nrow(freq)), function(i) {
  fisher.test(matrix(c(freq$liquid_mut[i], freq$liquid_n[i] - freq$liquid_mut[i],
                       freq$tissue_mut[i], freq$tissue_n[i] - freq$tissue_mut[i]),
                     nrow = 2, byrow = TRUE))$p.value
}, numeric(1))
freq <- freq %>% mutate(neg_log2_p = -log2(pmax(p_value, .Machine$double.xmin)))

# A: bubble plot. The original plotting workflow suppressed comparisons with
# <5 percentage-point differences and non-significant Fisher P values.
tumor_order <- c("Intestinal cancer", "Lung cancer", "Pancreatic cancer", "Breast cancer",
                 "Gallbladder cancer", "Gastric cancer", "Sarcoma", "Bladder cancer",
                 "Head and neck cancer", "Renal cell cancer", "Prostate cancer",
                 "Esophageal cancer", "Liver cancer")
fig3a_data <- freq %>% filter(abs(delta_pp) >= 5, p_value < 0.05) %>%
  mutate(Tumor = factor(Tumor, levels = rev(tumor_order)))
p3a <- ggplot(fig3a_data, aes(x = Hugo_Symbol, y = Tumor, size = neg_log2_p, fill = delta_pp)) +
  geom_point(shape = 21, color = "gray35", alpha = 0.85) +
  scale_fill_gradient2(low = "#B84D51", mid = "white", high = "#4C65B1", midpoint = 0) +
  scale_size_continuous(range = c(1.5, 7), name = "-log2(P)") +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, fill = "Tissue - liquid (%)") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.minor = element_blank())
ggsave(file.path(outdir, "Figure3A_mutation_frequency_bubble.pdf"), p3a, width = 14, height = 4.5)



# B: correlation across tumor type-gene pairs.
cor_result <- cor.test(freq$liquid_freq, freq$tissue_freq, method = "spearman", exact = FALSE)
fig3b_data <- freq %>% mutate(label = if_else(
  paste0(Tumor, "_", Hugo_Symbol) %in% c(
    "Esophageal cancer_TP53", "Pancreatic cancer_KRAS", "Intestinal cancer_TP53",
    "Breast cancer_TP53", "Gallbladder cancer_TP53", "Gastric cancer_TP53",
    "Breast cancer_PIK3CA", "Breast cancer_BRCA1", "Breast cancer_ERBB2",
    "Lung cancer_EGFR", "Breast cancer_ESR1"), Hugo_Symbol, ""))
p3b <- ggplot(fig3b_data, aes(x = liquid_freq, y = tissue_freq)) +
  geom_point(aes(fill = Tumor, size = liquid_n), shape = 21, color = "gray30", alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray65") +
  ggrepel::geom_text_repel(data = fig3b_data %>% filter(label != ""),
                           aes(label = paste(Tumor, label)), size = 2.7, max.overlaps = Inf) +
  scale_fill_manual(values = tumor_colors) + scale_size_continuous(range = c(1, 4)) +
  labs(x = "Liquid biopsy mutation frequency (%)", y = "Tissue mutation frequency (%)",
       subtitle = sprintf("Spearman R = %.2f, P = %.2g", unname(cor_result$estimate), cor_result$p.value)) +
  theme_classic()
ggsave(file.path(outdir, "Figure3B_mutation_frequency_correlation.pdf"), p3b, width = 7, height = 6)



# Figure 3C-D: oncogenic pathway profiles -------------------------------
# Original pathway definitions were restricted to Cell_paper == 1 and
# excluded the listed genes. NRF2 and Hippo were not plotted.
pathways <- read_xlsx(pathway_file) %>%
  filter(Cell_paper == 1,
         !Gene %in% c("ETV1", "ETV4", "ETV5", "FLT3", "GNAQ", "GNAS", "JAK1", "JAK2", "JAK3",
                      "MAP2K1", "MAP2K2", "MAP2K4", "MPL", "NPM1", "RUNX1", "SMARCB1"),
         !Pathway %in% c("NRF2", "Hippo")) %>%
  distinct(Pathway, Gene)

pathway_frequency <- function(maf, totals, total_col, frequency_col) {
  changed <- maf %>% distinct(Tumor, Tumor_Sample_Barcode, Hugo_Symbol) %>%
    inner_join(pathways, by = c("Hugo_Symbol" = "Gene")) %>%
    distinct(Tumor, Tumor_Sample_Barcode, Pathway) %>%
    dplyr::count(Tumor, Pathway, name = "mutated")
  changed %>% tidyr::complete(Tumor = tumors, Pathway = unique(pathways$Pathway), fill = list(mutated = 0L)) %>%
    left_join(totals, by = "Tumor") %>%
    mutate(frequency = 100 * mutated / .data[[total_col]]) %>%
    dplyr::select(Tumor, Pathway, frequency) %>% rename(!!frequency_col := frequency)
}

pathway_freq <- inner_join(
  pathway_frequency(axen, axen_totals, "liquid_n", "Liquid"),
  pathway_frequency(tissue, tissue_totals, "tissue_n", "Tissue"), by = c("Tumor", "Pathway")) %>%
  pivot_longer(c(Liquid, Tissue), names_to = "Source", values_to = "Frequency")


# C: IC and PRAD; D: LC and BC. Bars remain on the original percent scale.
plot_pathways <- function(tumor, panel) {
  dat <- pathway_freq %>% filter(Tumor == tumor) %>%
    mutate(Source = factor(Source, levels = c("Liquid", "Tissue")))
  p <- ggplot(dat, aes(x = Pathway, y = Frequency, fill = Source)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75, color = "white") +
    coord_polar() + scale_fill_manual(values = c("Liquid" = "#18499E", "Tissue" = "#E51813")) +
    scale_y_continuous(limits = c(0, 75), breaks = c(25, 50, 75)) +
    labs(title = tumor, x = NULL, y = NULL) + theme_minimal() +
    theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 9),
          legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  ggsave(file.path(outdir, paste0("Figure3", panel, "_", gsub(" ", "_", tumor), "_pathways.pdf")), p, width = 5, height = 5)
}
plot_pathways("Intestinal cancer", "C")
plot_pathways("Prostate cancer", "C")
plot_pathways("Lung cancer", "D")
plot_pathways("Breast cancer", "D")



# Figure 3E: patient-level OncoKB evidence levels ------------------------
# V5.0 workbook. The original revised 3E excluded KRAS/NRAS wild-type
# inference for colorectal cancer; other alteration rules are retained.
onco_e <- read_xlsx(oncokb_file, sheet = "v5.0") %>%
  mutate(Level = as.character(Level), Tumor = as.character(Tumor),
         Gene = as.character(Gene), Alterations = as.character(Alterations)) %>%
  filter(Level %in% c("1", "2", "3", "4"), Gene %in% panel_genes, !is.na(Tumor), Tumor != "Liquid",
         !(Gene %in% c("KRAS", "NRAS") & Alterations == "Wildtype" & Tumor == "Intestinal cancer" & Subtype == "CRC")) %>%
  separate_rows(Alterations, sep = ",\\s*") %>%
  distinct(Tumor, Gene, Alterations, Level, category, .keep_all = TRUE)
# The original Figure 3F did not annotate wild-type biomarkers.
onco_f <- onco_e %>% filter(!Alterations %in% "Wildtype")

add_subtype <- function(x, source) {
  lung_codes <- c("LCLC", "LUAD", "LUSC", "NSCLCPD")
  if (source == "Tissue") lung_codes <- c(lung_codes, "LUAS")
  x %>% add_annotations(clinical_subtype, c("Tumor_subtype_JW")) %>%
    filter((Tumor == "Intestinal cancer" & Tumor_subtype_JW == "COADREAD") |
             (Tumor == "Lung cancer" & Tumor_subtype_JW %in% lung_codes) |
             !Tumor %in% c("Intestinal cancer", "Lung cancer")) %>%
    mutate(HGVSp_Short_edit = sub("^p\\.", "", HGVSp_Short))
}
axen_onco <- add_subtype(axen, "Liquid")
tissue_onco <- add_subtype(tissue, "Tissue")

# Annotation lookup is restricted to matching tumor types or pan-cancer rules.
annotate_variants <- function(maf, onco) {
  alt <- onco %>% filter(!Alterations %in% c("Oncogenic Mutations", "Truncating Mutations", "Wildtype")) %>%
    dplyr::select(Onco_Tumor = Tumor, Gene, Alterations, Level)
  gene <- onco %>% filter(Alterations == "Oncogenic Mutations") %>%
    dplyr::select(Onco_Tumor = Tumor, Gene, Level)
  trunc <- onco %>% filter(Alterations == "Truncating Mutations") %>%
    dplyr::select(Onco_Tumor = Tumor, Gene, Level)
  match_alt <- maf %>% inner_join(alt, by = c("Hugo_Symbol" = "Gene", "HGVSp_Short_edit" = "Alterations"), relationship = "many-to-many") %>%
    mutate(annotation_type = "specific_mutation")
  match_gene <- maf %>% inner_join(gene, by = c("Hugo_Symbol" = "Gene"), relationship = "many-to-many") %>%
    mutate(annotation_type = "all_mutation")
  match_trunc <- maf %>% filter(Variant_Classification == "Nonsense_Mutation") %>%
    inner_join(trunc, by = c("Hugo_Symbol" = "Gene"), relationship = "many-to-many") %>%
    mutate(annotation_type = "truncating_mutation")
  # For each variant, prefer cancer-specific annotations over pan-cancer,
  # and specific alterations over gene-level annotations.
  bind_rows(match_alt, match_gene, match_trunc) %>%
    filter(Onco_Tumor == "Pan-cancer" | Onco_Tumor == Tumor) %>%
    mutate(Level = as.character(Level),
           priority = case_when(Onco_Tumor == Tumor & annotation_type == "specific_mutation" ~ 1L,
                                Onco_Tumor == Tumor & annotation_type == "all_mutation" ~ 2L,
                                Onco_Tumor == "Pan-cancer" & annotation_type == "specific_mutation" ~ 3L,
                                Onco_Tumor == "Pan-cancer" & annotation_type == "all_mutation" ~ 4L,
                                TRUE ~ 5L)) %>%
    arrange(priority, as.integer(Level)) %>%
    distinct(Tumor, `SUBJECT ID`, Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short_edit, .keep_all = TRUE) %>%
    dplyr::select(`SUBJECT ID`, Tumor, Tumor_Sample_Barcode, Hugo_Symbol,
                  HGVSp_Short_edit, Level, annotation_type)
}

level_colors <- c("1" = "#44B549", "2" = "#3768B2", "3" = "#6E2E91", "4" = "#626164", "none" = "#D7D7D7")

patient_levels <- function(maf, source) {
  all_patients <- maf %>% distinct(Tumor, `SUBJECT ID`)
  matched <- annotate_variants(maf, onco_e) %>%
    mutate(Level_num = as.integer(Level)) %>% arrange(`SUBJECT ID`, Level_num) %>%
    distinct(`SUBJECT ID`, .keep_all = TRUE) %>% dplyr::select(`SUBJECT ID`, Level)
  all_patients %>% left_join(matched, by = "SUBJECT ID") %>%
    mutate(Level = replace_na(Level, "none"), Source = source)
}
level_patients <- bind_rows(patient_levels(axen_onco, "Liquid"),
                            patient_levels(tissue_onco, "Tissue"))
level_summary <- level_patients %>% dplyr::count(Tumor, Source, Level, name = "n") %>%
  group_by(Tumor, Source) %>% mutate(percent = 100 * n / sum(n)) %>% ungroup() %>%
  tidyr::complete(Tumor, Source, Level = names(level_colors),
                  fill = list(n = 0, percent = 0))

level_test <- level_summary %>% group_by(Tumor) %>% group_modify(~{
  wide <- .x %>% dplyr::select(Source, Level, n) %>%
    pivot_wider(names_from = Level, values_from = n, values_fill = 0)
  mat <- as.matrix(wide[, -1, drop = FALSE]); rownames(mat) <- wide$Source
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  tibble(p_value = if (ncol(mat) > 1) fisher.test(mat)$p.value else NA_real_)
}) %>% ungroup() %>% mutate(FDR = p.adjust(p_value, method = "BH"))

row_order <- unlist(lapply(tumor_order, function(t) paste(t, c("Liquid", "Tissue"), sep = " - ")))
level_plot <- level_summary %>% mutate(Source = factor(Source, levels = c("Liquid", "Tissue")),
                                       Level = factor(Level, levels = c("none", "4", "3", "2", "1")),
                                       Tumor_Source = factor(paste(Tumor, Source, sep = " - "), levels = rev(row_order)))
signif_labels <- level_test %>% filter(!is.na(p_value), p_value < 0.05) %>%
  mutate(Tumor_Source = factor(paste(Tumor, "Liquid", sep = " - "), levels = rev(row_order)))
p3e <- ggplot(level_plot, aes(x = percent, y = Tumor_Source, fill = Level)) +
  geom_col(width = 0.8) + scale_fill_manual(values = level_colors) +
  geom_text(data = signif_labels, aes(x = 103, y = Tumor_Source, label = "*"),
            inherit.aes = FALSE, size = 4) + coord_cartesian(xlim = c(0, 108)) +
  labs(x = "Patients (%)", y = NULL, fill = "OncoKB level") + theme_minimal() +
  theme(panel.grid.major.y = element_blank())
ggsave(file.path(outdir, "Figure3E_OncoKB_levels.pdf"), p3e, width = 8, height = 9)



# Figure 3F: paired OncoKB alteration-composition donut plots ------------
# Select one sequencing sample per patient using the study's focus flags.
# Inner ring: highest OncoKB evidence level per patient.
# Outer ring: alterations at that level; a patient may contribute multiple
# alteration labels, matching the revised original Figure 3F workflow.
sample_info <- read_xlsx(sample_info_file, sheet = "KM00")

select_focus_samples <- function(maf, focus_col) {
  stopifnot(all(c("SUBJECT ID", "SGI ID", focus_col) %in% names(sample_info)))
  focus <- sample_info %>% filter(trimws(.data[[focus_col]]) == "O") %>%
    transmute(`SUBJECT ID` = as.character(`SUBJECT ID`),
              Tumor_Sample_Barcode = as.character(`SGI ID`)) %>% distinct()
  maf %>% mutate(`SUBJECT ID` = as.character(`SUBJECT ID`),
                 Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode)) %>%
    inner_join(focus, by = c("SUBJECT ID", "Tumor_Sample_Barcode"))
}

make_donut_data <- function(maf, source) {
  maf <- select_focus_samples(maf, if (source == "Liquid") "FINAL_AXEN_focus" else "FINAL")
  n_patients <- maf %>% distinct(Tumor, `SUBJECT ID`) %>% dplyr::count(Tumor, name = "total_n")
  matched <- annotate_variants(maf, onco_f) %>%
    mutate(AlterLabel = if_else(annotation_type == "specific_mutation",
                                paste(Hugo_Symbol, HGVSp_Short_edit), Hugo_Symbol),
           Level_num = as.integer(Level))
  matched <- matched %>% mutate(AlterLabel = case_when(
    Tumor == "Lung cancer" & AlterLabel %in% c("EGFR E746_A750del", "EGFR L747_S752del",
                                               "EGFR L747_T751del", "EGFR L747_P753delinsS", "EGFR L747_E749del") ~ "EGFR Exon 19 Del",
    Tumor == "Lung cancer" & AlterLabel %in% c("EGFR A767_V769dup", "EGFR S768_D770dup") ~ "EGFR Exon 20 Ins",
    TRUE ~ AlterLabel))
  top_level <- matched %>% group_by(Tumor, `SUBJECT ID`) %>%
    summarise(Level_num = min(Level_num), .groups = "drop")
  all_levels <- maf %>% distinct(Tumor, `SUBJECT ID`) %>%
    left_join(top_level, by = c("Tumor", "SUBJECT ID")) %>%
    mutate(Level = if_else(is.na(Level_num), "none", as.character(Level_num)))
  inner <- all_levels %>% dplyr::count(Tumor, Level, name = "n") %>%
    group_by(Tumor) %>% mutate(prop = 100 * n / sum(n), end = cumsum(prop), start = lag(end, default = 0)) %>% ungroup() %>%
    mutate(Source = source)
  outer <- matched %>% inner_join(top_level, by = c("Tumor", "SUBJECT ID", "Level_num")) %>%
    filter(Level_num %in% 1:3) %>% distinct(Tumor, `SUBJECT ID`, Level_num, AlterLabel) %>%
    dplyr::count(Tumor, Level_num, AlterLabel, name = "n") %>%
    left_join(n_patients, by = "Tumor") %>%
    mutate(Level = as.character(Level_num), freq = 100 * n / total_n) %>%
    arrange(Tumor, Level, AlterLabel) %>% group_by(Tumor, Level) %>%
    mutate(frac = n / sum(n), seg_end = cumsum(frac), seg_start = lag(seg_end, default = 0)) %>%
    ungroup() %>% left_join(inner %>% dplyr::select(Tumor, Level, level_start = start,
                                                    level_end = end), by = c("Tumor", "Level")) %>%
    mutate(start = level_start + seg_start * (level_end - level_start),
           end = level_start + seg_end * (level_end - level_start), Source = source)
  list(inner = inner, outer = outer)
}

liquid_donut <- make_donut_data(axen_onco, "Liquid")
tissue_donut <- make_donut_data(tissue_onco, "Tissue")
inner_all <- bind_rows(liquid_donut$inner, tissue_donut$inner)
outer_all <- bind_rows(liquid_donut$outer, tissue_donut$outer)

donut_tumors <- c("Breast cancer", "Lung cancer", "Prostate cancer", "Intestinal cancer", "Sarcoma", "Pancreatic cancer")
inner_all <- inner_all %>% filter(Tumor %in% donut_tumors)
outer_all <- outer_all %>% filter(Tumor %in% donut_tumors)

# Plotting source data are generated from patient-level evidence assignments.
make_donut_plot <- function(tumor) {
  inner <- inner_all %>% filter(Tumor == tumor) %>%
    mutate(Source = factor(Source, levels = c("Liquid", "Tissue")))
  outer <- outer_all %>% filter(Tumor == tumor) %>%
    mutate(Source = factor(Source, levels = c("Liquid", "Tissue")),
           angle = pi / 2 - pi * (start + end) / 100,
           label_x = 1.08 * cos(angle), label_y = 1.08 * sin(angle))
  alteration_cols <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(
    max(3, n_distinct(outer$AlterLabel))), sort(unique(outer$AlterLabel)))
  ring_bg <- data.frame(Source = factor(c("Liquid", "Tissue"),
                                        levels = c("Liquid", "Tissue")), start = 0, end = 2 * pi)
  p <- ggplot() +
    ggforce::geom_arc_bar(data = ring_bg, aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1,
                                              start = start, end = end), fill = "#D7D7D7", color = "white") +
    ggforce::geom_arc_bar(data = inner, aes(x0 = 0, y0 = 0, r0 = 0.4, r = 0.7,
                                            start = 2 * pi * start / 100, end = 2 * pi * end / 100, fill = Level), color = "white") +
    scale_fill_manual(values = level_colors) + ggnewscale::new_scale_fill() +
    ggforce::geom_arc_bar(data = outer, aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1,
                                            start = 2 * pi * start / 100, end = 2 * pi * end / 100, fill = AlterLabel), color = "white") +
    scale_fill_manual(values = alteration_cols) +
    ggrepel::geom_text_repel(data = outer %>% filter(freq >= 2),
                             aes(x = label_x, y = label_y, label = paste0(AlterLabel, "\n", sprintf("%.1f%%", freq))),
                             size = 2.4, max.overlaps = Inf) + facet_wrap(~Source, nrow = 1) +
    coord_fixed() + theme_void() + theme(legend.position = "none")
  ggsave(file.path(outdir, paste0("Figure3F_", gsub(" ", "_", tumor), "_OncoKB_donuts.pdf")),
         p, width = 9, height = 4.5)
}
invisible(lapply(donut_tumors, make_donut_plot))

# Original project data are controlled access. Verify panel-specific counts and
# plot statistics against the final manuscript before releasing this script.
