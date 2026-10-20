
# Figure 1: K-MASTER pan-cancer ctDNA genomic landscape
# Code accompanying: Pan-cancer circulating tumor DNA atlas establishes a
# precision oncology framework for prognostic and therapeutic stratification
# Input files are controlled-access or project-specific and are not distributed here.
# This script starts from a filtered K-MASTER AXEN MAF, not from raw sequencing data.



# Packages ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(readxl)
library(lubridate)
library(maftools)
library(data.table)
library(circlize)
library(ComplexHeatmap)
library(grid)
library(RColorBrewer)



# Input files ------------------------------------------------------------
# Filtered, annotated AXEN ctDNA MAF used in the original Figure 1 analysis.
# Includes Tumor, SUBJECT ID, SAMPLE ID, Tumor_Sample_Barcode, mutation calls, VAF.
maf_file <- "KM00_AXEN.maf.gz"

# K-MASTER clinical annotations (patient characteristics and tumor subtype).
clinical_file <- "KM00_CLINICAL.xlsx"

# Gene coverage list for the 87 genes shared by the K-MASTER panels.
panel_file <- "K-MASTER_panels_genelist.xlsx"

# Curated pathway-frequency table used in the original Figure 1B.
# Expected columns: Cancer, Pathway, AXEN_per, Tissue_per.
pathway_file <- "mutational_pathway.txt"

# Curated OncoKB v5.0 annotations used in the original Figure 1F.
# This is the study-specific annotation workbook, not a file downloaded by this script.
oncokb_file <- "oncokb_biomarker.txt"

outdir <- "results/Figure1"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)



# Cohort preparation -----------------------------------------------------
axen_final2 <- read.delim(maf_file, sep = "\t", check.names = FALSE)
axen_final2 <- axen_final2 %>%
  filter(!Tumor %in% c("Double cancer", "Other"))

# Include tumor types with >=10 unique patients (17 cancer types).
selected_tumors <- axen_final2 %>% distinct(Tumor, `SUBJECT ID`) %>%
  dplyr::count(Tumor, name = "n_patients") %>% filter(n_patients >= 10) %>% pull(Tumor)
axen_final3 <- axen_final2 %>% filter(Tumor %in% selected_tumors)

stopifnot(length(unique(axen_final3$Tumor)) == 17)
# The published cohort contains 1,225 patients and 1,243 ctDNA samples.
message("Patients: ", n_distinct(axen_final3$`SUBJECT ID`), "; samples: ", n_distinct(axen_final3$Tumor_Sample_Barcode))

# Supplement genomic records with clinical annotations, without duplicating Tumor.
cli <- readxl::read_xlsx(clinical_file, skip = 1)
if (!"SUBJECT ID" %in% names(cli) && "Subject ID" %in% names(cli)) {
  cli <- cli %>% rename(`SUBJECT ID` = `Subject ID`)
}
# Convert Excel serial dates before merging clinical data.
cli$Diagnosis_date_edit <- suppressWarnings(as.numeric(cli$Diagnosis_date_edit))
cli$Diagnosis_date_edit <- as.Date(cli$Diagnosis_date_edit, origin = "1899-12-30")

cli$Survival_250417 <- stringi::stri_replace_all_regex(cli$Survival_250417, vectorize = F,
                                                       pattern = c("사망", "생존 중", "FU loss", "생존", "연구탈락"),
                                                       replacement = c(1, 0, 0, 0, NA))

required_clinical <- c("SUBJECT ID", "Gender", "Tumor_subtype_JW", "NGS_stage",
                       "Birth_date", "Diagnosis_date_edit")
stopifnot(all(required_clinical %in% names(cli)))
# Prefer the clinical workbook for these fields if an existing MAF already has them.
axen_final3 <- axen_final3 %>%
  dplyr::select(-any_of(setdiff(required_clinical, "SUBJECT ID"))) %>%
  left_join(cli %>% dplyr::select(all_of(required_clinical)) %>%
              distinct(`SUBJECT ID`, .keep_all = TRUE), by = "SUBJECT ID")

# Convert birth and diagnosis dates; Excel numeric dates are supported.
to_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  suppressWarnings(as.Date(lubridate::ymd(as.character(x))))
}
axen_final3 <- axen_final3 %>%
  mutate(Birth_date = as.Date(Birth_date), Diagnosis_date_edit = as.Date(Diagnosis_date_edit),
    Age = floor(lubridate::time_length(lubridate::interval(Birth_date, Diagnosis_date_edit), "years")),
    NGS_stage_edit = as.character(NGS_stage), NGS_stage_edit = sub("^([1-4])[A-C]$", "\\1", NGS_stage_edit))

# Mutation counts follow the nonsynonymous MAF calls included by maftools.
maf_obj <- maftools::read.maf(maf = axen_final3, verbose = FALSE)
mutation_counts <- data.frame(Tumor_Sample_Barcode = as.character(maf_obj@variant.type.summary$Tumor_Sample_Barcode),
                              count = as.numeric(maf_obj@variant.type.summary$total))
axen_final3 <- axen_final3 %>%
  left_join(mutation_counts, by = "Tumor_Sample_Barcode") %>%
  mutate(count_edit = if_else(count >= 3, "mut-H", "mut-L", missing = NA_character_))

# Color definitions used across panels.
tumor_colors <- c("Lung cancer" = "#D55E00", "Breast cancer" = "#0072B2",
                  "Intestinal cancer" = "#009E73", "Gallbladder cancer" = "#F0E442",
                  "Pancreatic cancer" = "#56B4E9", "Gastric cancer" = "#E69F00",
                  "Prostate cancer" = "#B97B60", "Head and neck cancer" = "#D4A10A",
                  "Sarcoma" = "#6C8E3E", "Bladder cancer" = "#E1C16E",
                  "Liver cancer" = "#9C3B3B", "Renal cell cancer" = "#5582B0",
                  "Esophageal cancer" = "#F1B1B1", "Neuroendocrine tumor" = "#AD6F6F",
                  "Ovarian cancer" = "#AABF82", "Thyroid cancer" = "#7F5A6D",
                  "Cervical cancer" = "#4D8AB1", "Melanoma" = "#6D4583", "Thymoma" = "#5F4B8B")



# Figure 1A: patient and tumor characteristics --------------------------
vaf_df <- axen_final3 %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(max_vaf = max(as.numeric(VAF), na.rm = TRUE), Tumor = first(Tumor), Gender = first(Gender),
    Subtype = first(Tumor_subtype_JW), Stage = first(NGS_stage_edit),
    Age = first(Age), .groups = "drop"
  ) %>%
  mutate(
    Subtype = dplyr::recode(Subtype,
      "Nasopharyngeal Carcinoma (NPC)" = "NPC", "Kidney(RCC)" = "RCC",
      "Melanoma" = "MM", "Melanoma(no MEL)" = "MM",
      "Synovial Sarcoma (SYNS)" = "SYNS", "Urothelial_Carcinoma" = "Other"),
    Subtype = if_else(is.na(Subtype), "Other", Subtype)
  )

tumor_counts <- vaf_df %>% dplyr::count(Tumor, name = "n_samples") %>% arrange(desc(n_samples))
vaf_df <- vaf_df %>%
  mutate(sector = factor(Tumor, levels = tumor_counts$Tumor)) %>%
  group_by(sector) %>% mutate(x = row_number()) %>% ungroup()

subtype_palettes <- c("Accent", "Set1", "Set2", "Set3", "Paired", "Dark2", "Pastel1", "Pastel2")
subtype_palette <- unlist(lapply(subtype_palettes, function(p) {
  RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[p, "maxcolors"], p)
}))
subtype_levels <- c(sort(setdiff(unique(na.omit(vaf_df$Subtype)), "Other")), "Other")
subtype_colors <- setNames(c(subtype_palette[seq_len(length(subtype_levels) - 1)], "#CCCCCC"), subtype_levels)
sex_colors <- c("F" = "#E0CA4D", "M" = "#7C94E1")
stage_colors <- c("4" = "#AF6378", "3" = "#D7A9AB", "2" = "#E0C073", "1" = "#DEE2A5")

pdf(file.path(outdir, "Figure1A_circos.pdf"), width = 12, height = 12)
circos.clear()
circos.par(start.degree = 60, cell.padding = c(0.005, 0.2, 0.005, 0.2),
           track.margin = c(0.005, 0.005),
           gap.after = c(rep(1, nrow(tumor_counts) - 1), 30))
circos.initialize(factors = tumor_counts$Tumor,
                  xlim = cbind(rep(0, nrow(tumor_counts)), tumor_counts$n_samples))

# Track 1: tumor type.
circos.track(factors = vaf_df$sector, x = vaf_df$x, y = rep(1, nrow(vaf_df)),
             ylim = c(0, 1), track.height = 0.07, bg.border = NA,
             panel.fun = function(x, y) {
               s <- CELL_META$sector.index; xl <- CELL_META$xlim
               circos.rect(xl[1], 0, xl[2], 1,
                           col = unname(tumor_colors[s]), border = "black")
               circos.text(mean(xl), 0.5, s, facing = "clockwise",
                           niceFacing = TRUE, cex = 0.6, adj = c(0, 0.5))
             })

# Track 2: subtype.
circos.track(factors = vaf_df$sector, x = vaf_df$x, y = rep(1, nrow(vaf_df)),
             ylim = c(0, 1), track.height = 0.07, bg.border = NA,
             panel.fun = function(x, y) {
               ii <- which(as.character(vaf_df$sector) == CELL_META$sector.index)
               for (i in ii) {
                 circos.rect(vaf_df$x[i] - 0.5, 0, vaf_df$x[i] + 0.5, 1,
                             col = unname(subtype_colors[vaf_df$Subtype[i]]), border = NA)
               }
               xl <- CELL_META$xlim
               circos.rect(xl[1], 0, xl[2], 1, col = NA, border = "black")
             })

# Track 3: sex.
circos.track(factors = vaf_df$sector, x = vaf_df$x, y = rep(1, nrow(vaf_df)),
             ylim = c(0, 1), track.height = 0.07, bg.border = NA,
             panel.fun = function(x, y) {
               ii <- which(as.character(vaf_df$sector) == CELL_META$sector.index)
               for (i in ii) {
                 this_color <- unname(sex_colors[vaf_df$Gender[i]])
                 if (is.na(this_color)) this_color <- "#CCCCCC"
                 circos.rect(vaf_df$x[i] - 0.5, 0, vaf_df$x[i] + 0.5, 1,
                             col = this_color, border = NA)
               }
               xl <- CELL_META$xlim
               circos.rect(xl[1], 0, xl[2], 1, col = NA, border = "black")
             })


# Track 4: age at diagnosis.

# Exclude samples with missing or non-finite ages from this track.
age_plot <- vaf_df %>%
  filter(is.finite(x), is.finite(Age))

age_ylim <- range(age_plot$Age)

circos.track(
  factors = age_plot$sector,
  x = age_plot$x,
  y = age_plot$Age,
  ylim = age_ylim,
  track.height = 0.10,
  bg.border = NA,
  panel.fun = function(x, y) {
    
    # Draw age profiles within each tumor type.
    valid <- is.finite(x) & is.finite(y)
    
    if (sum(valid) >= 2) {
      x2 <- x[valid]
      y2 <- y[valid]
      od <- order(x2)
      
      circos.lines(
        x2[od], y2[od],
        type = "l",
        col = "#1B9E77",
        lwd = 1
      )
    }
    
    # Track border
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    
    circos.rect(
      xlim[1], ylim[1], xlim[2], ylim[2],
      border = "black", col = NA, lwd = 1
    )
    
    # Y-axis
    circos.yaxis(
      side = "left",
      at = pretty(age_ylim),
      labels.cex = 0.4,
      tick.length = circlize::convert_x(1, "mm")
    )
  }
)


# Track 5: disease stage at NGS.
circos.track(factors = vaf_df$sector, x = vaf_df$x, y = rep(1, nrow(vaf_df)),
             ylim = c(0, 1), track.height = 0.07, bg.border = NA,
             panel.fun = function(x, y) {
               ii <- which(as.character(vaf_df$sector) == CELL_META$sector.index)
               for (i in ii) {
                 this_color <- unname(stage_colors[vaf_df$Stage[i]])
                 if (is.na(this_color)) this_color <- "#CCCCCC"
                 circos.rect(vaf_df$x[i] - 0.5, 0, vaf_df$x[i] + 0.5, 1,
                             col = this_color, border = NA)
               }
               xl <- CELL_META$xlim
               circos.rect(xl[1], 0, xl[2], 1, col = NA, border = "black")
             })

# Track 6: maximum VAF per sample, log10(max VAF + 1).
vaf_df$log10_vaf <- log10(vaf_df$max_vaf + 1)
circos.track(factors = vaf_df$sector, x = vaf_df$x, y = vaf_df$log10_vaf,
             ylim = range(vaf_df$log10_vaf, na.rm = TRUE),
             track.height = 0.15, bg.border = NA,
             panel.fun = function(x, y) {
               ii <- which(as.character(vaf_df$sector) == CELL_META$sector.index)
               ord <- order(vaf_df$x[ii])
               circos.lines(vaf_df$x[ii][ord], vaf_df$log10_vaf[ii][ord],
                            type = "h", col = "#888888")
               circos.rect(CELL_META$xlim[1], CELL_META$ylim[1],
                           CELL_META$xlim[2], CELL_META$ylim[2], border = "black", col = NA)
             })

circos.clear()
dev.off()



# Figure 1B: mutation-frequency heatmap and pathway proportions ---------
# Multiple variants in the same gene and sample are counted only once.
f1c <- axen_final3 %>% distinct(`SAMPLE ID`, Hugo_Symbol, .keep_all = TRUE)
mutation_freq_df <- f1c %>%
  distinct(Tumor, `SAMPLE ID`, Hugo_Symbol) %>%
  dplyr::count(Tumor, Hugo_Symbol, name = "Mutation_Count") %>%
  left_join(f1c %>% distinct(Tumor, `SAMPLE ID`) %>%
              dplyr::count(Tumor, name = "Sample_Count"), by = "Tumor") %>%
  mutate(Mutation_Frequency = 100 * Mutation_Count / Sample_Count)

mutation_matrix <- mutation_freq_df %>%
  dplyr::select(Tumor, Hugo_Symbol, Mutation_Frequency) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = Mutation_Frequency,
              values_fill = 0) %>%
  as.data.frame()
rownames(mutation_matrix) <- mutation_matrix$Tumor
mutation_matrix$Tumor <- NULL
mutation_matrix <- as.matrix(mutation_matrix)

# The original figure orders cancer types by sample count, then ranks genes
# by their summed mutation frequencies across cancer types.
sample_per_tumor <- f1c %>% distinct(Tumor, `SAMPLE ID`) %>%
  dplyr::count(Tumor, name = "Sample_Count") %>% arrange(desc(Sample_Count))
mutation_matrix <- mutation_matrix[sample_per_tumor$Tumor, , drop = FALSE]
gene_order <- names(sort(colSums(mutation_matrix), decreasing = TRUE))
mutation_matrix_tt <- t(mutation_matrix)[head(gene_order, 20), , drop = FALSE]

TSG_list <- c("TP53", "ATM", "BRCA1", "BRCA2", "APC", "PTEN", "CDKN2A", "FBXW7", "CDH1", "RB1")
OCG_list <- c("KRAS", "PIK3CA", "NOTCH1", "MTOR", "CTNNB1", "ERBB4", "ERBB2", "ALK", "FGFR3", "MYCN", "AR", "RET", "EGFR")
mutation_colors <- c("Missense_Mutation" = "#FA8C6A", "Nonsense_Mutation" = "#F7D652",
                     "In_Frame_Ins" = "#52DEF7", "In_Frame_Del" = "#52DEF7",
                     "Frame_Shift_Ins" = "#854EED", "Frame_Shift_Del" = "#854EED", "Splice_Site" = "#69E056")
mutation_prop <- f1c %>%
  filter(Hugo_Symbol %in% rownames(mutation_matrix_tt)) %>%
  dplyr::count(Hugo_Symbol, Variant_Classification, name = "Count") %>%
  group_by(Hugo_Symbol) %>% mutate(Proportion = Count / sum(Count)) %>%
  ungroup() %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Proportion) %>%
  pivot_wider(names_from = Variant_Classification, values_from = Proportion,
              values_fill = 0) %>% as.data.frame()
rownames(mutation_prop) <- mutation_prop$Hugo_Symbol
mutation_prop$Hugo_Symbol <- NULL
mutation_prop <- as.matrix(mutation_prop[rownames(mutation_matrix_tt), , drop = FALSE])

row_anno <- rowAnnotation(Mutations = anno_barplot(mutation_prop,
  gp = gpar(fill = mutation_colors[colnames(mutation_prop)]), border = TRUE))
top_anno <- HeatmapAnnotation(Sample_Counts = anno_barplot(
  setNames(sample_per_tumor$Sample_Count, sample_per_tumor$Tumor),
  gp = gpar(fill = "gray"), border = TRUE))
gene_group <- setNames(ifelse(rownames(mutation_matrix_tt) %in% TSG_list, "TSG",
                       ifelse(rownames(mutation_matrix_tt) %in% OCG_list, "OG", "Other")),
                       rownames(mutation_matrix_tt))
tsg_col_fun <- colorRamp2(c(0, 30, 60, 100), c("white", "lightblue", "blue", "darkblue"))
og_col_fun <- colorRamp2(c(0, 30, 60, 100), c("white", "pink", "red", "darkred"))
other_col_fun <- colorRamp2(c(0, 30, 60, 100), c("white", "gray80", "gray50", "black"))

ht <- Heatmap(mutation_matrix_tt, name = NULL,
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
              left_annotation = row_anno, top_annotation = top_anno,
              cell_fun = function(j, i, x, y, width, height, fill) {
                gene <- rownames(mutation_matrix_tt)[i]
                col_fun <- switch(gene_group[[gene]], TSG = tsg_col_fun,
                                  OG = og_col_fun, Other = other_col_fun)
                grid.rect(x, y, width, height,
                          gp = gpar(fill = col_fun(mutation_matrix_tt[i, j]),
                                    col = "white", lwd = 2))
              })
pdf(file.path(outdir, "Figure1B_mutation_heatmap.pdf"), width = 12, height = 8)
draw(ht, annotation_legend_list = list(
  Legend(col_fun = tsg_col_fun, title = "TSG Mutation Frequency"),
  Legend(col_fun = og_col_fun, title = "OG Mutation Frequency")))
dev.off()


# The pathway table is an upstream, precomputed analysis input. This code plots
# AXEN_per values; it does not regenerate the gene-to-pathway annotations.
pathway_df <- read.delim(pathway_file, sep = "\t", check.names = FALSE)
stopifnot(all(c("Cancer", "Pathway", "AXEN_per") %in% names(pathway_df)))
ordered_cancer <- colnames(mutation_matrix_tt)  # exactly the heatmap column order
pathway_plot <- pathway_df %>%
  mutate(Cancer = factor(Cancer, levels = ordered_cancer)) %>%
  filter(!is.na(Cancer), !is.na(AXEN_per))

p1b_pathway <- ggplot(pathway_plot, aes(x = Cancer, y = AXEN_per, fill = Pathway)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c(
    "Cell Cycle" = "#F6AAA5", "MYC" = "#C0D9E5", "NOTCH" = "#FCC48F",
    "PI3K" = "#88C47B", "RTK-RAS" = "#D2C4D4", "TGF-Beta" = "#C9C86D",
    "TP53" = "#8A7DA9", "WNT" = "#FFEA85")) +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = "Mutation frequency (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, "Figure1B_pathways.pdf"), p1b_pathway, width = 12, height = 3)



# Figure 1C: sample-level mutation count -------------------------------
# Stage information is not used for filtering this panel.
plot_count <- mutation_counts %>%
  left_join(axen_final3 %>% distinct(Tumor_Sample_Barcode, Tumor),
            by = "Tumor_Sample_Barcode")
p1c <- ggplot(plot_count,
              aes(x = fct_reorder(Tumor, count), y = log2(count), fill = Tumor)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.3) +
  scale_fill_manual(values = tumor_colors) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_classic() +
  labs(x = NULL, y = "log2(Mutation count)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(outdir, "Figure1C_mutation_count.pdf"), p1c, width = 12, height = 4)



# Figure 1D: mut-H versus mut-L by tumor type ---------------------------
# Denominator: unique ctDNA samples within each tumor type.
df_ratio <- axen_final3 %>%
  distinct(Tumor, Tumor_Sample_Barcode, count_edit) %>%
  dplyr::count(Tumor, count_edit) %>%
  group_by(Tumor) %>% mutate(ratio = n / sum(n)) %>% ungroup()
mutH_order <- df_ratio %>% group_by(Tumor) %>%
  summarise(mutH_ratio = sum(ratio[count_edit == "mut-H"], na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(mutH_ratio)) %>% pull(Tumor)
df_ratio$Tumor <- factor(df_ratio$Tumor, levels = mutH_order)
p1d <- ggplot(df_ratio, aes(x = Tumor, y = ratio, fill = count_edit)) +
  geom_col(position = "fill") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = NULL, y = "Proportion") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(outdir, "Figure1D_mutation_count_groups.pdf"), p1d, width = 12, height = 4)



# Figure 1E: MATH score (minimum 3 nonsynonymous mutations) -------------
# Based on maftools::math.score, with the original default of 5 changed to 3.
# Avoid the original interactive command trace(math.score, edit = TRUE).
math_score_min3 <- function (maf, vafCol = NULL, sampleName = NULL, vafCutOff = 0.075) 
{
  if (is.null(sampleName)) {
    sampleName = as.character(getSampleSummary(maf)[, Tumor_Sample_Barcode])
  }
  dat = subsetMaf(maf, includeSyn = FALSE, tsb = sampleName, 
                  mafObj = FALSE)
  if (!"t_vaf" %in% colnames(dat)) {
    if (is.null(vafCol)) {
      if (all(c("t_ref_count", "t_alt_count") %in% colnames(dat))) {
        message("t_vaf field is missing, but found t_ref_count & t_alt_count columns. Estimating vaf..")
        dat[, `:=`(t_ref_count, as.numeric(as.character(t_ref_count)))]
        dat[, `:=`(t_alt_count, as.numeric(as.character(t_alt_count)))]
        dat[, `:=`(t_vaf, t_alt_count/(t_ref_count + 
                                         t_alt_count))]
      }
      else {
        print(colnames(maf))
        stop("t_vaf field is missing. Use vafCol to manually specify vaf column name.")
      }
    }
    else {
      colnames(dat)[which(colnames(dat) == vafCol)] = "t_vaf"
      dat[, `:=`(t_vaf, as.numeric(as.character(t_vaf)))]
    }
  }
  if (max(dat[, t_vaf], na.rm = TRUE) > 1) {
    dat[, `:=`(t_vaf, as.numeric(as.character(t_vaf))/100)]
  }
  dat = dat[!t_vaf < vafCutOff]
  dat = dat[, .(Tumor_Sample_Barcode, Hugo_Symbol, t_vaf)]
  dat.spl = split(dat, as.factor(as.character(dat$Tumor_Sample_Barcode)))
  math.dt = lapply(X = dat.spl, function(pat) {
    vaf = pat$t_vaf
    vaf = vaf[!is.na(vaf)]
    if (length(vaf) < 3) {
      message(paste("Not enough mutations in", unique(pat[, Tumor_Sample_Barcode]), ". Skipping.."), sep = "")
    }
    else {
      abs.med.dev = abs(vaf - median(vaf))
      pat.mad = median(abs.med.dev) * 100
      pat.math = pat.mad * 1.4826/median(vaf)
      d = data.table::data.table(Tumor_Sample_Barcode = as.character(unique(pat[, Tumor_Sample_Barcode])), MedianAbsoluteDeviation = pat.mad, 
                                 MATH = pat.math)
      return(d)
    }
  })
  math.dt = data.table::rbindlist(math.dt, use.names = TRUE, fill = TRUE)
  return(math.dt)
}

Axen_math <- math_score_min3(maf_obj, vafCutOff = 0) %>%
  left_join(axen_final3 %>% distinct(Tumor_Sample_Barcode, Tumor), by = "Tumor_Sample_Barcode")
message("Samples with MATH scores: ", nrow(Axen_math))
p1e <- ggplot(Axen_math, aes(x = fct_reorder(Tumor, MATH), y = MATH, fill = Tumor)) +
  geom_boxplot() + scale_fill_manual(values = tumor_colors) +
  labs(x = NULL, y = "MATH") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(outdir, "Figure1E_MATH.pdf"), p1e, width = 12, height = 4)



# Figure 1F: OncoKB v5.0 annotated alterations -------------------------
# Original analytical strategy: 87 genes; restricted LUAD/LUSC/LCLC/NSCLCPD
# for lung cancer and COADREAD for intestinal cancer; patient-level denominator.
onco <- readxl::read_xlsx(oncokb_file, sheet = "v5.0") %>%
  filter(is.na(annotation), !is.na(Tumor), Tumor != "Liquid") %>%
  mutate(Tumor = if_else(Tumor == "Skin cancer", "Melanoma", Tumor)) %>%
  distinct(Level, Gene, Alterations, `Cancer Types`, Tumor, category, .keep_all = TRUE)

genes <- read.delim(path = panel_file)
KCA = data.frame(KCA = matrix(genes$`K&C&A_87`, ncol = 1, nrow = 87)) # 87개 genes
onco <- onco %>% filter(Gene %in% KCA$KCA)
axen_final4 <- axen_final3 %>% filter(Hugo_Symbol %in% KCA$KCA)
onco_mut <- onco %>% filter(Alterations != "Wildtype")

onco_wt <- onco %>% filter(Gene %in% c("KRAS", "NRAS"), Alterations == "Wildtype")

fig1g <- axen_final4 %>%
  filter((Tumor == "Lung cancer" & Tumor_subtype_JW %in% c("LCLC", "LUAD", "LUSC", "NSCLCPD")) |
           (Tumor == "Intestinal cancer" & Tumor_subtype_JW == "COADREAD") |
           !(Tumor %in% c("Lung cancer", "Intestinal cancer"))) # 3125 mut
fig1g$HGVSp_Short_edit <- sub("^p\\.", "", fig1g$HGVSp_Short)

tumor_totals_axen <- fig1g %>% distinct(Tumor, `SUBJECT ID`) %>% group_by(Tumor) %>% summarise(AXEN_total = n(), .groups = "drop")

### 1) pan + all
onco_pan_onco <- onco_mut %>% filter(Tumor == "Pan-cancer", category == "all_mutation")

onco_pan_onco1 <- fig1g %>% filter(Hugo_Symbol %in% onco_pan_onco$Gene) %>% group_by(Tumor, `SUBJECT ID`, Hugo_Symbol) %>%
  summarise(n = n(), .groups = "drop") %>% group_by(Tumor, Hugo_Symbol) %>% summarise(AXEN_mut = n(), .groups = "drop") %>%
  left_join(tumor_totals_axen, by = "Tumor") %>%
  mutate(AXEN_per = round(AXEN_mut / AXEN_total * 100, 2), HGVSp_Short_edit = "all")

onco_pan_onco11 <- onco_pan_onco1 %>%
  left_join(onco %>% filter(Tumor == "Pan-cancer", category == "all_mutation") %>% dplyr::select(Gene, Level),
            by = c("Hugo_Symbol" = "Gene")) %>% mutate(Level = as.character(Level)) %>%
  arrange(Tumor, Hugo_Symbol, HGVSp_Short_edit, factor(Level, levels = c("1", "2", "3", "4", NA), ordered = TRUE)) %>%
  distinct(Tumor, Hugo_Symbol, HGVSp_Short_edit, .keep_all = TRUE)

### 2) pan + specific
onco_pan_spe <- onco_mut %>% filter(Tumor == "Pan-cancer", category != "all_mutation") %>%
  dplyr::select(Gene, Alterations, Level) %>% rename(Pan_Level = Level)

onco_pan_spe1 <- fig1g %>% filter(Hugo_Symbol %in% onco_pan_spe$Gene, HGVSp_Short_edit %in% onco_pan_spe$Alterations) %>%
  group_by(Tumor, `SUBJECT ID`, Hugo_Symbol, HGVSp_Short_edit) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(Tumor, Hugo_Symbol, HGVSp_Short_edit) %>% summarise(AXEN_mut = n(), .groups = "drop") %>%
  left_join(tumor_totals_axen, by = "Tumor") %>% mutate(AXEN_per = round(AXEN_mut / AXEN_total * 100, 2))

onco_pan_spe11 <- onco_pan_spe1 %>%
  left_join(onco %>% filter(category != "all_mutation") %>% dplyr::select(Tumor, Gene, Alterations, Level),
            by = c("Tumor", "Hugo_Symbol" = "Gene", "HGVSp_Short_edit" = "Alterations")) %>%
  left_join(onco_pan_spe, by = c("Hugo_Symbol" = "Gene", "HGVSp_Short_edit" = "Alterations")) %>%
  mutate(Level = coalesce(Level, Pan_Level), Level = as.character(Level)) %>%
  dplyr::select(-Pan_Level) %>%
  arrange(Tumor, Hugo_Symbol, HGVSp_Short_edit, factor(Level, levels = c("1", "2", "3", "4", NA), ordered = TRUE)) %>%
  distinct(Tumor, Hugo_Symbol, HGVSp_Short_edit, .keep_all = TRUE)

### 3) specific cancer + all mut
onco_spe_onco <- onco_mut %>%
  filter(Tumor != "Pan-cancer", category == "all_mutation")

onco_spe_onco1 <- fig1g %>%
  filter(Hugo_Symbol %in% onco_spe_onco$Gene) %>%
  semi_join(
    onco_spe_onco,
    by = c("Hugo_Symbol" = "Gene", "Tumor")
  ) %>%
  group_by(Tumor, `SUBJECT ID`, Hugo_Symbol) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Tumor, Hugo_Symbol) %>%
  summarise(AXEN_mut = n(), .groups = "drop") %>%
  left_join(tumor_totals_axen, by = "Tumor") %>%
  mutate(
    AXEN_per = round(AXEN_mut / AXEN_total * 100, 2),
    HGVSp_Short_edit = "all"
  )

onco_spe_onco11 <- onco_spe_onco1 %>%
  left_join(
    onco %>%
      filter(Tumor != "Pan-cancer", category == "all_mutation") %>%
      dplyr::select(Tumor, Gene, Level),
    by = c("Tumor", "Hugo_Symbol" = "Gene")
  ) %>%
  mutate(Level = as.character(Level)) %>%
  arrange(
    Tumor,
    Hugo_Symbol,
    HGVSp_Short_edit,
    factor(Level, levels = c("1", "2", "3", "4", NA), ordered = TRUE)
  ) %>%
  distinct(Tumor, Hugo_Symbol, HGVSp_Short_edit, .keep_all = TRUE)

### 4) specific cancer + specific mut
onco_spe_spe <- onco_mut %>%
  filter(Tumor != "Pan-cancer", category != "all_mutation")

onco_spe_spe1 <- fig1g %>%
  filter(
    Hugo_Symbol %in% onco_spe_spe$Gene,
    HGVSp_Short_edit %in% onco_spe_spe$Alterations
  ) %>%
  semi_join(
    onco_spe_spe,
    by = c(
      "Hugo_Symbol" = "Gene",
      "HGVSp_Short_edit" = "Alterations",
      "Tumor"
    )
  ) %>%
  group_by(Tumor, `SUBJECT ID`, Hugo_Symbol, HGVSp_Short_edit) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Tumor, Hugo_Symbol, HGVSp_Short_edit) %>%
  summarise(AXEN_mut = n(), .groups = "drop") %>%
  left_join(tumor_totals_axen, by = "Tumor") %>%
  mutate(AXEN_per = round(AXEN_mut / AXEN_total * 100, 2))

onco_spe_spe11 <- onco_spe_spe1 %>%
  left_join(
    onco %>%
      filter(category != "all_mutation") %>%
      dplyr::select(Tumor, Gene, Alterations, Level),
    by = c(
      "Tumor",
      "Hugo_Symbol" = "Gene",
      "HGVSp_Short_edit" = "Alterations"
    )
  ) %>%
  left_join(
    onco_pan_spe,
    by = c(
      "Hugo_Symbol" = "Gene",
      "HGVSp_Short_edit" = "Alterations"
    )
  ) %>%
  mutate(
    Level = coalesce(Level, Pan_Level),
    Level = as.character(Level)
  ) %>%
  dplyr::select(-Pan_Level) %>%
  arrange(
    Tumor,
    Hugo_Symbol,
    HGVSp_Short_edit,
    factor(Level, levels = c("1", "2", "3", "4", NA), ordered = TRUE)
  ) %>%
  distinct(Tumor, Hugo_Symbol, HGVSp_Short_edit, .keep_all = TRUE)


### K/RAS WT
# 전체 환자 universe
patient_universe <- fig1g %>%
  distinct(Tumor, `SUBJECT ID`)

# KRAS/NRAS mutation이 있는 환자
ras_mut_patients <- fig1g %>%
  filter(Hugo_Symbol %in% c("KRAS", "NRAS")) %>%
  distinct(Tumor, `SUBJECT ID`, Hugo_Symbol)

# onco_wt가 Pan-cancer annotation이면 모든 tumor에 적용
onco_wt_pan <- onco_wt %>%
  filter(Tumor == "Pan-cancer") %>%
  dplyr::select(Gene, Level) %>%
  distinct() %>%
  tidyr::crossing(Tumor = unique(patient_universe$Tumor)) %>%
  dplyr::select(Tumor, Gene, Level)

# onco_wt가 specific cancer annotation이면 해당 tumor에만 적용
onco_wt_spe <- onco_wt %>%
  filter(Tumor != "Pan-cancer") %>%
  dplyr::select(Tumor, Gene, Level) %>%
  distinct()

# 최종 WT annotation table
onco_wt_map <- bind_rows(onco_wt_pan, onco_wt_spe) %>%
  distinct(Tumor, Gene, .keep_all = TRUE)

wt_total <- onco_wt_map %>%
  left_join(patient_universe, by = "Tumor") %>%
  anti_join(
    ras_mut_patients,
    by = c("Tumor", "SUBJECT ID", "Gene" = "Hugo_Symbol")
  ) %>%
  group_by(Tumor, Gene, Level) %>%
  summarise(AXEN_mut = n_distinct(`SUBJECT ID`), .groups = "drop") %>%
  left_join(tumor_totals_axen, by = "Tumor") %>%
  mutate(
    AXEN_per = round(AXEN_mut / AXEN_total * 100, 2),
    Alt = paste0(Gene, " WT"),
    Level = as.character(Level)
  ) %>%
  dplyr::select(Alt, Tumor, AXEN_per, Level)

level_order <- c("1", "2", "3", "4", NA)

total2 <- bind_rows(
  onco_pan_onco11,
  onco_pan_spe11,
  onco_spe_onco11,
  onco_spe_spe11
) %>%
  mutate(
    Alt = paste0(Hugo_Symbol, " ", HGVSp_Short_edit),
    Tumor = as.character(Tumor),
    Level = as.character(Level),
    Level = factor(Level, levels = level_order)
  ) %>%
  dplyr::select(Alt, Tumor, AXEN_per, Level) %>%
  arrange(Alt, Tumor, Level) %>%
  distinct(Alt, Tumor, .keep_all = TRUE)

total2_fixed <- total2 %>%
  mutate(
    Tumor = as.character(Tumor),
    AXEN_per = as.numeric(AXEN_per),
    Level = as.character(Level),
    Level = ifelse(is.na(Level), "NA", Level)
  )

total2_fixed_wt <- bind_rows(
  total2_fixed,
  wt_total
) %>%
  mutate(
    Tumor = as.character(Tumor),
    AXEN_per = as.numeric(AXEN_per),
    Level = as.character(Level),
    Level = ifelse(is.na(Level), "NA", Level)
  ) %>%
  arrange(
    Alt,
    Tumor,
    factor(Level, levels = c("1", "2", "3", "4", "NA"), ordered = TRUE)
  ) %>%
  distinct(Alt, Tumor, .keep_all = TRUE)

group1 <- c(
  "EGFR E746_A750del",
  "EGFR L747_S752del",
  "EGFR L747_T751del",
  "EGFR L747_P753delinsS"
)

group2 <- c("EGFR A767_V769dup")

base_df <- total2_fixed_wt %>%
  filter(!Alt %in% c(group1, group2))

egfr_19 <- total2_fixed_wt %>%
  filter(Alt %in% group1) %>%
  group_by(Tumor) %>%
  summarise(AXEN_per = sum(AXEN_per, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Alt = "EGFR 19 in-frame del",
    Level = "1"
  )

egfr_20 <- total2_fixed_wt %>%
  filter(Alt %in% group2) %>%
  group_by(Tumor) %>%
  summarise(AXEN_per = sum(AXEN_per, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Alt = "EGFR 20 in-frame ins",
    Level = "1"
  )

total2_grouped <- bind_rows(base_df, egfr_19, egfr_20) %>%
  group_by(Alt, Tumor) %>%
  summarise(
    AXEN_per = sum(AXEN_per, na.rm = TRUE),
    Level = first(Level),
    .groups = "drop"
  ) %>%
  mutate(Level = ifelse(is.na(Level), "NA", Level))
total2_grouped <- total2_grouped %>%
  mutate(
    Alt = sub(" all$", "", Alt)
  )

level_colors <- c(
  "1" = "#44B549",
  "2" = "#3768B2",
  "3" = "#6E2E91",
  "4" = "#626164",
  "NA" = "lightgray"
)

tumor_order <- c(
  "Lung cancer",
  "Breast cancer",
  "Intestinal cancer",
  "Gallbladder cancer",
  "Pancreatic cancer",
  "Gastric cancer",
  "Prostate cancer",
  "Head and neck cancer",
  "Sarcoma",
  "Bladder cancer",
  "Liver cancer",
  "Renal cell cancer",
  "Esophageal cancer",
  "Melanoma",
  "Neuroendocrine tumor",
  "Ovarian cancer",
  "Thymoma"
)

tumor_labels <- c("NSCLC", "BC", "CRC", "HBC", "PC", "GC", "PCA", "HNSC", "SARC", "BLCA", "LIHC", "RCC", "EC", "SKCM", "NET", "OC", "THYM")
names(tumor_labels) <- tumor_order

alt_order <- total2_grouped %>% mutate(Level = factor(Level, levels = c("1", "2", "3", "4", "NA"), ordered = TRUE)) %>%
  arrange(Level, Alt) %>% distinct(Alt) %>% pull(Alt)

p1f <- ggplot(total2_grouped, aes(x = factor(Alt, levels = alt_order), y = factor(Tumor, levels = rev(tumor_order)))) +
  geom_point(shape = 21, aes(size = AXEN_per, fill = Level), alpha = 0.9, color = "black", stroke = 0.2) +
  scale_fill_manual(values = level_colors) +
  scale_size_continuous(range = c(2, 7), breaks = c(10, 20, 30, 40, 50), name = "Mutation frequency (%)") +
  scale_y_discrete(labels = tumor_labels) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title = element_blank()) +
  labs(fill = "OncoKB Level")

ggsave(file.path(outdir, "Figure1F_OncoKB_v5.pdf"), p1f, width = 16, height = 7)

# Save the plot data for verification (not the patient-level underlying MAF).
# Do NOT commit result tables derived from controlled-access patient data
# without an explicit data-release review.
# Code sharing and data sharing are subject to different permissions.
