# =============================================================================
# KEGG Module Differential Abundance Analysis — MaAsLin2
# Colleen Ahern | Last updated: 05/13/2026
#
# Three comparisons:
#   1. CA_S3_AS  vs. CA_AS
#   2. PHA_G1_AS vs. PHA_AS
#   3. PHA_G1_AS t2 vs. PHA_G1_AS t1
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(Maaslin2)
library(pheatmap)
library(ggplotify)
library(cowplot)
library(KEGGREST)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    <- "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics"
MT_DIR  <- file.path(BASE, "MT_analysis/MT_analysis")
OUT_DIR <- file.path(BASE, "MT_grouped_reads/maaslin2")
FIG_DIR <- "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures"
META    <- file.path(BASE, "bp_metadata_01132025.txt")
DATE    <- "05132026"

# =============================================================================
# 1.  Load and merge data
# =============================================================================

# ── Taxa annotations ──────────────────────────────────────────────────────────
taxa_gtdb <- read_tsv(file.path(MT_DIR, "processed_tables/tax_gtdb_metassembly.tsv")) %>%
  column_to_rownames("contig") %>%
  unite("taxa", domain, phylum, class, order, family, genus, species, sep = " ", remove = FALSE) %>%
  rownames_to_column("Contig")

# ── Contig-gene mapping ───────────────────────────────────────────────────────
contig_gene <- read.delim(
  file.path(MT_DIR, "input/bowtie2_metassembly.stranded2.counts.txt"),
  header = FALSE
)
contig_gene <- contig_gene[-1, ]
names(contig_gene) <- contig_gene[1, ]
contig_gene <- contig_gene[-1, ]
names(contig_gene) <- gsub("mapping/bowtie.meta_t\\.|.sorted.bam", "", names(contig_gene))
contig_gene <- contig_gene %>% rename(Contig = Chr)

contig_gene <- merge(contig_gene, taxa_gtdb, all.x = TRUE)

# ── Count matrix ──────────────────────────────────────────────────────────────
bpcoldata <- read.delim(META) %>%
  column_to_rownames("name")

bpcts <- contig_gene %>%
  column_to_rownames("Geneid") %>%
  select(-c(1:5, 26:33))

# align sample order
bpcoldata2 <- bpcoldata[rownames(bpcoldata) %in% colnames(bpcts), ]
bpcts      <- bpcts[, rownames(bpcoldata2)]
bpcts$GeneID <- rownames(bpcts)

# ── Gene and taxa annotations ─────────────────────────────────────────────────
gene_anno <- read_tsv(file.path(MT_DIR, "processed_tables/gene_annotation_metassembly.tsv"))
tax_anno  <- contig_gene[, c(2, 27:34)] %>%
  rename(GeneID = Geneid) %>%
  as.data.frame()

bpcts_geneanno <- merge(bpcts,      gene_anno, by = "GeneID", all = TRUE)
bpcts_geneanno <- merge(bpcts_geneanno, tax_anno,  by = "GeneID", all = TRUE)

# =============================================================================
# 2.  Aggregate counts by KEGG module
# =============================================================================

aa <- bpcts_geneanno %>%
  select(GeneID, everything()) %>%
  mutate(KEGG_Module = str_split(KEGG_Module, ",")) %>%
  unnest_wider(where(is.list), names_sep = "") %>%
  as.data.frame()

aa <- aa[, c(2:68, 1)]
aa$KEGG_Module1[is.na(aa$KEGG_Module1)] <- "-"

# ── Fast vectorized count aggregation by KEGG Module ──────────────────────
keggmod_raw <- aa %>%
  pivot_longer(cols = starts_with("KEGG_Module"),
               values_to = "KEGG_Module",
               names_to  = NULL) %>%
  filter(!is.na(KEGG_Module), KEGG_Module != "-") %>%
  group_by(KEGG_Module) %>%
  summarise(across(1:20, ~ sum(as.numeric(.), na.rm = TRUE))) %>%
  filter(KEGG_Module != "-")

# reshape aa to long format with one row per gene-module combination
sample_cols <- names(aa)[grepl("MT_PLANC", names(aa))]
sample_cols
aa_long <- aa %>%
  mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.)))) %>%
  pivot_longer(cols = starts_with("KEGG_Module"),
               values_to = "KEGG_Module",
               names_to  = NULL) %>%
  filter(!is.na(KEGG_Module), KEGG_Module != "-") %>%
  mutate(total_count = rowSums(across(all_of(sample_cols)), na.rm = TRUE)) %>%
  filter(total_count > 0)

# aggregate unique taxa per module
taxa_per_module <- aa_long %>%
  group_by(KEGG_Module) %>%
  summarise(taxa = paste(unique(na.omit(taxa)), collapse = "; "))

# merge into keggmod_raw
keggmod_raw <- left_join(keggmod_raw, taxa_per_module, by = "KEGG_Module")

write_tsv(keggmod_raw,
          file.path(dirname(OUT_DIR), paste0("KEGGmod_raw_taxaanno_", DATE, ".tsv")))

# =============================================================================
# 3.  Prepare for MaAsLin2
# =============================================================================

keggmod_rawrt <- keggmod_raw %>%
  column_to_rownames("KEGG_Module") %>%
  select(all_of(sample_cols))  # keeps only the 20 count columns, drops taxa

taxa_annot <- keggmod_raw %>%
  select(KEGG_Module, taxa)

# drop outlier samples identified by PCA
drop_samples <- c("I6_CA14_1_MT_PLANC", "I8_CA58_1_MT_PLANC")
keggmod_rawrt <- keggmod_rawrt[, !names(keggmod_rawrt) %in% drop_samples]

bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(keggmod_rawrt), ] %>%
  mutate(
    group         = gsub(" + ", "_", group,         fixed = TRUE),
    condition_cba = gsub(" + ", "_", condition_cba, fixed = TRUE),
    condition_cba = gsub(" ",   "_", condition_cba),
    name          = rownames(.)
  )

# rows = samples, cols = KEGG modules
keggmod_t <- as.data.frame(t(keggmod_rawrt))
stopifnot(all(rownames(keggmod_t) == rownames(bpcoldata3)))

# ── MaAsLin2 helper ───────────────────────────────────────────────────────────
run_maaslin <- function(data, meta, output, fixed_effect, reference) {
  Maaslin2(
    input_data      = data,
    input_metadata  = meta,
    output          = output,
    fixed_effects   = fixed_effect,
    reference       = reference,
    normalization   = "TSS",
    transform       = "LOG",
    analysis_method = "LM",
    min_prevalence  = 0.1,
    min_abundance   = 0.0,
    cores           = 1
  )
}

# =============================================================================
# 4.  MaAsLin2 — three comparisons
# =============================================================================

# ── 1. CA_S3_AS vs. CA_AS ─────────────────────────────────────────────────────
meta1 <- bpcoldata3 %>% filter(group %in% c("CA_S3_AS", "CA_AS"))
fit1  <- run_maaslin(
  data         = keggmod_t[rownames(meta1), ],
  meta         = meta1,
  output       = file.path(OUT_DIR, paste0("KEGG_module_maaslin2_CA_", DATE)),
  fixed_effect = "group",
  reference    = "group,CA_AS"
)

# ── 2. PHA_G1_AS vs. PHA_AS ───────────────────────────────────────────────────
meta2 <- bpcoldata3 %>% filter(group %in% c("PHA_G1_AS", "PHA_AS"))
fit2  <- run_maaslin(
  data         = keggmod_t[rownames(meta2), ],
  meta         = meta2,
  output       = file.path(OUT_DIR, paste0("KEGG_module_maaslin2_PHA_", DATE)),
  fixed_effect = "group",
  reference    = "group,PHA_AS"
)

# ── 3. PHA_G1_AS t2 vs. PHA_G1_AS t1 ─────────────────────────────────────────
meta3 <- bpcoldata3 %>% filter(condition_cba %in% c("PHA_G1_AS_t2", "PHA_G1_AS_t1"))
fit3  <- run_maaslin(
  data         = keggmod_t[rownames(meta3), ],
  meta         = meta3,
  output       = file.path(OUT_DIR, paste0("KEGG_module_maaslin2_PHAG1ASt2t1_", DATE)),
  fixed_effect = "condition_cba",
  reference    = "condition_cba,PHA_G1_AS_t1"
)

# =============================================================================
# 5.  Combine results
# =============================================================================

combine_results <- function(fit, feature_col = "feature") {
  fit$results %>% select(feature, coef, qval)
}

res_CA   <- combine_results(fit1)
res_PHA  <- combine_results(fit2)
res_PHAb <- combine_results(fit3)

all_modules <- union(res_CA$feature, union(res_PHA$feature, res_PHAb$feature))

combined_mat <- data.frame(
  KEGG_Module = all_modules,
  lfc_CA      = setNames(res_CA$coef,   res_CA$feature)[all_modules],
  q_CA        = setNames(res_CA$qval,   res_CA$feature)[all_modules],
  lfc_PHA     = setNames(res_PHA$coef,  res_PHA$feature)[all_modules],
  q_PHA       = setNames(res_PHA$qval,  res_PHA$feature)[all_modules],
  lfc_PHAb    = setNames(res_PHAb$coef, res_PHAb$feature)[all_modules],
  q_PHAb      = setNames(res_PHAb$qval, res_PHAb$feature)[all_modules],
  row.names   = all_modules
) %>%
  mutate(across(starts_with("lfc"), ~ replace_na(., 0)),
         across(starts_with("q"),   ~ replace_na(., 1)))

# =============================================================================
# 6.  Fetch KEGG module annotations
# =============================================================================

# Only fetch if not already saved
kegg_anno_path <- file.path(BASE, paste0("KEGG_modules_list_with_categories_", DATE, ".tsv"))

if (!file.exists(kegg_anno_path)) {
  module_ids <- names(keggList("module"))
  chunks     <- split(module_ids, ceiling(seq_along(module_ids) / 10))
  results    <- list()
  
  for (i in seq_along(chunks)) {
    cat("Fetching chunk", i, "of", length(chunks), "\n")
    results <- c(results, keggGet(chunks[[i]]))
    Sys.sleep(0.3)
  }
  
  module_table <- do.call(rbind, lapply(results, function(x) {
    class_split <- strsplit(x$CLASS, "; ")[[1]]
    data.frame(
      KEGG_Module = x$ENTRY,
      module_name = x$NAME,
      category1   = ifelse(length(class_split) >= 1, class_split[1], NA),
      category2   = ifelse(length(class_split) >= 2, class_split[2], NA),
      category3   = ifelse(length(class_split) >= 3, class_split[3], NA),
      stringsAsFactors = FALSE
    )
  }))
  
  write_tsv(module_table, kegg_anno_path)
} else {
  module_table <- read_tsv(kegg_anno_path)
}

# Merge annotations
combined_mat_anno <- merge(combined_mat, module_table, by = "KEGG_Module", all.x = TRUE)
write_tsv(combined_mat_anno,
          file.path(OUT_DIR, paste0("KEGG_module_maaslin2_CAPHAPHAb_coeffpadjanno_", DATE, ".tsv")))

# =============================================================================
# 7.  Heatmaps
# =============================================================================

# ── Helper: build plot and annotation matrices ────────────────────────────────
build_heatmap_mats <- function(df, lfc_cols, q_cols, col_labels) {
  plot_mat <- as.data.frame(df[, lfc_cols])
  padj_mat <- as.data.frame(df[, q_cols])
  rownames(plot_mat) <- rownames(padj_mat) <- paste0(df$KEGG_Module, " - ", df$module_name)
  
  plot_masked <- plot_mat
  plot_masked[padj_mat > 0.05] <- 0
  colnames(plot_masked) <- col_labels
  plot_masked
}

# ── Color palette for KEGG categories ────────────────────────────────────────
cat2_palette <- c(
  "Carbohydrate metabolism"                     = "#009E73",  # teal green
  "Amino acid metabolism"                       = "#0072B2",  # blue
  "Nucleotide metabolism"                       = "#E69F00",  # yellow orange
  "Lipid metabolism"                            = "#CC79A7",  # pink purple
  "Glycan metabolism"                           = "#56B4E9",  # sky blue
  "Energy metabolism"                           = "#D55E00",  # vermillion
  "Metabolism of cofactors and vitamins"        = "#F0E442",  # yellow
  "Xenobiotics biodegradation"                  = "#648FFF",  # periwinkle
  "Biosynthesis of other secondary metabolites" = "#FE6100",  # orange
  "Biosynthesis of terpenoids and polyketides"  = "#785EF0",  # violet
  "Gene set"                                    = "#808080"   # grey
)

make_heatmap <- function(df, lfc_cols, q_cols, col_labels, cellwidth = 120) {
  plot_masked <- build_heatmap_mats(df, lfc_cols, q_cols, col_labels)
  
  cat2_vals   <- unique(na.omit(df$category2))
  cat2_colors <- cat2_palette[names(cat2_palette) %in% cat2_vals]
  missing     <- setdiff(cat2_vals, names(cat2_palette))
  if (length(missing) > 0)
    cat2_colors <- c(cat2_colors, setNames(rainbow(length(missing)), missing))
  
  row_annot <- data.frame(
    "KEGG Category" = as.character(df$category2),
    row.names       = paste0(df$KEGG_Module, " - ", df$module_name),
    check.names     = FALSE
  )
  
  pheatmap(
    plot_masked,
    cluster_cols         = FALSE,
    cluster_rows         = TRUE,
    color                = colorRampPalette(c("blue", "black", "red"))(100),
    breaks               = seq(-2, 2, length.out = 101),
    border_color         = "grey80",
    annotation_row       = row_annot,
    annotation_colors    = list("KEGG Category" = cat2_colors),
    annotation_names_row = FALSE,
    fontsize             = 12,
    fontsize_row         = 10,
    angle_col            = 0,
    cellwidth            = cellwidth,
    cellheight           = 13,
    treeheight_row       = 50
  )
}

# ── heatmap ────────────────────────────────────────────────────────────────
df_all <- combined_mat_anno %>%
  filter(q_CA < 0.05 | q_PHA < 0.05 | q_PHAb < 0.05, !is.na(module_name))

modph_all <- make_heatmap(
  df         = df_all,
  lfc_cols   = c("lfc_CA", "lfc_PHA", "lfc_PHAb"),
  q_cols     = c("q_CA",   "q_PHA",   "q_PHAb"),
  col_labels = c("CA + S3 + AS vs.\nCA + AS",
                 "PHA + G1 + AS vs.\nPHA + AS",
                 "PHA + G1 + AS t2 vs.\nPHA + G1 + AS t1")
)

# =============================================================================
# 8.  Save figures
# =============================================================================
# Convert pheatmaps to ggplot objects

cog_gg <- as.ggplot(cogph) +
  theme(plot.margin = margin(t = 19, b = 0, l = -575, r = 0))

mod_gg <- as.ggplot(modph_all) +
  theme(plot.margin = margin(t = 12, b = 10, l = 16, r = 0))

# Stack and label
combo_gg <- plot_grid(cog_gg, mod_gg,
                      labels = c("A)", "B)"),
                      ncol = 1,
                      rel_heights = c(1, 2.6))  # adjust ratio 

combo_gg

ggsave(plot     = combo_gg,
       filename = file.path(FIG_DIR, paste0("KEGGmodcog_maaslin2_CAPHAPHAb_heatmap_", DATE, ".pdf")),
       width = 20.5, height = 18, units = "in", dpi = 300, bg = "white")

message("Done. Outputs written to: ", OUT_DIR)