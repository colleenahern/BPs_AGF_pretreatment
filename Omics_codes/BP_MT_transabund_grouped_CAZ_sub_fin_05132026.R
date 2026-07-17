# =============================================================================
# CAZyme Differential Abundance Analysis — MaAsLin2
# Colleen Ahern | Last updated: 05/13/2026
#
# Two comparisons:
#   1. CA_S3_AS  vs. CA_AS
#   2. PHA_G1_AS vs. PHA_AS
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(DRnaSeq)
library(Maaslin2)
library(ggrepel)
library(ggpubr)
library(pheatmap)

# ── Constants ─────────────────────────────────────────────────────────────────
DATE    <- "05132026"
BASE    <- "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics"
OUT_DIR <- file.path(BASE, "MT_grouped_reads/maaslin2")
FIG_DIR <- "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots"

# ── Load taxonomy and contig/gene table ───────────────────────────────────────
taxa_gtdb <- read_tsv(file.path(BASE, "MT_analysis/MT_analysis/processed_tables/tax_gtdb_metassembly.tsv"))

contig_gene <- read.delim(
  file.path(BASE, "MT_analysis/MT_analysis/input/bowtie2_metassembly.stranded2.counts.txt"),
  header = FALSE
)
contig_gene <- contig_gene[-1, ]
colnames(contig_gene) <- unlist(contig_gene[1, ])
contig_gene <- contig_gene[-1, ]
names(contig_gene) <- gsub("mapping/bowtie.meta_t\\.|.sorted\\.bam", "", names(contig_gene))

# ── Merge contig/gene with taxonomy ───────────────────────────────────────────
taxa_gtdb2 <- taxa_gtdb |>
  column_to_rownames("contig") |>
  unite("taxa", domain, phylum, class, order, family, genus, species, sep = " ", remove = FALSE) |>
  rownames_to_column("Contig")

contig_gene2 <- contig_gene |>
  rename(Contig = Chr) |>
  merge(taxa_gtdb2, all.x = TRUE)

# ── Load metadata and gene annotations ────────────────────────────────────────
bpcoldata <- read.delim(file.path(BASE, "bp_metadata_01132025.txt")) |>
  column_to_rownames("name")

gene_anno <- read_tsv(file.path(BASE, "MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv"))

# ── Build count/annotation tables ─────────────────────────────────────────────
tax_anno <- contig_gene2[, c(2, 27:34)] |>
  rename(GeneID = Geneid) |>
  as.data.frame()

bpcts <- contig_gene2 |>
  column_to_rownames("Geneid") |>
  select(-c(1:5, 26:33))

bpcoldata2 <- bpcoldata[rownames(bpcoldata) %in% colnames(bpcts), ]
bpcts      <- bpcts[, rownames(bpcoldata2)]
stopifnot(all(rownames(bpcoldata2) == colnames(bpcts)))
bpcts$GeneID <- rownames(bpcts)

bpcts_geneanno <- bpcts |>
  merge(gene_anno, by = "GeneID", all = TRUE) |>
  merge(tax_anno,  by = "GeneID", all = TRUE)

# ── Build CAZyme grouped read table ───────────────────────────────────────────
# Expand dbcan_annotations (pipe-separated) into separate columns
aa <- bpcts_geneanno |>
  select(GeneID, everything()) |>
  mutate(dbcan_annotations = str_split(dbcan_annotations, "\\s*\\|\\s*")) |>
  unnest_wider(where(is.list), names_sep = "") |>
  as.data.frame()

aa <- aa[, c(2:54, 1)]

# Replace NA in first annotation column with "-" so all genes are counted
aa$dbcan_annotations1[is.na(aa$dbcan_annotations1)] <- "-"

# Identify annotation and sample columns
dbcan_cols  <- names(aa)[grepl("dbcan_annotations", names(aa))]
sample_cols <- colnames(aa)[1:20]

# Get unique CAZyme IDs
all_cazymes <- unique(as.vector(as.matrix(aa[, dbcan_cols])))
all_cazymes <- all_cazymes[!is.na(all_cazymes)]

# Pivot to long format: one row per gene × CAZyme annotation
aa_long <- aa |>
  pivot_longer(cols = all_of(dbcan_cols), values_to = "CAZyme", names_to = NULL) |>
  filter(!is.na(CAZyme), CAZyme %in% all_cazymes)

# ── Count matrix ──────────────────────────────────────────────────────────────
dbcantable_gen_raw <- aa_long |>
  pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "count") |>
  mutate(count = as.numeric(count)) |>
  group_by(CAZyme, sample) |>
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = sample, values_from = count, values_fill = 0) |>
  column_to_rownames("CAZyme")

# Ensure column order matches original
dbcantable_gen_raw <- dbcantable_gen_raw[, sample_cols]

# ── Taxonomy string per CAZyme ────────────────────────────────────────────────
taxaf <- aa_long |>
  filter(!is.na(taxa)) |>
  group_by(CAZyme) |>
  summarise(taxaf = paste(unique(taxa), collapse = "; "), .groups = "drop")

dbcantable_gen_raw$taxaf <- taxaf$taxaf[match(rownames(dbcantable_gen_raw), taxaf$CAZyme)]

# ── Save ──────────────────────────────────────────────────────────────────────
dbcantable_gen_raw$CAZID <- rownames(dbcantable_gen_raw)
write_tsv(
  dbcantable_gen_raw,
  file.path(BASE, paste0("MT_grouped_reads/dbcantable_raw_CORRECT_taxa_", DATE, ".tsv"))
)

dbcantable_gen_rawr <- read_tsv(file.path(BASE, paste0("MT_grouped_reads/dbcantable_raw_CORRECT_taxa_", DATE, ".tsv"))) %>%
  column_to_rownames("CAZID")
dbcantable <- dbcantable_gen_rawr[rownames(dbcantable_gen_rawr) != "-", ]

# ── Drop outlier samples ───────────────────────────────────────────────────────
drops      <- c("I6_CA14_1_MT_PLANC", "I8_CA58_1_MT_PLANC", "taxaf")
dbcantable <- dbcantable[, !names(dbcantable) %in% drops]

bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(dbcantable), ]
bpcoldata3$group <- gsub(" + ", "_", bpcoldata3$group, fixed = TRUE)
bpcoldata3$name  <- rownames(bpcoldata3)

# ── Taxa contribution to CAZymes, per condition ───────────────────────────────
conditions <- c("CA_S3_AS", "CA_AS", "PHA_G1_AS", "PHA_AS")

# Pre-build the long CAZyme x taxa x sample table once (all samples, all conditions)
cazyme_taxa_long_all <- aa_long |>
  filter(!is.na(taxa), CAZyme != "-") |>   # <-- drop unannotated genes
  pivot_longer(
    cols      = all_of(sample_cols),
    names_to  = "sample",
    values_to = "count"
  ) |>
  mutate(count = as.numeric(count))

# Map sample -> condition/group once
sample_group_map <- bpcoldata3 |>
  select(name, group)

cazyme_taxa_long_all <- cazyme_taxa_long_all |>
  left_join(sample_group_map, by = c("sample" = "name"))

# ── Loop over conditions, summarise contribution within each ─────────────────
cazyme_taxa_contrib_list <- map(conditions, function(cond) {
  cazyme_taxa_long_all |>
    filter(group == cond) |>
    group_by(CAZyme, taxa) |>
    summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") |>
    filter(total_count > 0) |>
    group_by(CAZyme) |>
    mutate(pct_of_cazyme = total_count / sum(total_count) * 100) |>
    ungroup() |>
    mutate(condition = cond)
})

names(cazyme_taxa_contrib_list) <- conditions

cazyme_taxa_contrib_all <- bind_rows(cazyme_taxa_contrib_list) |>
  arrange(condition, CAZyme, desc(total_count))

top_taxa_per_cazyme_all <- cazyme_taxa_contrib_all |>
  group_by(condition, CAZyme) |>
  slice_max(order_by = total_count, n = 1) |>
  ungroup()

write_tsv(
  cazyme_taxa_contrib_all,
  file.path(BASE, paste0("MT_grouped_reads/CAZyme_taxa_contrib_allconditions_", DATE, ".tsv"))
)


# ── MaAsLin2 helper ───────────────────────────────────────────────────────────
# Transpose: rows = samples, columns = features
dbcantable_t <- as.data.frame(t(dbcantable))
stopifnot(all(rownames(dbcantable_t) == rownames(bpcoldata3)))

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

# ── Comparison 1: CA_S3_AS vs CA_AS ───────────────────────────────────────────
meta1 <- bpcoldata3 |> filter(group %in% c("CA_S3_AS", "CA_AS"))
fit1  <- run_maaslin(
  data         = dbcantable_t[rownames(meta1), ],
  meta         = meta1,
  output       = file.path(OUT_DIR, paste0("CAZ_module_maaslin2_CA_", DATE)),
  fixed_effect = "group",
  reference    = "group,CA_AS"
)

# ── Comparison 2: PHA_G1_AS vs PHA_AS ────────────────────────────────────────
meta2 <- bpcoldata3 |> filter(group %in% c("PHA_G1_AS", "PHA_AS"))
fit2  <- run_maaslin(
  data         = dbcantable_t[rownames(meta2), ],
  meta         = meta2,
  output       = file.path(OUT_DIR, paste0("CAZ_module_maaslin2_PHA_", DATE)),
  fixed_effect = "group",
  reference    = "group,PHA_AS"
)

# ── Build combined results matrix ─────────────────────────────────────────────
res_CA  <- fit1$results
res_PHA <- fit2$results

all_cogs <- union(res_CA$feature, res_PHA$feature)

combined_mat <- data.frame(
  lfc_CA  = setNames(res_CA$coef,  res_CA$feature)[all_cogs],
  q_CA    = setNames(res_CA$qval,  res_CA$feature)[all_cogs],
  lfc_PHA = setNames(res_PHA$coef, res_PHA$feature)[all_cogs],
  q_PHA   = setNames(res_PHA$qval, res_PHA$feature)[all_cogs],
  row.names = all_cogs
) |>
  mutate(
    lfc_CA  = replace_na(lfc_CA,  0),
    lfc_PHA = replace_na(lfc_PHA, 0),
    q_CA    = replace_na(q_CA,    1),
    q_PHA   = replace_na(q_PHA,   1),
    CAZyme  = all_cogs
  )

combined_mat <- combined_mat[,c(5,1:4)]
out_tsv <- file.path(OUT_DIR, paste0("CAZ_maaslin2_CAPHA_coeffpadj_", DATE, ".tsv"))
write_tsv(combined_mat, out_tsv)

# ── Volcano plot helper ───────────────────────────────────────────────────────
make_volcano <- function(df, lfc_col, q_col, title, n_up = 10, n_down = 10) {
  
  df$diffexpressed <- case_when(
    df[[lfc_col]] > 0 & df[[q_col]] < 0.05 ~ "UP",
    df[[lfc_col]] < 0 & df[[q_col]] < 0.05 ~ "DOWN",
    TRUE ~ "NO"
  )
  
  top_up   <- df |> filter(diffexpressed == "UP")   |> slice_min(order_by = .data[[q_col]], n = n_up)
  top_down <- df |> filter(diffexpressed == "DOWN")  |> slice_min(order_by = .data[[q_col]], n = n_down)
  df$delabel <- ifelse(df$CAZyme %in% c(top_up$CAZyme, top_down$CAZyme), df$CAZyme, NA)
  
  ggplot(df, aes(x = .data[[lfc_col]], y = -log10(.data[[q_col]]),
                 col = diffexpressed, label = delabel)) +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
    geom_point(size = 1) +
    scale_color_manual(
      values = c("blue", "grey", "red"),
      labels = c("Downregulated", "Not significant", "Upregulated")
    ) +
    coord_cartesian(ylim = c(0, 4), xlim = c(-6, 6)) +
    scale_x_continuous(breaks = seq(-6, 6, 2)) +
    labs(title = title, x = "Coefficient", y = expression("-log"[10]*"q-value")) +
    geom_label_repel(size = 3, show.legend = FALSE,
                     min.segment.length = unit(0, "lines"), max.overlaps = Inf) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text       = element_text(size = 12),
      axis.title      = element_text(size = 12),
      legend.title    = element_blank(),
      legend.text     = element_text(size = 12),
      legend.position = "bottom",
      plot.margin     = margin(5, 5, 5, 5)
    )
}

# ── Plot ──────────────────────────────────────────────────────────────────────
# Combine CAZ_sub volcano plots and combine that with Chord plots in Powerpoint for final figure
df <- read_tsv(out_tsv) |> as.data.frame()

np1 <- make_volcano(df, "lfc_CA",  "q_CA",  "CA + S3 + AS vs. CA + AS")
np2 <- make_volcano(df, "lfc_PHA", "q_PHA", "PHA + G1 + AS vs. PHA + AS", n_up = 20)

np3 <- ggarrange(np1, np2, ncol = 2, nrow = 1,
                 common.legend = TRUE, legend = "bottom",
                 labels = c("A)", "B)"))
np3

ggsave(
  plot     = np3,
  filename = file.path(FIG_DIR, paste0("CAZ_volcfig_", DATE, ".pdf")),
  width = 12, height = 6, units = "in", dpi = 300, bg = "white"
)