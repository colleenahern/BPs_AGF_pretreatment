# =============================================================================
# COG Category Differential Abundance Analysis — MaAsLin2
# Colleen Ahern | Last updated: 05/13/2026
#
# Two comparisons:
#   1. CA_S3_AS  vs. CA_AS
#   2. PHA_G1_AS vs. PHA_AS
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(Maaslin2)
library(pheatmap)
library(ggplotify)
library(cowplot)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    <- "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics"
MT_DIR  <- file.path(BASE, "MT_analysis/MT_analysis")
OUT_DIR <- file.path(BASE, "MT_grouped_reads/maaslin2")
FIG_DIR <- "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures"
META    <- file.path(BASE, "bp_metadata_01132025.txt")
DATE    <- "05132026"

# =============================================================================
# 1.  Load and prepare data
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
# 2.  Aggregate counts by COG category
# =============================================================================

aa <- bpcts_geneanno %>%
  select(GeneID, everything()) %>%
  mutate(COG_category = str_split(COG_category, "")) %>%
  unnest_wider(where(is.list), names_sep = "") %>%
  as.data.frame()

aa <- aa[, c(2:ncol(aa), 1)]
aa$COG_category1[is.na(aa$COG_category1)] <- "-"
sample_cols <- names(aa)[grepl("MT_PLANC", names(aa))]

# fast vectorized aggregation
aa[, sample_cols] <- lapply(aa[, sample_cols], function(x) suppressWarnings(as.numeric(x)))

cogtable_gen <- aa %>%
  pivot_longer(cols = starts_with("COG_category"),
               values_to = "COG_category",
               names_to  = NULL) %>%
  filter(!is.na(COG_category), COG_category != "-") %>%
  group_by(COG_category) %>%
  summarise(across(all_of(sample_cols), ~ sum(., na.rm = TRUE)))

# reshape aa to long format with one row per gene-COG combination
aa_long_cog <- aa %>%
  mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.)))) %>%
  pivot_longer(cols = starts_with("COG_category"),
               values_to = "COG_category",
               names_to  = NULL) %>%
  filter(!is.na(COG_category), COG_category != "-") %>%
  mutate(total_count = rowSums(across(all_of(sample_cols)), na.rm = TRUE)) %>%
  filter(total_count > 0)

# aggregate unique taxa per COG category
taxa_per_cog <- aa_long_cog %>%
  group_by(COG_category) %>%
  summarise(taxa = paste(unique(na.omit(taxa)), collapse = "; "))

# merge into cogtable_gen
cogtable_gen <- left_join(cogtable_gen, taxa_per_cog, by = "COG_category")

write_tsv(cogtable_gen,
          file.path(BASE, "MT_grouped_reads", paste0("COGtable_gen_raw_", DATE, ".tsv")))

# =============================================================================
# 3.  Prepare for MaAsLin2
# =============================================================================

cogtable_gen <- cogtable_gen %>%
  column_to_rownames("COG_category")

# drop outlier samples and unwanted COG categories
drop_samples <- c("I6_CA14_1_MT_PLANC", "I8_CA58_1_MT_PLANC", "taxa")

cogtable_gen_filt <- cogtable_gen[!names(cogtable_gen) %in% drop_samples]

bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(cogtable_gen_filt), ] %>%
  mutate(
    condition_cba = gsub("\\s+\\+\\s+", "_", condition_cba),
    condition_cba = gsub(" ", "_", condition_cba),
    group         = gsub("\\s+\\+\\s+", "_", group),
    name          = rownames(.)
  )

stopifnot(all(names(cogtable_gen_filt) == rownames(bpcoldata3)))

# rows = samples, cols = COG categories
cogtable_t <- as.data.frame(t(cogtable_gen_filt))

# =============================================================================
# 4.  MaAsLin2
# =============================================================================

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

# ── 1. CA_S3_AS vs. CA_AS ─────────────────────────────────────────────────────
meta1 <- bpcoldata3 %>% filter(group %in% c("CA_S3_AS", "CA_AS"))
fit1  <- run_maaslin(
  data         = cogtable_t[rownames(meta1), ],
  meta         = meta1,
  output       = file.path(OUT_DIR, paste0("maaslin2_COG_CA_", DATE)),
  fixed_effect = "group",
  reference    = "group,CA_AS"
)

# ── 2. PHA_G1_AS vs. PHA_AS ───────────────────────────────────────────────────
meta2 <- bpcoldata3 %>% filter(group %in% c("PHA_G1_AS", "PHA_AS"))
fit2  <- run_maaslin(
  data         = cogtable_t[rownames(meta2), ],
  meta         = meta2,
  output       = file.path(OUT_DIR, paste0("maaslin2_COG_PHA_", DATE)),
  fixed_effect = "group",
  reference    = "group,PHA_AS"
)

# =============================================================================
# 5.  Build heatmap matrices
# =============================================================================

res_CA  <- fit1$results
res_PHA <- fit2$results

res_CA <- read.table(
  file.path(OUT_DIR, paste0("maaslin2_COG_CA_", DATE), "all_results.tsv"),
  header = TRUE, sep = "\t"
)
res_PHA <- read.table(
  file.path(OUT_DIR, paste0("maaslin2_COG_PHA_", DATE), "all_results.tsv"),
  header = TRUE, sep = "\t"
)

all_cogs <- union(res_CA$feature, res_PHA$feature)

lfc_CA  <- setNames(res_CA$coef, res_CA$feature)
lfc_PHA <- setNames(res_PHA$coef, res_PHA$feature)
q_CA    <- setNames(res_CA$qval, res_CA$feature)
q_PHA   <- setNames(res_PHA$qval, res_PHA$feature)

cogcombo <- data.frame(
  lfc_CA  = lfc_CA[all_cogs],
  q_CA    = q_CA[all_cogs],
  lfc_PHA = lfc_PHA[all_cogs],
  q_PHA   = q_PHA[all_cogs],
  row.names = all_cogs
)

plot_mat <- data.frame(
  CA  = replace_na(lfc_CA[all_cogs],  0),
  PHA = replace_na(lfc_PHA[all_cogs], 0),
  row.names = all_cogs
)

padj_mat <- data.frame(
  CA  = replace_na(q_CA[all_cogs],  1),
  PHA = replace_na(q_PHA[all_cogs], 1),
  row.names = all_cogs
)

plot_mat_masked <- plot_mat
plot_mat_masked[padj_mat > 0.05] <- 0

# =============================================================================
# 6.  Annotations
# =============================================================================

cog_definitions <- c(
  J = "J - Translation, ribosomal structure and biogenesis",
  A = "A - RNA processing and modification",
  K = "K - Transcription",
  L = "L - Replication, recombination and repair",
  B = "B - Chromatin structure and dynamics",
  D = "D - Cell cycle control, cell division, chromosome partitioning",
  V = "V - Defense mechanisms",
  T = "T - Signal transduction mechanisms",
  M = "M - Cell wall/membrane/envelope biogenesis",
  N = "N - Cell motility",
  Z = "Z - Cytoskeleton",
  W = "W - Extracellular structures",
  U = "U - Intracellular trafficking, secretion, and vesicular transport",
  O = "O - Posttranslational modification, protein turnover, chaperones",
  X = "X - Mobilome: prophages, transposons",
  C = "C - Energy production and conversion",
  G = "G - Carbohydrate transport and metabolism",
  E = "E - Amino acid transport and metabolism",
  F = "F - Nucleotide transport and metabolism",
  H = "H - Coenzyme transport and metabolism",
  I = "I - Lipid transport and metabolism",
  P = "P - Inorganic ion transport and metabolism",
  Q = "Q - Secondary metabolites biosynthesis, transport and catabolism",
  R = "R - General function prediction only",
  S = "S - Function unknown",
  Y = "Y - Nuclear structure"
)

cog_domain_map <- c(
  J = "Information storage & processing",
  A = "Information storage & processing",
  K = "Information storage & processing",
  L = "Information storage & processing",
  B = "Information storage & processing",
  D = "Cellular processes & signaling",
  V = "Cellular processes & signaling",
  T = "Cellular processes & signaling",
  M = "Cellular processes & signaling",
  N = "Cellular processes & signaling",
  Z = "Cellular processes & signaling",
  W = "Cellular processes & signaling",
  U = "Cellular processes & signaling",
  O = "Cellular processes & signaling",
  X = "Cellular processes & signaling",
  C = "Metabolism",
  G = "Metabolism",
  E = "Metabolism",
  F = "Metabolism",
  H = "Metabolism",
  I = "Metabolism",
  P = "Metabolism",
  Q = "Metabolism",
  R = "Poorly characterized",
  S = "Poorly characterized",
  Y = "Nuclear structure"
)

cogcombo$cog_definition <- cog_definitions[rownames(cogcombo)]
cogcombo$cog_domain     <- cog_domain_map[rownames(cogcombo)]
cogcombo <- cogcombo %>%
  rownames_to_column("COG")
write_tsv(cogcombo,
          file.path(BASE, "MT_grouped_reads", paste0("COG_Maaslin2res_CAS3ASPHAG1AS_", DATE, ".tsv")))

# rename rows using full definitions
rownames(plot_mat_masked) <- cog_definitions[rownames(plot_mat_masked)]

cog_letters <- rownames(plot_mat_masked)  # already single letters before renaming
# extract letter from the renamed rows
cog_letters <- substr(rownames(plot_mat_masked), 1, 1)

row_annot <- data.frame(
  "COG Category" = cog_domain_map[cog_letters],
  row.names      = rownames(plot_mat_masked),
  check.names    = FALSE
)

domain_colors <- list(
  "COG Category" = c(
    "Information storage & processing" = "#0072B2",  # blue
    "Cellular processes & signaling"   = "#CC79A7",  # yellow-orange
    "Metabolism"                        = "#009E73",  # teal-green
    "Poorly characterized"              = "#E69F00",  # pink-purple
    "Nuclear structure"                 = "#D55E00"   # vermillion
  )
)

colnames(plot_mat_masked) <- c("CA + S3 + AS vs.\nCA + AS",
                               "PHA + G1 + AS vs.\nPHA + AS")

# =============================================================================
# 7.  Heatmap
# =============================================================================

cogph <- pheatmap(
  plot_mat_masked,
  cluster_cols         = FALSE,
  cluster_rows         = TRUE,
  color                = colorRampPalette(c("blue", "black", "red"))(100),
  breaks               = seq(-1, 1, length.out = 101),
  border_color         = "grey",
  annotation_row       = row_annot,
  annotation_colors    = domain_colors,
  annotation_names_row = FALSE,
  fontsize             = 12,
  fontsize_row         = 10,
  angle_col            = 0,
  cellwidth            = 120,
  cellheight           = 13
)

message("Done. Outputs written to: ", OUT_DIR)
