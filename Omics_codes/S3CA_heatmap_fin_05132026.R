# =============================================================================
# S3CA DESeq2 — Volcano Plot and CAZyme Heatmap
# Colleen Ahern | Last updated: 05/13/2026
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(DRnaSeq)
library(ape)
library(ggrepel)
library(pheatmap)
library(grid)
library(Biostrings)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    <- "/Users/colleenahern/Documents/Magda_BPs_experiment"
RNA_DIR <- file.path(BASE, "RNAseq_S3CA")
FIG_DIR <- file.path(BASE, "Figures")
DATE    <- "05132026"

# =============================================================================
# 1.  Load count data and apply TPM filter
# =============================================================================

s3seqdata  <- read.delim(file.path(RNA_DIR, "featcount/s3ca_countmat_07302025.txt"),
                         stringsAsFactors = FALSE)
s3sampleinfo <- read.delim(file.path(RNA_DIR, "DESeq2_analysis/s3sampleinfo.txt"),
                           stringsAsFactors = TRUE)

# set rownames and drop non-count columns
s3countdata <- s3seqdata %>%
  column_to_rownames("Geneid") %>%
  dplyr::select(-(1:5))

stopifnot(all(colnames(s3countdata) %in% s3sampleinfo$libraryName))
colnames(s3countdata) <- s3sampleinfo$sampleName

# drop S3CAAS samples and S3_5 outlier
s3countdata  <- s3countdata[, c(1:4, 6:10)]
s3sampleinfo <- s3sampleinfo[c(1:4, 6:10), ]

# TPM normalization
gl  <- s3seqdata$Length
tpm <- as.data.frame(tpm(s3countdata, gl))
tpm$GeneID <- rownames(tpm)

# filter: keep genes where mean TPM >= 2 in at least one condition
tpm <- tpm %>%
  filter(
    rowMeans(select(., S3_1, S3_2, S3_3, S3_4))                      >= 2 |
      rowMeans(select(., S3CA_1, S3CA_2, S3CA_3, S3CA_4, S3CA_5))      >= 2
  )
# reduces from 27677 to ~13515 genes

# =============================================================================
# 2.  Load DESeq2 results and filter
# =============================================================================

vary_43 <- read_csv(file.path(RNA_DIR, "DESeq2_analysis/s3_DESeq_ALL_08222025.csv")) %>%
  select(GeneID = 8, S3CA_S3_log2FoldChange = 2, S3CA_S3_padj = 6, 7) %>%
  filter(!is.na(S3CA_S3_padj)) %>%                          # 25498 -> 22493
  inner_join(select(tpm, GeneID), by = "GeneID") %>%        # 22493 -> 13501
  arrange(S3CA_S3_padj)

# =============================================================================
# 3.  Add JGI gene annotations
# =============================================================================

gene_anno <- read.gff(file.path(BASE, "S3_JGI_files/Neolan1_GeneCatalog_20200610.gff3"),
                      GFF3 = TRUE) %>%
  filter(type == "gene") %>%
  dplyr::select(seqid, everything()) %>%
  mutate(attributes = str_split(attributes, ";")) %>%
  unnest_wider(where(is.list), names_sep = "") %>%
  select(9, 12, 13, 14) %>%
  rename(GeneID = 1, Protein_name = 2, ProteinID = 3, TranscriptID = 4) %>%
  mutate(
    GeneID       = gsub("ID=", "", GeneID),
    Protein_name = gsub("product_name=", "", Protein_name),
    ProteinID    = gsub("proteinId=", "", ProteinID),
    TranscriptID = gsub("transcriptId=", "", TranscriptID)
  )

vary_43 <- vary_43 %>%
  inner_join(gene_anno, by = "GeneID") %>%
  arrange(S3CA_S3_padj) %>%
  mutate(Name = paste0(Protein_name, " (proteinID = ", ProteinID, ")"))

write_csv(vary_43,
          file.path(RNA_DIR, "DESeq2_analysis/Tables/s3_vary43_s3cavs3_anno_05132026.csv"))

# =============================================================================
# 4.  Add CAZyme annotations from dbCAN
# =============================================================================

cazlist <- read.delim(file.path(RNA_DIR, "dbcan/overview.txt")) %>%
  filter(X.ofTools > 1) %>%
  mutate(
    Gene.ID = gsub("jgi\\|Neolan1\\|", "", Gene.ID),
    Gene.ID = gsub("\\|.*", "", Gene.ID),
    DIAMOND = gsub(".*fasta\\+", "", DIAMOND)
  )

rownames(cazlist) <- NULL

# split HMMER into separate columns and deduplicate
cazlist <- cazlist %>%
  mutate(HMMER = str_split(HMMER, "\\+")) %>%
  unnest_wider(where(is.list), names_sep = "")

# remove confidence scores from HMMER columns
cazlist <- cazlist %>%
  mutate(across(3:9, ~ gsub("\\(.*\\)", "", .)))

# build deduplicated HMMER annotation
cazlist$HMMER_u <- apply(cazlist[, 3:9], 1, function(row) {
  vals <- unique(na.omit(as.character(row)))
  paste(vals, collapse = "+")
})

cazlist <- cazlist %>%
  select(1:9, HMMER_u, 10:12) %>%
  mutate(
    HMMER_u = str_split(HMMER_u, "\\+"),
    eCAMI   = str_split(eCAMI,   "\\+"),
    DIAMOND = str_split(DIAMOND, "\\+")
  ) %>%
  unnest_wider(HMMER_u, names_sep = "") %>%
  unnest_wider(eCAMI,   names_sep = "") %>%
  unnest_wider(DIAMOND, names_sep = "")

# summarise CAZyme family counts per gene
cazlist$CAZ <- apply(cazlist[, 10:22], 1, function(row) {
  z <- as.data.frame(table(as.character(row)))
  paste(paste0(z$Var1, " x", z$Freq), collapse = "; ")
})

# keep only families detected by >= 2 tools, clean up annotation string
caz_anno <- cazlist %>%
  mutate(
    CAZf     = str_split(CAZ, "; "),
  ) %>%
  unnest_wider(CAZf, names_sep = "") %>%
  mutate(across(starts_with("CAZf"), ~ if_else(grepl("x1", .), NA_character_, .))) %>%
  unite(CAZ_name, starts_with("CAZf"), na.rm = TRUE, sep = "; ") %>%
  mutate(
    CAZ_anno = gsub(" x[0-9]+", "", CAZ_name)
  ) %>%
  filter(CAZ_anno != "") %>%
  rename(ProteinID = Gene.ID) %>%
  select(ProteinID, CAZ_anno)

write_csv(caz_anno, file.path(RNA_DIR, "dbcan/cazlist_05132026.csv"))

vary_43caz <- vary_43 %>%
  left_join(caz_anno, by = "ProteinID")

write_csv(vary_43caz,
          file.path(RNA_DIR, "DESeq2_analysis/Tables/s3_vary43_s3cavss3_anno_cazanno_05132026.csv"))

# =============================================================================
# 4b.  Import Wilken and Lawson datasets to compare to my RNAseq data
# =============================================================================

vary_43caz <- read_csv(file.path(RNA_DIR, "DESeq2_analysis/Tables/s3_vary43_s3cavss3_anno_cazanno_05132026.csv"))
stelmocell <- read.delim(file.path(RNA_DIR, "dbcan/cazyreport.txt"))
gh_doc_genes <- stelmocell %>%
  filter(grepl("GH", description) & grepl("DOC", description)) %>%
  rename(ProteinID = proteinId) %>%
  rename_with(~ paste0("Wilken_", .), .cols = -ProteinID)

scaffold_names <- names(readAAStringSet(file.path(BASE, "S3_JGI_files/Neolan1_GeneCatalog_proteins_20200610.aa.fasta")))
claw_table <- tibble(header = scaffold_names) %>%
  separate(header, into = c("jgi", "genome", "ProteinID", "orf"), sep = "\\|", extra = "merge") %>%
  select(ProteinID, orf)

cldata <- read_csv(file.path(RNA_DIR, "DESeq2_analysis/media-3.csv")) %>%
  filter(genome == "Neolan1") %>%
  select(-...1)

cltabdat <- merge(claw_table, cldata, by = "orf", all.y = TRUE) %>%
  rename_with(~ paste0("Lawson_", .), .cols = -ProteinID)

vary_43sigpos <- vary_43caz %>%
  filter(S3CA_S3_log2FoldChange > 0, S3CA_S3_padj < 0.05)

table(vary_43sigpos$ProteinID %in% cltabdat$ProteinID) # 285 of my N. lanati genes that were significantly upregulated are expressed in the Lawson et al dataset
table(vary_43sigpos$ProteinID %in% gh_doc_genes$ProteinID) # 108 of my N. lanati genes that were significantly upregulated are were identified as cellulosome-related GHs by Wilken et al

vary_43caz_lit <- merge(vary_43caz, gh_doc_genes, by = "ProteinID", all.x = TRUE)
vary_43caz_lit <- merge(vary_43caz_lit, cltabdat, by = "ProteinID", all.x = TRUE)
write_csv(vary_43caz_lit,
          file.path(RNA_DIR, "DESeq2_analysis/Tables/vary_43caz_lit_07212026.csv"))

# =============================================================================
# 5.  Volcano plot
# =============================================================================

df <- vary_43caz %>%
  mutate(diffexpressed = case_when(
    S3CA_S3_log2FoldChange >  1 & S3CA_S3_padj < 0.05 ~ "UP",
    S3CA_S3_log2FoldChange < -1 & S3CA_S3_padj < 0.05 ~ "DOWN",
    TRUE ~ "NO"
  ))

# label top 10 up and top 10 down by padj
top_up   <- df %>% filter(diffexpressed == "UP")   %>% slice_min(S3CA_S3_padj, n = 10) %>% pull(GeneID)
top_down <- df %>% filter(diffexpressed == "DOWN")  %>% slice_min(S3CA_S3_padj, n = 10) %>% pull(GeneID)
df$delabel <- ifelse(df$GeneID %in% c(top_up, top_down), df$GeneID, NA)

volcano_theme <- theme_classic(base_size = 20) +
  theme(
    axis.title.y    = element_text(face = "bold", margin = margin(0, 20, 0, 0), size = rel(1.1), color = "black"),
    axis.title.x    = element_text(hjust = 0.5, face = "bold", margin = margin(20, 0, 0, 0), size = rel(1.1), color = "black"),
    plot.title      = element_text(hjust = 0.5),
    legend.title    = element_blank(),
    legend.position = "bottom"
  )

volcano_plot <- ggplot(df, aes(x = S3CA_S3_log2FoldChange, y = -log10(S3CA_S3_padj),
                               col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
  geom_point(size = 1) +
  scale_color_manual(
    values = c(DOWN = "blue", NO = "grey", UP = "red"),
    labels = c("Downregulated", "Not significant", "Upregulated")
  ) +
  coord_cartesian(ylim = c(0, 40), xlim = c(-4, 7)) +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"q-value")) +
  scale_x_continuous(breaks = seq(-4, 7, 2)) +
  geom_label_repel(size = 3, show.legend = FALSE, min.segment.length = unit(0, "lines"),
                   max.overlaps = Inf) +
  theme_classic() +
  theme(
    text         = element_text(size = 12),
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 12),
    legend.text  = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

volcano_plot

ggsave(plot     = volcano_plot,
       filename = file.path(FIG_DIR, paste0("S3CA_volcano_", DATE, ".pdf")),
       width = 8, height = 8, units = "in", dpi = 300, bg = "white")

ggsave(plot     = volcano_plot,
       filename = file.path(FIG_DIR, paste0("S3CA_volcano_", DATE, ".tiff")),
       width = 8, height = 8, units = "in", dpi = 200, bg = "white")

# =============================================================================
# 6.  CAZyme heatmap
# =============================================================================

# filter to high LFC, significant genes with CAZyme annotation or alphafold hits
vary_43filt <- vary_43caz %>%
  filter(abs(S3CA_S3_log2FoldChange) > 3, S3CA_S3_padj < 0.01)

alp <- read_csv(file.path(RNA_DIR, "DESeq2_analysis/Alphafold/Alphafold_genes_05082026.csv")) %>%
  mutate(across(where(is.character), ~ gsub(",(?! )", ", ", ., perl = TRUE)))

pids <- c(259404, 987317, 997245, 1036142, 1048004, 1056260, 1289410,
          1293984, 1298002, 1433957, 1653640, 1675787, 1725095, 1751806)

vary_44 <- vary_43filt %>%
  mutate(ProteinID = as.character(ProteinID)) %>%
  left_join(alp %>% mutate(ProteinID = as.character(ProteinID)), by = "ProteinID") %>%
  filter(!is.na(CAZ_anno) | ProteinID %in% as.character(pids)) %>%
  mutate(
    Namealpha = if_else(!is.na(Foldseekanno),
                        paste0(Protein_name, " (Foldseek: ", Foldseekanno, ")"),
                        Protein_name),
    Namefin   = if_else(!is.na(CAZ_anno),
                        paste0(Namealpha, "; CAZyme family: ", CAZ_anno),
                        Namealpha),
    Namefin   = gsub("\t", " ", Namefin),
    Namefin   = sub("^([a-z])", "\\U\\1", Namefin, perl = TRUE),
    Namefin   = paste0(ProteinID, " - ", Namefin)
  )

# heatmap
my_palette <- colorRampPalette(c("blue", "black", "red"))(299)
newnames   <- "CA + S3 vs. S3"

s3_heatmap <- pheatmap(
  as.data.frame(vary_44$S3CA_S3_log2FoldChange),
  col          = my_palette,
  cellwidth    = 50,
  cellheight   = 15,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  fontsize_row  = 10,
  fontsize  = 12,
  labels_col   = newnames,
  labels_row   = vary_44$Namefin,
  border_color = "white",
  breaks       = seq(-5, 7, length.out = 300),
  angle_col    = 0
)

# =============================================================================
# 7.  Monoculture biogas plots
# =============================================================================

library(ggpubr)
library(cowplot)

biogas_theme <- theme_classic() +
  theme(
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 12),
    plot.title   = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold", vjust = -2),
    plot.margin  = margin(23, 0, 5, 5)
  )

# ── S3CA monoculture ──────────────────────────────────────────────────────────
s3ca2 <- read_csv(file.path(BASE, "biogas_data/S3CA_mono_gas_org_05062026.csv"))
s3ca2 <- s3ca2[, colSums(is.na(s3ca2)) == 0]
s3ca2vals <- s3ca2[, grepl("Day|Cumulative|avg|std", names(s3ca2))]

colnames(s3ca2vals) <- c("Day",
                         "S3_media_6", "S3_media_7", "S3_media_8",
                         "S3_avg", "S3_std",
                         "S3CA_media_22", "S3CA_media_23", "S3CA_media_24",
                         "S3CA_avg", "S3CA_std",
                         "media_40", "media_41", "media_42",
                         "media_avg", "media_std")

cols_s3ca <- c("CA + S3" = "#DC267F", "S3" = "#785EF0", "Media" = "#648FFF")

v1b <- ggplot(s3ca2vals, aes(x = Day)) +
  geom_line(aes(y = S3_avg,    color = "S3"))     +
  geom_point(aes(y = S3_avg,   color = "S3"))     +
  geom_errorbar(aes(ymax = S3_avg + S3_std, ymin = S3_avg - S3_std, color = "S3"), width = 0.3) +
  geom_line(aes(y = S3CA_avg,  color = "CA + S3")) +
  geom_point(aes(y = S3CA_avg, color = "CA + S3")) +
  geom_errorbar(aes(ymax = S3CA_avg + S3CA_std, ymin = S3CA_avg - S3CA_std, color = "CA + S3"), width = 0.3) +
  geom_line(aes(y = media_avg,  color = "Media"))  +
  geom_point(aes(y = media_avg, color = "Media"))  +
  geom_errorbar(aes(ymax = media_avg + media_std, ymin = media_avg - media_std, color = "Media"), width = 0.3) +
  scale_color_manual(name = "Condition", values = cols_s3ca, breaks = c("CA + S3", "S3", "Media")) +
  scale_x_continuous(breaks = c(0, 7, 13), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.7, 10), breaks = seq(0, 10, by = 2)) +
  labs(x = "Day", y = "Cumulative Biogas Production (psi)") +
  biogas_theme

v1b

# ── G1PHA monoculture ─────────────────────────────────────────────────────────
g1pha2 <- read_csv(file.path(BASE, "biogas_data/G1PHA_mono_gas_org_05062026.csv"))
g1pha2 <- g1pha2[, colSums(is.na(g1pha2)) == 0]
g1pha2vals <- g1pha2[, grepl("Day|Cumulative|avg|std", names(g1pha2))]

colnames(g1pha2vals) <- c("Day",
                          "G1PHA_media_1", "G1PHA_media_2", "G1PHA_media_5",
                          "G1PHA_avg", "G1PHA_std",
                          "G1_media_10", "G1_media_11", "G1_media_14",
                          "G1_avg", "G1_std",
                          "media_17", "media_18", "media_19",
                          "media_avg", "media_std")

cols_g1pha <- c("PHA + G1" = "#DC267F", "G1" = "#785EF0", "Media" = "#648FFF")

w1b <- ggplot(g1pha2vals, aes(x = Day)) +
  geom_line(aes(y = G1_avg,    color = "G1"))      +
  geom_point(aes(y = G1_avg,   color = "G1"))      +
  geom_errorbar(aes(ymax = G1_avg + G1_std, ymin = G1_avg - G1_std, color = "G1"), width = 0.2) +
  geom_line(aes(y = G1PHA_avg,  color = "PHA + G1")) +
  geom_point(aes(y = G1PHA_avg, color = "PHA + G1")) +
  geom_errorbar(aes(ymax = G1PHA_avg + G1PHA_std, ymin = G1PHA_avg - G1PHA_std, color = "PHA + G1"), width = 0.2) +
  geom_line(aes(y = media_avg,  color = "Media"))   +
  geom_point(aes(y = media_avg, color = "Media"))   +
  geom_errorbar(aes(ymax = media_avg + media_std, ymin = media_avg - media_std, color = "Media"), width = 0.2) +
  scale_color_manual(name = "Condition", values = cols_g1pha, breaks = c("PHA + G1", "G1", "Media")) +
  scale_x_continuous(breaks = c(0, 4, 8)) +
  scale_y_continuous(limits = c(-0.7, 10), breaks = seq(0, 10, by = 2)) +
  labs(x = "Day", y = "Cumulative Biogas Production (psi)") +
  biogas_theme

w1b

# =============================================================================
# 8.  Combined figure: biogas (A + B) and heatmap (C)
# =============================================================================
v1b <- v1b + theme(plot.margin = margin(23, 0, 5, 10))
w1b <- w1b + theme(plot.margin = margin(23, 10, 5, 5))

biogas_panel <- ggarrange(v1b, w1b, ncol = 2, nrow = 1, labels = c("A)", "B)"))

dp <- plot_grid(
  biogas_panel, s3_heatmap$gtable,
  nrow = 2, ncol = 1,
  rel_heights = c(1, 2),
  labels = c("", "C)")
)

dp

ggsave(plot     = dp,
       filename = file.path(FIG_DIR, paste0("Monocultures_S3CA_heatmap_fig_", DATE, ".pdf")),
       width = 13.5, height = 14, units = "in", dpi = 300, bg = "white")

ggsave(plot     = dp,
       filename = file.path(FIG_DIR, paste0("Monocultures_S3CA_heatmap_fig_", DATE, ".tiff")),
       width = 13.5, height = 14, units = "in", dpi = 200, bg = "white")

message("Done. Figures written to: ", FIG_DIR)