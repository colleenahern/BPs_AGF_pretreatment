# C. B. Ahern
# 04/20/2026
# Making cord plot for BPs MT data

library(readxl)
library(readr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(circlize)
library(grid)
library(cowplot)

taxa_gtdb <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/tax_gtdb_metassembly.tsv")
contig_gene <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/input/bowtie2_metassembly.stranded2.counts.txt", header = FALSE)
contig_gene <- contig_gene[-c(1),]
names(contig_gene) <- contig_gene[1,]
contig_gene <- contig_gene[-c(1),]
names(contig_gene) <- gsub("mapping/bowtie.meta_t.","",names(contig_gene))
names(contig_gene) <- gsub(".sorted.bam","",names(contig_gene))
contig_gene <- contig_gene %>%
  rename("contig" = "Chr")

taxam <- merge(contig_gene, taxa_gtdb, by = "contig", all = TRUE)
taxam <- taxam[,-c(3:6)]
taxam <- taxam %>%
  rename("GeneID" = "Geneid")

bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
bpcoldata <- bpcoldata %>%
  column_to_rownames("name")
bpcoldata2 <- bpcoldata[rownames(bpcoldata) %in% colnames(taxam[3:22]),]

# Input gene annotation info
gene_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")
gene_annof <- gene_anno[,c(1,7:23)]
  
taxamg <- merge(taxam, gene_annof, by = "GeneID", all = TRUE)

# CA + S3 + AS
# taxam_cas3as <- taxamg[,c(1,2,21,22,3,11,16,23:46)]
cas3asrn <- rownames(bpcoldata2[bpcoldata2$group == "CA + S3 + AS",])
taxam_cas3as <- taxamg[, c(colnames(taxamg)[1:2], cas3asrn, colnames(taxamg)[23:46])]

# Check to make sure of right groups
bpcoldata2[rownames(bpcoldata2) %in% cas3asrn, 7]
bpcoldata2[rownames(bpcoldata2) %in% cas3asrn, 8]

taxam_cas3as_long <- taxam_cas3as %>%
  pivot_longer(
    cols = 3:7,   # your sample columns
    names_to = "Sample",
    values_to = "Counts"
  )

taxam_cas3as_long3b <- taxam_cas3as_long %>%
  mutate(
    phylum = ifelse(is.na(phylum) | phylum == "", "Others", phylum),
    CAZ_class = str_extract_all(dbcan_annotations, "\\b(GH|PL|CE|AA|GT|CBM)\\d+")) %>%
  filter(!is.na(CAZ_class), !is.na(phylum))

taxam_cas3as_long3b$Counts <- as.numeric(taxam_cas3as_long3b$Counts) # convert from character to numeric
taxam_cas3as_long3b$Counts 

matb <- taxam_cas3as_long3b %>%
  group_by(Sample) %>%
  mutate(rel = Counts / sum(Counts, na.rm = TRUE)) 

write_tsv(matb, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/CAS3AS_chordplotdata_04282026.tsv")

# check to make sure I am normalizing correctly - divide each count in a sample by the total CAZyme counts in that sample
I8_CA58_1_MT_PLANC_cazcount <- sum(as.numeric(taxam_cas3as$I8_CA58_1_MT_PLANC[!is.na(taxam_cas3as$dbcan_annotations)]))
I8_CA58_1_MT_PLANC_cazcount
I9_CA61_1_MT_PLANC_cazcount <- sum(as.numeric(taxam_cas3as$I9_CA61_1_MT_PLANC[!is.na(taxam_cas3as$dbcan_annotations)]))
I9_CA61_1_MT_PLANC_cazcount
I10_CA62_1_MT_PLANC_cazcount <- sum(as.numeric(taxam_cas3as$I10_CA62_1_MT_PLANC[!is.na(taxam_cas3as$dbcan_annotations)]))
I10_CA62_1_MT_PLANC_cazcount
I1_CA61_2_MT_PLANC_cazcount <- sum(as.numeric(taxam_cas3as$I1_CA61_2_MT_PLANC[!is.na(taxam_cas3as$dbcan_annotations)]))
I1_CA61_2_MT_PLANC_cazcount
`I15-8_MT_PLANC_BSE_S15_cazcount` <- sum(as.numeric(taxam_cas3as$`I15-8_MT_PLANC_BSE_S15`[!is.na(taxam_cas3as$dbcan_annotations)]))
`I15-8_MT_PLANC_BSE_S15_cazcount`
  
mat2b <- matb %>%
  unnest(CAZ_class) %>%
  mutate(CAZ_class = str_extract(CAZ_class, "GH|PL|CE|AA|GT|CBM")) %>%
  group_by(Sample, phylum, CAZ_class) %>%
  summarise(rep_value = sum(rel), .groups = "drop") %>%
  group_by(phylum, CAZ_class) %>%
  summarise(value = mean(rep_value), .groups = "drop")

top_phylab <- mat2b %>%
  group_by(phylum) %>%
  summarise(counts = sum(value, na.rm = TRUE), .groups = "drop") %>%
  slice_max(counts, n = 10) %>%
  pull(phylum)

mat3b <- mat2b %>%
  mutate(phylum = ifelse(phylum %in% top_phylab, phylum, "Others")) %>%
  mutate(phylum = str_remove(phylum, "^p_"))

# draw cord plot
# install.packages("circlize")
library(circlize)

phylab <- unique(mat3b$phylum)
phylab
cazymesb <- unique(mat3b$CAZ_class)
cazymesb

phylum_colorsb <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(phylab)), phylab)
cazyme_colorsb <- setNames(colorRampPalette(brewer.pal(8, "Set1"))(length(cazymesb)), cazymesb)
grid.colb <- c(phylum_colorsb, cazyme_colorsb)

# png("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/CAS3AS_chord_plot_04212026.png", width = 1900, height = 1800, res = 300)
pdf("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/CAS3AS_chord_plot_04212026.pdf", width = 6, height = 6, bg = "white")

circlize::circos.clear()
circos.par(canvas.xlim = c(-1.2, 1.2),
           canvas.ylim = c(-1.6, 1.1))

chordDiagram(mat3b,
             grid.col = grid.colb,
             transparency = 0.3,
             directional = 1,
             annotationTrack = "grid",
             preAllocateTracks = 1)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1],
                CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 1)
  },
  bg.border = NA
)

text(0, 1, "CA + S3 + AS", cex = 1.2, font = 2, adj = c(0.5, 0.5))

dev.off()


# CA +  AS
caasrn <- rownames(bpcoldata2[bpcoldata2$group == "CA + AS",])
taxam_caas <- taxamg[, c(colnames(taxamg)[1:2], caasrn, colnames(taxamg)[23:46])]

# Check to make sure of right groups
bpcoldata2[rownames(bpcoldata2) %in% caasrn, 7]
bpcoldata2[rownames(bpcoldata2) %in% caasrn, 8]

taxam_caas_long <- taxam_caas %>%
  pivot_longer(
    cols = 3:7,   # your sample columns
    names_to = "Sample",
    values_to = "Counts"
  )

taxam_caas_long3b <- taxam_caas_long %>%
  mutate(
    phylum = ifelse(is.na(phylum) | phylum == "", "Others", phylum),
    CAZ_class = str_extract_all(dbcan_annotations, "\\b(GH|PL|CE|AA|GT|CBM)\\d+")) %>%
  filter(!is.na(CAZ_class), !is.na(phylum))

taxam_caas_long3b$Counts <- as.numeric(taxam_caas_long3b$Counts) # convert from character to numeric
taxam_caas_long3b$Counts 

# average counts across replicates
mat_caas_b <- taxam_caas_long3b %>%
  group_by(Sample) %>%
  mutate(rel = Counts / sum(Counts, na.rm = TRUE)) 

write_tsv(mat_caas_b, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/CAAS_chordplotdata_04282026.tsv")

# check to make sure I am normalizing correctly - divide each count in a sample by the total CAZyme counts in that sample
cazcount3 <- sum(as.numeric(taxam_caas[!is.na(taxam_caas$dbcan_annotations),3]))
cazcount3
cazcount4 <- sum(as.numeric(taxam_caas[!is.na(taxam_caas$dbcan_annotations),4]))
cazcount4
cazcount5 <- sum(as.numeric(taxam_caas[!is.na(taxam_caas$dbcan_annotations),5]))
cazcount5
cazcount6 <- sum(as.numeric(taxam_caas[!is.na(taxam_caas$dbcan_annotations),6]))
cazcount6
cazcount7 <- sum(as.numeric(taxam_caas[!is.na(taxam_caas$dbcan_annotations),7]))
cazcount7

mat2_caas_b <- mat_caas_b %>%
  unnest(CAZ_class) %>%
  mutate(CAZ_class = str_extract(CAZ_class, "GH|PL|CE|AA|GT|CBM")) %>%
  group_by(Sample, phylum, CAZ_class) %>%
  summarise(rep_value = sum(rel), .groups = "drop") %>%
  group_by(phylum, CAZ_class) %>%
  summarise(value = mean(rep_value), .groups = "drop")

top_phylab_caas <- mat2_caas_b %>%
  group_by(phylum) %>%
  summarise(counts = sum(value, na.rm = TRUE), .groups = "drop") %>%
  slice_max(counts, n = 10) %>%
  pull(phylum)

mat3_caas_b <- mat2_caas_b %>%
  mutate(phylum = ifelse(phylum %in% top_phylab_caas, phylum, "Others")) %>%
  mutate(phylum = str_remove(phylum, "^p_"))

# draw cord plot
library(circlize)

phylab_caas <- unique(mat3_caas_b$phylum)
phylab_caas
cazymesb_caas <- unique(mat3_caas_b$CAZ_class)
cazymesb_caas

phylum_colorsb_caas <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(phylab_caas)), phylab_caas)
cazyme_colorsb_caas <- setNames(colorRampPalette(brewer.pal(8, "Set1"))(length(cazymesb_caas)), cazymesb_caas)
grid.colb_caas <- c(phylum_colorsb_caas, cazyme_colorsb_caas)

grid.colb_caas <- c(
  # Matching taxa from grid.colb — same colors
  "Others"            = "#1B9E77",
  "Bacillota A"       = "#AE6D1C",
  "Bacteroidota"      = "#9B58A5",
  "Chloroflexota"     = "#D8367D",
  "Desulfobacterota"  = "#749829",
  "Halobacteriota"    = "#C9930D",
  "Pseudomonadota"    = "#97722D",
  "Thermotogota"      = "#666666",
  "AA"                = "#E41A1C",
  "CBM"               = "#3F918B",
  "CE"                = "#896191",
  "GH"                = "#FF980A",
  "GT"                = "#C9992C",
  "PL"                = "#F781BF",
  # Unique to grid.colb_caas
  "Atribacterota"     = "#A16864",  # reusing Bacillota B color
  "Verrucomicrobiota" = "#BBA90B"   # reusing Fibrobacterota color
)

# png("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/caas_chord_plot_04222026.png", width = 1900, height = 1800, res = 300)
pdf("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/caas_chord_plot_04222026.pdf", width = 6, height = 6, bg = "white")
pdf("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/caas_chord_plot_05062026.pdf", width = 6, height = 6, bg = "white")


circlize::circos.clear()
# circos.par(canvas.xlim = c(-1.2, 1), canvas.ylim = c(-.9, 1.1))
# circos.par(canvas.xlim = c(0, 0), canvas.ylim = c(0, 0))
# circos.par(canvas.xlim = c(-1, 1), canvas.ylim = c(-1, 1))
circos.par(canvas.xlim = c(-1.2, 1.2),
           canvas.ylim = c(-1.6, 1.1),
           start.degree = 2)

chordDiagram(mat3_caas_b,
             grid.col = grid.colb_caas,
             transparency = 0.3,
             directional = 1,
             annotationTrack = "grid",
             preAllocateTracks = 1)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1],
                CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 1)
  },
  bg.border = NA
)


text(0, 1, "CA + AS", cex = 1.2, font = 2, adj = c(0.5, 0.5))
# Add panel label
# text(-1.5, 1.5, "D)", cex = 1.21, font = 2)
dev.off()

# PHA + G1 +  AS
phag1asrn <- rownames(bpcoldata2[bpcoldata2$group == "PHA + G1 + AS",])
taxam_phag1as <- taxamg[, c(colnames(taxamg)[1:2], phag1asrn, colnames(taxamg)[23:46])]

# Check to make sure of right groups
bpcoldata2[rownames(bpcoldata2) %in% phag1asrn, 7]
bpcoldata2[rownames(bpcoldata2) %in% phag1asrn, 8]

taxam_phag1as_long <- taxam_phag1as %>%
  pivot_longer(
    cols = 3:8,   # your sample columns
    names_to = "Sample",
    values_to = "Counts"
  )

taxam_phag1as_long3b <- taxam_phag1as_long %>%
  mutate(
    phylum = ifelse(is.na(phylum) | phylum == "", "Others", phylum),
    CAZ_class = str_extract_all(dbcan_annotations, "\\b(GH|PL|CE|AA|GT|CBM)\\d+")) %>%
  filter(!is.na(CAZ_class), !is.na(phylum))

taxam_phag1as_long3b$Counts <- as.numeric(taxam_phag1as_long3b$Counts) # convert from character to numeric
taxam_phag1as_long3b$Counts 

# average counts across replicates
mat_phag1as_b <- taxam_phag1as_long3b %>%
  group_by(Sample) %>%
  mutate(rel = Counts / sum(Counts, na.rm = TRUE)) 

write_tsv(mat_phag1as_b, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/PHAG1AS_chordplotdata_04282026.tsv")

# check to make sure I am normalizing correctly - divide each count in a sample by the total CAZyme counts in that sample
cazcount3 <- sum(as.numeric(taxam_phag1as[!is.na(taxam_phag1as$dbcan_annotations),3]))
cazcount3
cazcount4 <- sum(as.numeric(taxam_phag1as[!is.na(taxam_phag1as$dbcan_annotations),4]))
cazcount4
cazcount5 <- sum(as.numeric(taxam_phag1as[!is.na(taxam_phag1as$dbcan_annotations),5]))
cazcount5
cazcount6 <- sum(as.numeric(taxam_phag1as[!is.na(taxam_phag1as$dbcan_annotations),6]))
cazcount6
cazcount7 <- sum(as.numeric(taxam_phag1as[!is.na(taxam_phag1as$dbcan_annotations),7]))
cazcount7
cazcount8 <- sum(as.numeric(taxam_phag1as[!is.na(taxam_phag1as$dbcan_annotations),8]))
cazcount8

mat2_phag1as_b <- mat_phag1as_b %>%
  unnest(CAZ_class) %>%
  mutate(CAZ_class = str_extract(CAZ_class, "GH|PL|CE|AA|GT|CBM")) %>%
  group_by(Sample, phylum, CAZ_class) %>%
  summarise(rep_value = sum(rel), .groups = "drop") %>%
  group_by(phylum, CAZ_class) %>%
  summarise(value = mean(rep_value), .groups = "drop")

top_phylab_phag1as <- mat2_phag1as_b %>%
  group_by(phylum) %>%
  summarise(counts = sum(value, na.rm = TRUE), .groups = "drop") %>%
  slice_max(counts, n = 10) %>%
  pull(phylum)

mat3_phag1as_b <- mat2_phag1as_b %>%
  mutate(phylum = ifelse(phylum %in% top_phylab_phag1as, phylum, "Others")) %>%
  mutate(phylum = str_remove(phylum, "^p_"))

# draw cord plot
library(circlize)

phylab_phag1as <- unique(mat3_phag1as_b$phylum)
phylab_phag1as
cazymesb_phag1as <- unique(mat3_phag1as_b$CAZ_class)
cazymesb_phag1as

phylum_colorsb_phag1as <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(phylab_phag1as)), phylab_phag1as)
cazyme_colorsb_phag1as <- setNames(colorRampPalette(brewer.pal(8, "Set1"))(length(cazymesb_phag1as)), cazymesb_phag1as)
grid.colb_phag1as <- c(phylum_colorsb_phag1as, cazyme_colorsb_phag1as)

grid.colb_phag1as <- c(
  # Matching taxa from grid.colb — same colors
  "Others"           = "#1B9E77",
  "Bacillota A"      = "#AE6D1C",
  "Bacillota B"      = "#A16864",
  "Bacteroidota"     = "#9B58A5",
  "Chloroflexota"    = "#D8367D",
  "Desulfobacterota" = "#749829",
  "Halobacteriota"   = "#C9930D",
  "Pseudomonadota"   = "#97722D",
  "Thermotogota"     = "#666666",
  "AA"               = "#E41A1C",
  "CBM"              = "#3F918B",
  "CE"               = "#896191",
  "GH"               = "#FF980A",
  "GT"               = "#C9992C",
  "PL"               = "#F781BF",
  # Unique to grid.colb_phag1as
  "Atribacterota"    = "#BBA90B"  # reusing Fibrobacterota color
)

pdf("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/PHAG1AS_chord_plot_04232026.pdf", width = 6, height = 6, bg = "white")
pdf("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/PHAG1AS_chord_plot_05062026.pdf", width = 6, height = 6, bg = "white")


circlize::circos.clear()
circos.par(canvas.xlim = c(-1.2, 1.2),
           canvas.ylim = c(-1.6, 1.1)) #,
          # start.degree = 5)

chordDiagram(mat3_phag1as_b,
             grid.col = grid.colb_phag1as,
             transparency = 0.3,
             directional = 1,
             annotationTrack = "grid",
             preAllocateTracks = 1)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1],
                CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 1)
  },
  bg.border = NA
)


text(0, 1, "PHA + G1 + AS", cex = 1.2, font = 2, adj = c(0.5, 0.5))
# Add panel label
# text(-1.5, 1.5, "D)", cex = 1.21, font = 2)
dev.off()

# PHA + AS
phaasrn <- rownames(bpcoldata2[bpcoldata2$group == "PHA + AS",])
taxam_phaas <- taxamg[, c(colnames(taxamg)[1:2], phaasrn, colnames(taxamg)[23:46])]

# Check to make sure of right groups
bpcoldata2[rownames(bpcoldata2) %in% phaasrn, 7]
bpcoldata2[rownames(bpcoldata2) %in% phaasrn, 8]

taxam_phaas_long <- taxam_phaas %>%
  pivot_longer(
    cols = 3:6,   # your sample columns
    names_to = "Sample",
    values_to = "Counts"
  )

taxam_phaas_long3b <- taxam_phaas_long %>%
  mutate(
    phylum = ifelse(is.na(phylum) | phylum == "", "Others", phylum),
    CAZ_class = str_extract_all(dbcan_annotations, "\\b(GH|PL|CE|AA|GT|CBM)\\d+")) %>%
  filter(!is.na(CAZ_class), !is.na(phylum))

taxam_phaas_long3b$Counts <- as.numeric(taxam_phaas_long3b$Counts) # convert from character to numeric
taxam_phaas_long3b$Counts 

# average counts across replicates
mat_phaas_b <- taxam_phaas_long3b %>%
  group_by(Sample) %>%
  mutate(rel = Counts / sum(Counts, na.rm = TRUE)) 

write_tsv(mat_phaas_b, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/PHAAS_chordplotdata_04282026.tsv")

# check to make sure I am normalizing correctly - divide each count in a sample by the total CAZyme counts in that sample
cazcount3 <- sum(as.numeric(taxam_phaas[!is.na(taxam_phaas$dbcan_annotations),3]))
cazcount3
cazcount4 <- sum(as.numeric(taxam_phaas[!is.na(taxam_phaas$dbcan_annotations),4]))
cazcount4
cazcount5 <- sum(as.numeric(taxam_phaas[!is.na(taxam_phaas$dbcan_annotations),5]))
cazcount5
cazcount6 <- sum(as.numeric(taxam_phaas[!is.na(taxam_phaas$dbcan_annotations),6]))
cazcount6

mat2_phaas_b <- mat_phaas_b %>%
  unnest(CAZ_class) %>%
  mutate(CAZ_class = str_extract(CAZ_class, "GH|PL|CE|AA|GT|CBM")) %>%
  group_by(Sample, phylum, CAZ_class) %>%
  summarise(rep_value = sum(rel), .groups = "drop") %>%
  group_by(phylum, CAZ_class) %>%
  summarise(value = mean(rep_value), .groups = "drop")

top_phylab_phaas <- mat2_phaas_b %>%
  group_by(phylum) %>%
  summarise(counts = sum(value, na.rm = TRUE), .groups = "drop") %>%
  slice_max(counts, n = 10) %>%
  pull(phylum)

mat3_phaas_b <- mat2_phaas_b %>%
  mutate(phylum = ifelse(phylum %in% top_phylab_phaas, phylum, "Others")) %>%
  mutate(phylum = str_remove(phylum, "^p_"))

# draw cord plot
library(circlize)

phylab_phaas <- unique(mat3_phaas_b$phylum)
phylab_phaas
cazymesb_phaas <- unique(mat3_phaas_b$CAZ_class)
cazymesb_phaas

phylum_colorsb_phaas <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(phylab_phaas)), phylab_phaas)
cazyme_colorsb_phaas <- setNames(colorRampPalette(brewer.pal(8, "Set1"))(length(cazymesb_phaas)), cazymesb_phaas)
grid.colb_phaas <- c(phylum_colorsb_phaas, cazyme_colorsb_phaas)
grid.colb_phaas <- c(
  # Matching taxa from grid.colb — same colors
  "Others"            = "#1B9E77",
  "Bacillota A"       = "#AE6D1C",
  "Bacteroidota"      = "#9B58A5",
  "Chloroflexota"     = "#D8367D",
  "Desulfobacterota"  = "#749829",
  "Halobacteriota"    = "#C9930D",
  "Pseudomonadota"    = "#97722D",
  "Thermotogota"      = "#666666",
  "AA"                = "#E41A1C",
  "CBM"               = "#3F918B",
  "CE"                = "#896191",
  "GH"                = "#FF980A",
  "GT"                = "#C9992C",
  "PL"                = "#F781BF",
  # Unique to grid.colb_phaas
  "Atribacterota"     = "#A16864",  # reusing Bacillota B color
  "Verrucomicrobiota" = "#BBA90B"   # reusing Fibrobacterota color
)

pdf("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/PHAAS_chord_plot_04232026.pdf", width = 6, height = 6, bg = "white")
pdf("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/PHAAS_chord_plot_05062026.pdf", width = 6, height = 6, bg = "white")

circlize::circos.clear()
circos.par(canvas.xlim = c(-1.2, 1.2),
           canvas.ylim = c(-1.6, 1.1))

chordDiagram(mat3_phaas_b,
             grid.col = grid.colb_phaas,
             transparency = 0.3,
             directional = 1,
             annotationTrack = "grid",
             preAllocateTracks = 1)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1],
                CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 1)
  },
  bg.border = NA
)


text(0, 1, "PHA + AS", cex = 1.2, font = 2, adj = c(0.5, 0.5))
# Add panel label
# text(-1.5, 1.5, "D)", cex = 1.21, font = 2)
dev.off()


# 04/21/2026 combine CAZ_sub volcano plots with Chord plots for figure
library(png)
library(ggpubr)
library(magick)
library(pdftools)
library(gridGraphics)
library(grid)
library(ggpubr)

cas3as_chord <- image_read_pdf("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/CAS3AS_chord_plot_04212026.pdf") %>%
  rasterGrob(interpolate = TRUE) %>%
  as_ggplot() #+ theme(plot.margin = margin(20, 0, 5, 0))
cas3as_chord

cas3as_chord <- image_read("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/CAS3AS_chord_plot_04212026.png") %>%
  rasterGrob(interpolate = TRUE) %>%
  as_ggplot() #+ theme(plot.margin = margin(20, 0, 5, 0))
cas3as_chord

phag1as_chord <- image_read("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/phag1as_chord_plot_04222026.png") %>%
  rasterGrob(interpolate = TRUE) %>%
  as_ggplot() #+ theme(plot.margin = margin(20, 0, 5, 0))
phag1as_chord


volcfig <- ggarrange(np1, np2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom", labels = c("A)", "B)"))
volcfig
ggsave(plot = volcfig, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/CAZ_volcfig_04222026.pdf", width = 12, height = 6, units = "in", dpi = 300)


chordfig <- ggarrange(cas3as_chord, phag1as_chord, ncol=2, nrow=2, labels = c("C)", "D)"))
chordfig

volcchordfig <- ggarrange(volcfig, cas3as_chord, phag1as_chord, ncol=1, nrow=2)
volcchordfig

ggsave(plot = volcchordfig, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/CAZ_volcchord_fig_04222026.pdf", width = 10, height = 10, units = "in", dpi = 300)
ggsave(plot = volcchordfig, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/CAZ_volcchord_fig_04222026.tiff", width = 10, height = 10, units = "in", dpi = 200, bg="white")
