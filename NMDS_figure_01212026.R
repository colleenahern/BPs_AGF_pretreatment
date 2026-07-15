# 01/21/2024
# Colleen Ahern
# Final figure version: Using qiime2 output from Magda_BPs_experiment to create NMDS plot


# This is for ARCHAEA data (I rarefied 3,000 cutoff)

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

library(qiime2R)
library(phyloseq)
metadata_a <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast/metadata/metadata_a_12122024.txt", header=T, stringsAsFactors=F, sep="\t")

setwd("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast")

pswd <- qza_to_phyloseq(features="caqiime2b/arc2analysis/16S-table-nocmb-rarefied-3000_filtered.qza",
                        tree="caqiime2b/arc2analysis/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "output_arc2/qiime2/input/taxonomy.qza", 
                        metadata= "metadata/metadata_a_12122024.txt")

sample_variables(pswd)
names(pswd@sam_data) <- c("Target", "Experiment", "Test", "TestCond", "AF", "Timepoint", "Time",  "Replicate", "Group", "Condition")
pswd@sam_data

# For film
pswd_film <- subset_samples(pswd, Experiment == "BP_Film") 
pswd_film@sam_data
# pswd_film@sam_data$Name <- row.names(pswd_film@sam_data)
pswd_film@sam_data$Condition <- gsub(" film", "", pswd_film@sam_data$Condition)
pswd_film@sam_data$Condition
pswd_film@sam_data

# drop <- c("A114_62_3CA_S3", "A108_62_2CA_S3", "A117_62_FCA_S3")
# pswd_film = prune_samples(!(sample_names(pswd_film) %in% drop), pswd_film)

pswd.prop <- transform_sample_counts(pswd_film, function(otu) otu/sum(otu))

ord.nmds.bray_wd <- ordinate(pswd.prop, method="NMDS", distance="bray")

library(ggforce)
# deg_plot <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color="Substrate", shape = "Day", title="Bray NMDS") + scale_shape_binned(name = "Day", limits=c(0,90), breaks=c(0,45,90)) + geom_point(size = 4) 
# deg_plot <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color = "Group", title="Bray NMDS") + geom_point(size = 3) 

deg_plot <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color = "Condition", shape = "Time", title="Bray NMDS: Archaea") + 
  geom_point(size = 3) + 
  # ggforce::geom_mark_ellipse(aes(fill = Testcond,
  #                                color = Testcond)) + doesn't work - separates Condition by timepoint
  theme_classic() + 
  theme(text = element_text(size = 12), plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.title = element_text(size = 12), 
        legend.title = element_text(size = 12, face = "bold", vjust = -2), legend.text = element_text(size = 12), axis.text = element_text(size = 12))

# deg_plot <- plot_ordination(pswd.prop, ord.nmds.bray_wd, label = "Name") 
# deg_plot <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color = "Condition", title="Bray NMDS: Archaea") + 
#   geom_point(size = 4) + 
#   geom_point(size = 4) + 
#   stat_ellipse() +
#   theme_bw() + 
#   theme(text = element_text(size = 14)) 

deg_plot 

########################################################################################################################
### This is for BACTERIA data (I rarefied 10,000 cutoff)

setwd("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast")

pswd <- qza_to_phyloseq(features="caqiime2b/bac2analysis/16S-table-nocma-rarefied-10000_filtered.qza",
                        tree="caqiime2b/bac2analysis/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "output_bac2/qiime2/input/taxonomy.qza", 
                        metadata= "metadata/metadata_b_12122024.txt")

sample_variables(pswd)
names(pswd@sam_data) <- c("Target", "Experiment", "Test", "TestCond", "AF", "Timepoint", "Time",  "Replicate", "Group", "Condition")
sample_variables(pswd)

# For film
pswd_film <- subset_samples(pswd, Experiment == "BP_Film") 
pswd_film@sam_data$Condition <- gsub(" film", "", pswd_film@sam_data$Condition)
pswd_film@sam_data$Condition
print(pswd_film)

pswd.prop <- transform_sample_counts(pswd_film, function(otu) otu/sum(otu))

ord.nmds.bray_wd <- ordinate(pswd.prop, method="NMDS", distance="bray")

library(ggforce)
# deg_plot <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color="Substrate", shape = "Day", title="Bray NMDS") + scale_shape_binned(name = "Day", limits=c(0,90), breaks=c(0,45,90)) + geom_point(size = 4) 
# deg_plot <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color = "Group", title="Bray NMDS") + geom_point(size = 3) 

deg_plot2 <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color = "Condition", shape = "Time", title="Bray NMDS: Bacteria") + 
  geom_point(size = 3) + 
  theme_classic() + 
  theme(text = element_text(size = 12), plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.title = element_text(size = 12), 
        legend.title = element_text(size = 12, face = "bold", vjust = -2), legend.text = element_text(size = 12), axis.text = element_text(size = 12))

deg_plot2

# install.packages("ggpubr")
library(ggpubr)
deg_plot3 <- ggarrange(deg_plot, deg_plot2, ncol=2, nrow=1, common.legend = TRUE, legend="right", labels = c("A)", "B)"))
deg_plot3

# deg_plot3b <- deg_plot3 + theme(plot.margin = margin(0, 0, 50, 0)) 
# deg_plot3b
# library(gridExtra)
# grid.arrange(deg_plot3, ph, ncol=1, nrow=2)
# grid.arrange(heatmap_plot, scatter_plot, ncol = 2)
# 
# 
# library(patchwork)
# deg_plot3 + ph

library(cowplot)
nmdsancaligned_plot <- align_patches(deg_plot3, ancombcphgg)
# nmdsancplot <- plot_grid(nmdsancaligned_plot[[1]], nmdsancaligned_plot[[2]], nrow = 2, ncol = 1, rel_heights = c(1.5, 1))
nmdsancplot <- plot_grid(deg_plot3, ancombcphgg, nrow = 2, ncol = 1, labels = c("", "C)"))
nmdsancplot
# ggsave(plot = nmdsancplot, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/16S_figure_02052026.pdf", width = 10, height = 8, units = "in", dpi = 300)
# ggsave(plot = nmdsancplot, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/16S_figure_02052026.tiff", width = 10, height = 8, units = "in", dpi = 200, bg = "white")

# changed color gradient for heatmap so upregulated is shown better
ggsave(plot = nmdsancplot, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/16S_figure_02232026.pdf", width = 10, height = 8, units = "in", dpi = 300)
ggsave(plot = nmdsancplot, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/16S_figure_02232026.tiff", width = 10, height = 8, units = "in", dpi = 200, bg = "white")


# Permanova to check statistical significance
install.packages("vegan")  
library(vegan)

##### ARCHAEA: subset phyloseq object for CA + S3 + AS vs CA + AS/PHA + AS comparison only
psarc_sub <- subset_samples(pswd_film, Condition %in% c("CA + S3 + AS", "CA + AS", "PHA + AS"))
psarc_sub@sam_data
psarc_sub@sam_data$Subcond <- c("A", "A", "A", "A", "B", "B", "A", "A", "A", "A", "B")
psarc_sub@sam_data

# subset phyloseq object for PHA + G1 + AS vs CA + AS/PHA + AS comparison only
psarc_sub <- subset_samples(pswd_film, Condition %in% c("PHA + G1 + AS", "CA + AS", "PHA + AS"))
psarc_sub@sam_data
psarc_sub@sam_data$Subcond <- c("B", "B", "A", "A", "A", "A", "B", "A", "A", "A")
psarc_sub@sam_data

# OTU table
otu <- otu_table(psarc_sub)
otu_mat <- as(otu, "matrix")

# Ensure taxa are columns
if (taxa_are_rows(psarc_sub)) {
  otu_mat <- t(otu_mat)
}

# Metadata
meta <- data.frame(sample_data(psarc_sub))

dist_mat <- vegdist(otu_mat, method = "bray")
adonis2(dist_mat ~ Subcond, data = meta)

dist_mat <- vegdist(otu_mat, method = "bray")
adonis2(dist_mat ~ Subcond, data = meta)


##### BACTERIA: subset phyloseq object for CA + S3 + AS vs CA + AS/PHA + AS comparison only
pswd_film@sam_data
psbac_sub <- subset_samples(pswd_film, Condition %in% c("CA + S3 + AS", "CA + AS"))
psbac_sub@sam_data

# subset phyloseq object for PHA + G1 + AS vs CA + AS/PHA + AS comparison only
psbac_sub <- subset_samples(pswd_film, Condition %in% c("PHA + G1 + AS", "PHA + AS"))
psbac_sub@sam_data

# OTU table
otu <- otu_table(psbac_sub)
otu_mat <- as(otu, "matrix")

# Ensure taxa are columns
if (taxa_are_rows(psbac_sub)) {
  otu_mat <- t(otu_mat)
}

# Metadata
meta <- data.frame(sample_data(psbac_sub))

dist_mat <- vegdist(otu_mat, method = "bray")
adonis2(dist_mat ~ Condition, data = meta)

bd <- betadisper(dist_mat, meta$Condition)
anova(bd)




## Results:
# Archaea CA + S3 + AS vs CA + AS/PHA + AS Pr(>F) = 0.01; R2 = 0.28431; F = 3.5752; dispersion =  0.02523
# Archaea PHA + G1 + AS vs CA + AS/PHA + AS Pr(>F) = 0.016; R2 = 0.47321; F = 7.1863; dispersion =  0.007
# Bacteria CA + S3 + AS vs CA + AS/PHA + AS Pr(>F) = 0.002; R2 = 0.27185; F = 4.4802; dispersion =  0.01869
# Bacteria PHA + G1 + AS vs CA + AS/PHA + AS Pr(>F) = 0.008; R2 = 0.26322; F = 3.5725; dispersion =  0.4655