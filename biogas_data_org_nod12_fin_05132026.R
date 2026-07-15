# =============================================================================
# Biogas Data — Line Plots and Bar Plots with Significance
# Colleen Ahern | Last updated: 05/13/2026
#
# Notes:
#   - Day 12 removed from CA and PLA samples due to bad S3 + CA data
#   - PHA uses G1 not S3 so Day 12 is kept for PHA
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(ggpubr)
library(cowplot)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    <- "/Users/colleenahern/Documents/Magda_BPs_experiment"
FIG_DIR <- file.path(BASE, "Figures")
DATE    <- "05132026"

# ── Shared theme ──────────────────────────────────────────────────────────────
biogas_line_theme <- theme_classic() +
  theme(
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 12),
    plot.title   = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold", vjust = -2)
  )

biogas_bar_theme <- theme_classic() +
  theme(
    axis.text.y  = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(size = 12),
    plot.title   = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# ── Vertical lines for timepoints ─────────────────────────────────────────────
timepoint_lines <- list(
  geom_vline(xintercept = 21, linetype = "dashed", color = "grey50", linewidth = 0.4),
  geom_vline(xintercept = 27, linetype = "dashed", color = "grey50", linewidth = 0.4),
  geom_vline(xintercept = 31, linetype = "dashed", color = "grey50", linewidth = 0.4),
  geom_vline(xintercept = 47, linetype = "dashed", color = "grey50", linewidth = 0.4),
  annotate("text", x = 21, y = -Inf, label = "t1", size = 10 / .pt, color = "grey50",
           hjust = -0.2, vjust = -0.5),
  annotate("text", x = 27, y = -Inf, label = "t2", size = 10 / .pt, color = "grey50",
           hjust = -0.2, vjust = -0.5),
  annotate("text", x = 31, y = -Inf, label = "t3", size = 10 / .pt, color = "grey50",
           hjust = -0.2, vjust = -0.5),
  annotate("text", x = 47, y = -Inf, label = "t4", size = 10 / .pt, color = "grey50",
           hjust = -0.2, vjust = -0.5)
)

# ── Helper: build one line + errorbar layer ────────────────────────────────────
biogas_layer <- function(avg_col, sd_col, color_label, width = 1.5) {
  list(
    geom_line(aes(y = .data[[avg_col]], color = color_label, group = 1)),
    geom_point(aes(y = .data[[avg_col]], color = color_label)),
    geom_errorbar(aes(ymax = .data[[avg_col]] + .data[[sd_col]],
                      ymin = .data[[avg_col]] - .data[[sd_col]],
                      color = color_label), width = width)
  )
}

# =============================================================================
# 1.  Load data
# =============================================================================

# PHA data (all days)
data <- read_csv(file.path(BASE, "biogas_data/BP_expII_gas_organized_forR_01192026.csv")) %>%
  rename_with(~ gsub("PLA_70C", "tPLA", .) %>%
                gsub("PHA_70C", "tPHA", .) %>%
                gsub("CA_70C",  "tCA",  .))

# CA and PLA data (no Day 12)
datanod12 <- read_csv(file.path(BASE, "biogas_data/BP_expII_gas_organized_forR_noD12_CAPLA_04162026.csv")) %>%
  dplyr::select(-matches("ind"))

# =============================================================================
# 2.  PHA line plot
# =============================================================================

cols_pha <- c("AS" = "#0072B2", "PHA + AS" = "#E69F00", "tPHA + AS" = "#56B4E9",
              "PHA + G1 + AS" = "#D55E00", "tPHA + G1 + AS" = "#CC79A7", "G1 + AS" = "#009E73")

phalin <- ggplot(data, aes(x = Days)) +
  biogas_layer("averageAS",       "sdAS",          "AS") +
  biogas_layer("averagePHAAS",    "sdPHAAS",       "PHA + AS") +
  biogas_layer("averagePHA70CAS", "sdPHA70CAS",    "tPHA + AS") +
  biogas_layer("averagePHAG1AS",  "sdPHAG1AS",     "PHA + G1 + AS") +
  biogas_layer("averagePHA70CG1AS", "sdPHA70CG1AS","tPHA + G1 + AS") +
  biogas_layer("averageG1AS",     "sdG1AS",        "G1 + AS") +
  timepoint_lines +
  scale_color_manual(name = "Condition", values = cols_pha,
                     breaks = c("AS", "PHA + AS", "PHA + G1 + AS", "G1 + AS", "tPHA + AS", "tPHA + G1 + AS")) +
  scale_x_continuous(breaks = c(4, 32, 60), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 250)) +
  labs(title = "PHA", x = "Day", y = "Cumulative Biogas Production (mL)") +
  biogas_line_theme

phalin
# =============================================================================
# 3.  PHA bar plot
# =============================================================================

data_subphag1 <- data[13, grepl("PHA \\+ AS|PHA \\+ G1 \\+ AS|tPHA \\+ G1 \\+ AS", colnames(data))]
data_subphag1 <- cbind(data[13, c(39:41, 2:4)], data_subphag1) %>%
  t() %>% as.data.frame() %>%
  mutate(group = substr(rownames(.), 1, nchar(rownames(.)) - 2),
         group = factor(group, levels = c("PHA + G1 + AS", "PHA + AS", "G1 + AS",
                                          "AS", "tPHA + AS", "tPHA + G1 + AS")))

phabar <- ggbarplot(data_subphag1, x = "group", y = "V1",
                    add = "mean_se", fill = "group",
                    ylab = "Total Biogas Production (mL)") +
  labs(x = NULL, title = "PHA") +
  ylim(0, 400) +
  scale_fill_manual(values = c("AS" = "#0072B2", "PHA + AS" = "#E69F00", "tPHA + AS" = "#56B4E9",
                               "G1 + AS" = "#009E73", "PHA + G1 + AS" = "#D55E00", "tPHA + G1 + AS" = "#CC79A7")) +
  scale_x_discrete(labels = c("AS" = "AS", "PHA + AS" = "PHA +\nAS", "tPHA + AS" = "tPHA +\nAS",
                              "G1 + AS" = "G1 +\nAS", "PHA + G1 + AS" = "PHA +\nG1 +\nAS",
                              "tPHA + G1 + AS" = "tPHA +\nG1 +\nAS")) +
  stat_compare_means(method = "t.test", label = "p.signif",
                     comparisons = list(c("PHA + G1 + AS", "PHA + AS"), c("PHA + G1 + AS", "G1 + AS"),
                                        c("PHA + G1 + AS", "AS"), c("tPHA + G1 + AS", "tPHA + AS"),
                                        c("tPHA + G1 + AS", "AS"), c("tPHA + G1 + AS", "G1 + AS")),
                     label.y = c(240, 270, 300, 330, 360, 390)) +
  biogas_bar_theme

phabar
# =============================================================================
# 4.  CA line plot (replicates 1, 2, 3)
# =============================================================================

cols_ca <- c("AS" = "#0072B2", "CA + AS" = "#E69F00", "tCA + AS" = "#56B4E9",
             "S3 + AS" = "#009E73", "CA + S3 + AS" = "#D55E00", "tCA + S3 + AS" = "#CC79A7")

calin <- ggplot(datanod12, aes(x = Day)) +
  biogas_layer("averageAS",        "sdAS",         "AS") +
  biogas_layer("averageCAAS",      "sdCAAS",       "CA + AS") +
  biogas_layer("averagetCAAS",     "sdtCAAS",      "tCA + AS") +
  biogas_layer("averageS3AS",      "sdS3AS",       "S3 + AS") +
  biogas_layer("averageCAS3AS123", "sdCAS3AS123",  "CA + S3 + AS") +
  biogas_layer("averagetCAS3AS",   "sdtCAS3AS",    "tCA + S3 + AS") +
  timepoint_lines +
  scale_color_manual(name = "Condition", values = cols_ca) +
  scale_x_continuous(breaks = c(4, 32, 60), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 250)) +
  labs(title = "CA", x = "Day", y = "Cumulative Biogas Production (mL)") +
  biogas_line_theme

calin
# =============================================================================
# 5.  CA bar plot (replicates 1, 2, 3)
# =============================================================================

data_subcas3 <- datanod12[12, grepl("CA \\+ AS|CA \\+ S3 \\+ AS|tCA \\+ S3 \\+ AS", colnames(datanod12))]
data_subcas3 <- cbind(datanod12[12, c(2, 3, 6:8)], data_subcas3) %>%
  t() %>% as.data.frame() %>%
  mutate(group = substr(rownames(.), 1, nchar(rownames(.)) - 2),
         group = factor(group, levels = c("CA + S3 + AS", "CA + AS", "S3 + AS",
                                          "AS", "tCA + AS", "tCA + S3 + AS")))

cabar <- ggbarplot(data_subcas3, x = "group", y = "V1",
                   add = "mean_se", fill = "group",
                   ylab = "Total Biogas Production (mL)") +
  labs(x = NULL, title = "CA") +
  ylim(0, 400) +
  scale_fill_manual(values = c("AS" = "#0072B2", "CA + AS" = "#E69F00", "tCA + AS" = "#56B4E9",
                               "S3 + AS" = "#009E73", "CA + S3 + AS" = "#D55E00", "tCA + S3 + AS" = "#CC79A7")) +
  scale_x_discrete(labels = c("AS" = "AS", "CA + AS" = "CA +\nAS", "tCA + AS" = "tCA +\nAS",
                              "S3 + AS" = "S3 +\nAS", "CA + S3 + AS" = "CA +\nS3 +\nAS",
                              "tCA + S3 + AS" = "tCA +\nS3 +\nAS")) +
  stat_compare_means(method = "t.test", label = "p.signif",
                     comparisons = list(c("CA + S3 + AS", "CA + AS"), c("CA + S3 + AS", "S3 + AS"),
                                        c("CA + S3 + AS", "AS"), c("tCA + S3 + AS", "tCA + AS"),
                                        c("tCA + S3 + AS", "AS"), c("tCA + S3 + AS", "S3 + AS")),
                     label.y = c(240, 270, 300, 330, 360, 390)) +
  biogas_bar_theme

cabar
# =============================================================================
# 6.  PLA line plot
# =============================================================================

cols_pla <- c("AS" = "#0072B2", "PLA + AS" = "#E69F00", "tPLA + AS" = "#56B4E9",
              "S3 + AS" = "#009E73", "PLA + S3 + AS" = "#D55E00", "tPLA + S3 + AS" = "#CC79A7")

plalin <- ggplot(datanod12, aes(x = Day)) +
  biogas_layer("averageAS",      "sdAS",       "AS") +
  biogas_layer("averagePLAAS",   "sdPLAAS",    "PLA + AS") +
  biogas_layer("averagetPLAAS",  "sdtPLAAS",   "tPLA + AS") +
  biogas_layer("averageS3AS",    "sdS3AS",     "S3 + AS") +
  biogas_layer("averagePLAS3AS", "sdPLAS3AS",  "PLA + S3 + AS") +
  biogas_layer("averagetPLAS3AS","sdtPLAS3AS", "tPLA + S3 + AS") +
  timepoint_lines +
  scale_color_manual(name = "Condition", values = cols_pla) +
  scale_x_continuous(breaks = c(4, 32, 60), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 250)) +
  labs(title = "PLA", x = "Day", y = "Cumulative Biogas Production (mL)") +
  biogas_line_theme

plalin
# =============================================================================
# 7.  PLA bar plot
# =============================================================================

data_subplas3 <- datanod12[12, grepl("PLA \\+ AS|PLA \\+ S3 \\+ AS|tPLA \\+ S3 \\+ AS", colnames(datanod12))]
data_subplas3 <- cbind(datanod12[12, c(2, 3, 6:8)], data_subplas3) %>%
  t() %>% as.data.frame() %>%
  mutate(group = substr(rownames(.), 1, nchar(rownames(.)) - 2),
         group = factor(group, levels = c("PLA + S3 + AS", "PLA + AS", "S3 + AS",
                                          "AS", "tPLA + AS", "tPLA + S3 + AS")))

plabar <- ggbarplot(data_subplas3, x = "group", y = "V1",
                    add = "mean_se", fill = "group",
                    ylab = "Total Biogas Production (mL)") +
  labs(x = NULL, title = "PLA") +
  ylim(0, 400) +
  scale_fill_manual(values = c("AS" = "#0072B2", "PLA + AS" = "#E69F00", "tPLA + AS" = "#56B4E9",
                               "S3 + AS" = "#009E73", "PLA + S3 + AS" = "#D55E00", "tPLA + S3 + AS" = "#CC79A7")) +
  scale_x_discrete(labels = c("AS" = "AS", "PLA + AS" = "PLA +\nAS", "tPLA + AS" = "tPLA +\nAS",
                              "S3 + AS" = "S3 +\nAS", "PLA + S3 + AS" = "PLA +\nS3 +\nAS",
                              "tPLA + S3 + AS" = "tPLA +\nS3 +\nAS")) +
  stat_compare_means(method = "t.test", label = "p.signif",
                     comparisons = list(c("PLA + S3 + AS", "PLA + AS"), c("PLA + S3 + AS", "S3 + AS"),
                                        c("PLA + S3 + AS", "AS"), c("tPLA + S3 + AS", "tPLA + AS"),
                                        c("tPLA + S3 + AS", "AS"), c("tPLA + S3 + AS", "S3 + AS")),
                     label.y = c(240, 270, 300, 330, 360, 390)) +
  biogas_bar_theme

plabar

# ============================================================
# 8. % Difference Calculations
# ============================================================

means_pha <- tapply(data_subphag1$V1, data_subphag1$group, mean, na.rm = TRUE)
means_ca  <- tapply(data_subcas3$V1,  data_subcas3$group,  mean, na.rm = TRUE)
means_pla <- tapply(data_subplas3$V1, data_subplas3$group, mean, na.rm = TRUE)

pct_diff <- function(a, b) round((a - b) / b * 100, 1)

# --- PHA ---
pha_results <- data.frame(
  Comparison = c("PHA+G1+AS vs PHA+AS", "PHA+G1+AS vs G1+AS", "PHA+G1+AS vs AS",
                 "tPHA+G1+AS vs tPHA+AS", "tPHA+G1+AS vs AS", "tPHA+G1+AS vs G1+AS"),
  Pct_Diff = c(
    pct_diff(means_pha["PHA + G1 + AS"],  means_pha["PHA + AS"]),
    pct_diff(means_pha["PHA + G1 + AS"],  means_pha["G1 + AS"]),
    pct_diff(means_pha["PHA + G1 + AS"],  means_pha["AS"]),
    pct_diff(means_pha["tPHA + G1 + AS"], means_pha["tPHA + AS"]),
    pct_diff(means_pha["tPHA + G1 + AS"], means_pha["AS"]),
    pct_diff(means_pha["tPHA + G1 + AS"], means_pha["G1 + AS"])
  )
)

# --- CA ---
ca_results <- data.frame(
  Comparison = c("CA+S3+AS vs CA+AS", "CA+S3+AS vs S3+AS", "CA+S3+AS vs AS",
                 "tCA+S3+AS vs tCA+AS", "tCA+S3+AS vs AS", "tCA+S3+AS vs S3+AS"),
  Pct_Diff = c(
    pct_diff(means_ca["CA + S3 + AS"],  means_ca["CA + AS"]),
    pct_diff(means_ca["CA + S3 + AS"],  means_ca["S3 + AS"]),
    pct_diff(means_ca["CA + S3 + AS"],  means_ca["AS"]),
    pct_diff(means_ca["tCA + S3 + AS"], means_ca["tCA + AS"]),
    pct_diff(means_ca["tCA + S3 + AS"], means_ca["AS"]),
    pct_diff(means_ca["tCA + S3 + AS"], means_ca["S3 + AS"])
  )
)

# --- PLA ---
pla_results <- data.frame(
  Comparison = c("PLA+S3+AS vs PLA+AS", "PLA+S3+AS vs S3+AS", "PLA+S3+AS vs AS",
                 "tPLA+S3+AS vs tPLA+AS", "tPLA+S3+AS vs AS", "tPLA+S3+AS vs S3+AS"),
  Pct_Diff = c(
    pct_diff(means_pla["PLA + S3 + AS"],  means_pla["PLA + AS"]),
    pct_diff(means_pla["PLA + S3 + AS"],  means_pla["S3 + AS"]),
    pct_diff(means_pla["PLA + S3 + AS"],  means_pla["AS"]),
    pct_diff(means_pla["tPLA + S3 + AS"], means_pla["tPLA + AS"]),
    pct_diff(means_pla["tPLA + S3 + AS"], means_pla["AS"]),
    pct_diff(means_pla["tPLA + S3 + AS"], means_pla["S3 + AS"])
  )
)

# --- Print results ---
cat("=== PHA ===\n");  print(pha_results, row.names = FALSE)
cat("\n=== CA ===\n"); print(ca_results,  row.names = FALSE)
cat("\n=== PLA ===\n"); print(pla_results, row.names = FALSE)

# =============================================================================
# 9. Combined figures
# =============================================================================

# ── Version 1: CA replicates 1, 2, 3 ─────────────────────────────────────────
finbiogas <- plot_grid(phalin, phabar, calin, cabar, plalin, plabar,
                       nrow = 3, ncol = 2,
                       labels = c("A)", "B)", "C)", "D)", "E)", "F)"))
finbiogas

ggsave(plot     = finbiogas,
       filename = file.path(FIG_DIR, paste0("Biogas_fig_noD12CAPLA_123_", DATE, ".pdf")),
       width = 10, height = 10, units = "in", dpi = 300, bg = "white")

ggsave(plot     = finbiogas,
       filename = file.path(FIG_DIR, paste0("Biogas_fig_noD12CAPLA_123_", DATE, ".tiff")),
       width = 10, height = 10, units = "in", dpi = 200, bg = "white")

message("Done. Figures written to: ", FIG_DIR)