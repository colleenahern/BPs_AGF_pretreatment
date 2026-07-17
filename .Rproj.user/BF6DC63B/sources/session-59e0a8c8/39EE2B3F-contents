# Colleen Ahern
# 02/24/2026
# Supplemental biogas figures

library(ggplot2)
library(readr)
library(cowplot)

datasup <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/biogas_data/BP_expI_gas_organized_forR_supp_02232026.csv")

# Make growth curves
palette.colors(palette = "Okabe-Ito")

cols <- c("AS"="#0072B2","PLA + AS"="#E69F00","S3 + AS"="#009E73", "PLA + S3 + AS"="#D55E00","tPLA + S3 + AS"="#CC79A7")
p1sup <- ggplot(datasup, aes(x = Day)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averagePLAAS, group = 1, color = "PLA + AS")) + 
  geom_point(aes(y = averagePLAAS, color = "PLA + AS")) + 
  geom_errorbar(aes(ymax = averagePLAAS + sdPLAAS, 
                    ymin = averagePLAAS - sdPLAAS, color = "PLA + AS")) +
  geom_line(aes(y = averageS3AS, group = 1, color = "S3 + AS")) + 
  geom_point(aes(y = averageS3AS, color = "S3 + AS")) + 
  geom_errorbar(aes(ymax = averageS3AS + sdS3AS, 
                    ymin = averageS3AS - sdS3AS, color = "S3 + AS")) +
  geom_line(aes(y = averagePLAS3AS, group = 1, color = "PLA + S3 + AS")) + 
  geom_point(aes(y = averagePLAS3AS, color = "PLA + S3 + AS")) + 
  geom_errorbar(aes(ymax = averagePLAS3AS + sdPLAS3AS, 
                    ymin = averagePLAS3AS - sdPLAS3AS, color = "PLA + S3 + AS")) +
  geom_line(aes(y = averagetPLAS3AS, group = 1, color = "tPLA + S3 + AS")) + 
  geom_point(aes(y = averagetPLAS3AS, color = "tPLA + S3 + AS")) +
  geom_errorbar(aes(ymax = averagetPLAS3AS + sdtPLAS3AS, 
                    ymin = averagetPLAS3AS - sdtPLAS3AS, color = "tPLA + S3 + AS")) +
  labs(
    title = "S3-PLA pairing",
    x = "Day",
    y = "Cumulative Biogas \nProduction (mL)") +
  ylim(0, 250) +
  scale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 12),
         plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
         legend.text = element_text(size=12),
         legend.title = element_text(size=12, face = "bold", vjust = -2))

p1sup

cols <- c("AS"="#0072B2","PLA + AS"="#E69F00","G1 + AS"="#009E73", "PLA + G1 + AS"="#D55E00","tPLA + G1 + AS"="#CC79A7")
p1gsup <- ggplot(datasup, aes(x = Day)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averagePLAAS, group = 1, color = "PLA + AS")) + 
  geom_point(aes(y = averagePLAAS, color = "PLA + AS")) + 
  geom_errorbar(aes(ymax = averagePLAAS + sdPLAAS, 
                    ymin = averagePLAAS - sdPLAAS, color = "PLA + AS")) +
  geom_line(aes(y = averageG1AS, group = 1, color = "G1 + AS")) + 
  geom_point(aes(y = averageG1AS, color = "G1 + AS")) + 
  geom_errorbar(aes(ymax = averageG1AS + sdG1AS, 
                    ymin = averageG1AS - sdG1AS, color = "G1 + AS")) +
  geom_line(aes(y = averagePLAG1AS, group = 1, color = "PLA + G1 + AS")) + 
  geom_point(aes(y = averagePLAG1AS, color = "PLA + G1 + AS")) + 
  geom_errorbar(aes(ymax = averagePLAG1AS + sdPLAG1AS, 
                    ymin = averagePLAG1AS - sdPLAG1AS, color = "PLA + G1 + AS")) +
  geom_line(aes(y = averagetPLAG1AS, group = 1, color = "tPLA + G1 + AS")) + 
  geom_point(aes(y = averagetPLAG1AS, color = "tPLA + G1 + AS")) +
  geom_errorbar(aes(ymax = averagetPLAG1AS + sdtPLAG1AS, 
                    ymin = averagetPLAG1AS - sdtPLAG1AS, color = "tPLA + G1 + AS")) +
  labs(
    title = "G1-PLA pairing",
    x = "Day",
    y = "Cumulative Biogas \nProduction (mL)") +
  ylim(0, 250) +
  scale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 12),
         plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
         legend.text = element_text(size=12),
         legend.title = element_text(size=12, face = "bold", vjust = -2))

p1gsup

cols <- c("AS"="#0072B2","PHA + AS"="#E69F00", "PHA + G1 + AS"="#D55E00","tPHA + G1 + AS"="#CC79A7", "G1 + AS"="#009E73")
q1sup <- ggplot(datasup, aes(x = Day)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averagePHAAS, group = 1, color = "PHA + AS")) + 
  geom_point(aes(y = averagePHAAS, color = "PHA + AS")) + 
  geom_errorbar(aes(ymax = averagePHAAS + sdPHAAS, 
                    ymin = averagePHAAS - sdPHAAS, color = "PHA + AS")) +
  geom_line(aes(y = averagePHAG1AS, group = 1, color = "PHA + G1 + AS")) + 
  geom_point(aes(y = averagePHAG1AS, color = "PHA + G1 + AS")) + 
  geom_errorbar(aes(ymax = averagePHAG1AS + sdPHAG1AS, 
                    ymin = averagePHAG1AS - sdPHAG1AS, color = "PHA + G1 + AS")) +
  geom_line(aes(y = averagetPHAG1AS, group = 1, color = "tPHA + G1 + AS")) + 
  geom_point(aes(y = averagetPHAG1AS, color = "tPHA + G1 + AS")) +
  geom_errorbar(aes(ymax = averagetPHAG1AS + sdtPHAG1AS, 
                    ymin = averagetPHAG1AS - sdtPHAG1AS, color = "tPHA + G1 + AS")) +
  geom_line(aes(y = averageG1AS, group = 1, color = "G1 + AS")) + 
  geom_point(aes(y = averageG1AS, color = "G1 + AS")) + 
  geom_errorbar(aes(ymax = averageG1AS + sdG1AS, 
                    ymin = averageG1AS - sdG1AS, color = "G1 + AS")) +
  labs(
    title = "G1-PHA pairing",
    x = "Day",
    y = "Cumulative Biogas \nProduction (mL)") +
  ylim(0, 250) +
  scale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols, breaks = c("AS", "PHA + AS", "PHA + G1 + AS", "G1 + AS", "tPHA + AS", "tPHA + G1 + AS")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2))

q1sup

cols <- c("AS"="#0072B2","PHA + AS"="#E69F00", "PHA + S3 + AS"="#D55E00","tPHA + S3 + AS"="#CC79A7", "S3 + AS"="#009E73")
q1ssup <- ggplot(datasup, aes(x = Day)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averagePHAAS, group = 1, color = "PHA + AS")) + 
  geom_point(aes(y = averagePHAAS, color = "PHA + AS")) + 
  geom_errorbar(aes(ymax = averagePHAAS + sdPHAAS, 
                    ymin = averagePHAAS - sdPHAAS, color = "PHA + AS")) +
  geom_line(aes(y = averagePHAS3AS, group = 1, color = "PHA + S3 + AS")) + 
  geom_point(aes(y = averagePHAS3AS, color = "PHA + S3 + AS")) + 
  geom_errorbar(aes(ymax = averagePHAS3AS + sdPHAS3AS, 
                    ymin = averagePHAS3AS - sdPHAS3AS, color = "PHA + S3 + AS")) +
  geom_line(aes(y = averagetPHAS3AS, group = 1, color = "tPHA + S3 + AS")) + 
  geom_point(aes(y = averagetPHAS3AS, color = "tPHA + S3 + AS")) +
  geom_errorbar(aes(ymax = averagetPHAS3AS + sdtPHAS3AS, 
                    ymin = averagetPHAS3AS - sdtPHAS3AS, color = "tPHA + S3 + AS")) +
  geom_line(aes(y = averageS3AS, group = 1, color = "S3 + AS")) + 
  geom_point(aes(y = averageS3AS, color = "S3 + AS")) + 
  geom_errorbar(aes(ymax = averageS3AS + sdS3AS, 
                    ymin = averageS3AS - sdS3AS, color = "S3 + AS")) +
  labs(
    title = "S3-PHA pairing",
    x = "Day",
    y = "Cumulative Biogas \nProduction (mL)") +
  ylim(0, 250) +
  scale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols, breaks = c("AS", "PHA + AS", "PHA + S3 + AS", "S3 + AS", "tPHA + AS", "tPHA + S3 + AS")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2))

q1ssup

cols <- c("AS"="#0072B2","CA + AS"="#E69F00","S3 + AS"="#009E73", "CA + S3 + AS"="#D55E00","tCA + S3 + AS"="#CC79A7")
r1sup <- ggplot(datasup, aes(x = Day)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averageCAAS, group = 1, color = "CA + AS")) + 
  geom_point(aes(y = averageCAAS, color = "CA + AS")) + 
  geom_errorbar(aes(ymax = averageCAAS + sdCAAS, 
                    ymin = averageCAAS - sdCAAS, color = "CA + AS")) +
  geom_line(aes(y = averageS3AS, group = 1, color = "S3 + AS")) + 
  geom_point(aes(y = averageS3AS, color = "S3 + AS")) + 
  geom_errorbar(aes(ymax = averageS3AS + sdS3AS, 
                    ymin = averageS3AS - sdS3AS, color = "S3 + AS")) +
  geom_line(aes(y = averageCAS3AS, group = 1, color = "CA + S3 + AS")) + 
  geom_point(aes(y = averageCAS3AS, color = "CA + S3 + AS")) + 
  geom_errorbar(aes(ymax = averageCAS3AS + sdCAS3AS, 
                    ymin = averageCAS3AS - sdCAS3AS, color = "CA + S3 + AS")) +
  geom_line(aes(y = averagetCAS3AS, group = 1, color = "tCA + S3 + AS")) + 
  geom_point(aes(y = averagetCAS3AS, color = "tCA + S3 + AS")) +
  geom_errorbar(aes(ymax = averagetCAS3AS + sdtCAS3AS, 
                    ymin = averagetCAS3AS - sdtCAS3AS, color = "tCA + S3 + AS")) +
  labs(
    title = "S3-CA pairing",
    x = "Day",
    y = "Cumulative Biogas \nProduction (mL)") +
  ylim(0, 250) +
  scale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2))
r1sup

cols <- c("AS"="#0072B2","CA + AS"="#E69F00","G1 + AS"="#009E73", "CA + G1 + AS"="#D55E00","tCA + G1 + AS"="#CC79A7")
r1gsup <- ggplot(datasup, aes(x = Day)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averageCAAS, group = 1, color = "CA + AS")) + 
  geom_point(aes(y = averageCAAS, color = "CA + AS")) + 
  geom_errorbar(aes(ymax = averageCAAS + sdCAAS, 
                    ymin = averageCAAS - sdCAAS, color = "CA + AS")) +
  geom_line(aes(y = averageG1AS, group = 1, color = "G1 + AS")) + 
  geom_point(aes(y = averageG1AS, color = "G1 + AS")) + 
  geom_errorbar(aes(ymax = averageG1AS + sdG1AS, 
                    ymin = averageG1AS - sdG1AS, color = "G1 + AS")) +
  geom_line(aes(y = averageCAG1AS, group = 1, color = "CA + G1 + AS")) + 
  geom_point(aes(y = averageCAG1AS, color = "CA + G1 + AS")) + 
  geom_errorbar(aes(ymax = averageCAG1AS + sdCAG1AS, 
                    ymin = averageCAG1AS - sdCAG1AS, color = "CA + G1 + AS")) +
  geom_line(aes(y = averagetCAG1AS, group = 1, color = "tCA + G1 + AS")) + 
  geom_point(aes(y = averagetCAG1AS, color = "tCA + G1 + AS")) +
  geom_errorbar(aes(ymax = averagetCAG1AS + sdtCAG1AS, 
                    ymin = averagetCAG1AS - sdtCAG1AS, color = "tCA + G1 + AS")) +
  labs(
    title = "G1-CA pairing",
    x = "Day",
    y = "Cumulative Biogas \nProduction (mL)") +
  ylim(0, 250) +
  scale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2))
r1gsup

# Combine
tp1sup <- plot_grid(p1sup, p1gsup, r1sup, r1gsup, q1ssup, q1sup, nrow = 3, ncol = 2, labels = c("A)", "B)", "C)", "D)", "E)", "F)"))
tp1sup

ggsave(plot = tp1sup, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Biogas_suppfig_07152026.pdf", width = 10, height = 10, units = "in", dpi = 300)
ggsave(plot = tp1sup, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Biogas_suppfig_07152026.tiff", width = 10, height = 10, units = "in", dpi = 200, bg="white")
