# Colleen Ahern
# Final manuscript heatmap from my S3CA DESeq2 analysis 

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("tximport")
# BiocManager::install("tximportData")
# BiocManager::install("pasilla")

# NOTE: skip to line 383 for saved vary_43caz table for heatmap making

library(tximport)
library(readr)
library(tximportData)
library(tidyverse)
library(pasilla)

# Heatmap for all genes

# Import the count data and sample info
s3seqdata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/featcount/s3ca_countmat_07302025.txt", stringsAsFactors = FALSE)
s3sampleinfo <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/s3sampleinfo.txt", stringsAsFactors = TRUE)
head(s3seqdata)
dim(s3seqdata)

s3sampleinfo

# Store GeneID as rownames and then remove non-count columns
s3countdata <- s3seqdata
rownames(s3countdata) <- s3countdata$Geneid
s3countdata <- s3countdata[,-c(1:6)]
head(s3countdata)

# Make sure that the column names are the same as libraryName in the s3sampleinfo file
colnames(s3countdata) # using substr, you can change the characters of the colnames
all(colnames(s3countdata) %in% s3sampleinfo$libraryName)
all(colnames(s3countdata)==s3sampleinfo$libraryName)
# s3countdata <- s3countdata[,as.character(s3sampleinfo$libraryName)]
# table(colnames(s3countdata)==s3sampleinfo$libraryName) # Now they match
colnames(s3countdata) <- s3sampleinfo$sampleName

# my S3CAAS samples suck so remove those, also remove S3_5 cus it's an outlier
s3countdata <- s3countdata[,c(1:4,6:10)]
s3sampleinfo <- s3sampleinfo[c(1:4,6:10),]

# Apply tpm cutoff to the gene count matrix
# Make tpm gene count matrix
# install.packages("remotes")
# remotes::install_github("davidrequena/drfun")
library(DRnaSeq)
all(rownames(s3countdata) == s3seqdata$Geneid)
gl <- s3seqdata$Length
tpm <- as.data.frame(tpm(s3countdata, gl))

# Import/reformat sample information 
asm <- s3sampleinfo
all(colnames(tpm) == asm$sampleName)

# Make sure the sample names in the count matrix and sample info match
colnames(tpm) == asm$sampleName

tpm$rem <- 0

# Apply tpm cutoff: the average TPM of the  biological replicates for any of the four substrates must be >2
for (i in 1:nrow(tpm)) {
  tpm$rem[i] <- ifelse(mean(c(tpm$S3_1[i],tpm$S3_2[i],tpm$S3_3[i],tpm$S3_4[i])) < 2 && mean(c(tpm$S3CA_1[i],tpm$S3CA_2[i],tpm$S3CA_3[i],tpm$S3CA_4[i],tpm$S3CA_5[i])) < 2, "REMOVE", "KEEP")
}

tpm <- tpm[(tpm$rem=="KEEP"),] # reduces from 27677 genes to 13515
tpm <- tpm[,-c(10)]
tpm$GeneID <- rownames(tpm)

# Import list of differentially regulated genes for heatmap
# Remove genes with a p-value = NA for all three comparitive conditions

vary_43 <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/s3_DESeq_ALL_08222025.csv") # My DESeq
vary_43 <- vary_43[,c(8,2,6,7)]

vary_43 <- vary_43[is.na(vary_43$S3CA_S3_padj) == FALSE,] # reduces from 25498 to 22493

# Merge tpm matrix with vary_43 matrix - to make a matrix that contains only genes that pass both of those filters
a <- merge(vary_43,tpm,"GeneID") #13501 genes 
vary_43 <- a[,c(1:4)]

# Sort the matrix from lowest padj to largest
vary_43 <- vary_43[order(vary_43$S3CA_S3_padj,decreasing=FALSE),]
vary_43_sig <- vary_43[vary_43$S3CA_S3_padj < 0.05,]

# Import list of JGI gene annotations
# BiocManager::install("ape")
library(ape)
library(tidyverse)
my_gene_col_44 <- read.gff("/Users/colleenahern/Documents/Magda_BPs_experiment/S3_JGI_files/Neolan1_GeneCatalog_20200610.gff3", GFF3 = TRUE)
my_gene_col_44 <- my_gene_col_44[my_gene_col_44$type == "gene",] # 27677 genes
my_gene_col_44b <- my_gene_col_44 %>%
  dplyr::select(seqid, everything())%>%
  mutate(attributes=str_split(attributes, ";")) %>% 
  unnest_wider(where(is.list), names_sep = "")

my_gene_col_44 <- my_gene_col_44b[,c(9,12,13,14)]
my_gene_col_44
my_gene_col_44$attributes1 <- gsub("ID=","",my_gene_col_44$attributes1)
my_gene_col_44$attributes4 <- gsub("product_name=","",my_gene_col_44$attributes4)
my_gene_col_44$attributes5 <- gsub("proteinId=","",my_gene_col_44$attributes5)
my_gene_col_44$attributes6 <- gsub("transcriptId=","",my_gene_col_44$attributes6)
colnames(my_gene_col_44) <- c("GeneID","Protein_name","ProteinID","TranscriptID")

write_csv(vary_43, "/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/Tables/s3_vary43_s3cavs3_beforeanno_08132025.csv")
# Add annotations to the data frame and reorder padjmin
vary_43 <- merge(vary_43, my_gene_col_44, "GeneID")
vary_43 <- vary_43[order(vary_43$S3CA_S3_padj,decreasing=FALSE),]

for (i in 1:nrow(vary_43)) {
  vary_43$Name[i] <- paste(vary_43$Protein_name[i], " ", "(proteinID = ", vary_43$ProteinID[i], ")", sep = "")
}

# write_csv(vary_43, "/Users/colleenahern/Documents/s3hry/DESeq2_analysis/Tables/vary_43_07242025.csv")
# vary_43 <- read_csv("/Users/colleenahern/Documents/s3hry/DESeq2_analysis/Tables/vary_43_07242025.csv")
vary_43 <- vary_43[order(vary_43$S3CA_S3_padj,decreasing=FALSE),]
write_csv(vary_43, "/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/Tables/s3_vary43_s3cavs3_anno_08132025.csv")
vary_43 <- read.csv("/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/Tables/s3_vary43_s3cavs3_anno_08132025.csv")

# 08/14/2025 - import list of CAZymes from dbcan
cazlist <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/dbcan/overview.txt")
cazlist2 <- cazlist[cazlist$X.ofTools >1,]
cazlist_f <- cazlist2
cazlist_f$Gene.ID <- gsub("jgi\\|Neolan1\\|","",cazlist_f$Gene.ID)
cazlist_f$Gene.ID <- gsub("\\|.*", "", cazlist_f$Gene.ID)

cazlist_f$DIAMOND <- gsub(".*fasta\\+","",cazlist_f$DIAMOND)
# cazlist_f$HMMER <- gsub("\\(.*\\)", "", cazlist_f$HMMER)

rownames(cazlist_f) <- 1:nrow(cazlist_f)

cazlist_ff <- cazlist_f %>%
  dplyr::select(Gene.ID, everything())%>%
  mutate(HMMER=str_split(HMMER, "\\+")) %>% 
  unnest_wider(where(is.list), names_sep = "")

for (i in 1:nrow(cazlist_ff)) {
  for (j in 3:9)
    cazlist_ff[i,j] <- gsub("\\(.*\\)", "", cazlist_ff[i,j])
}

cazlist_ff$HMMER_u <- NA
for (i in 1:nrow(cazlist_ff)) {
    gg <- unique(as.vector((as.matrix(cazlist_ff[i,3:9]))))
    gg <- gg[!is.na(gg)]
    for (j in 1:length(gg)) {
      cazlist_ff$HMMER_u[i] <- paste(cazlist_ff$HMMER_u[i],gg[j], sep = "+")
    }
}
cazlist_ff$HMMER_u <- gsub("NA\\+","",cazlist_ff$HMMER_u)
cazlist_ff <- cazlist_ff[,c(1:9,13,10:12)]

cazlist_ff <- cazlist_ff %>%
  dplyr::select(Gene.ID, everything())%>%
  mutate(HMMER_u=str_split(HMMER_u, "\\+")) %>% 
  unnest_wider(where(is.list), names_sep = "")
cazlist_ff <- cazlist_ff %>%
  dplyr::select(Gene.ID, everything())%>%
  mutate(eCAMI=str_split(eCAMI, "\\+")) %>% 
  unnest_wider(where(is.list), names_sep = "")
cazlist_ff <- cazlist_ff %>%
  dplyr::select(Gene.ID, everything())%>%
  mutate(DIAMOND=str_split(DIAMOND, "\\+")) %>% 
  unnest_wider(where(is.list), names_sep = "")


cazlist_fff <- cazlist_ff
cazlist_fff$CAZ <- NA

for (i in 1:nrow(cazlist_fff)) {
  z <- as.data.frame(table(as.character(cazlist_fff[i,10:22])))
  cazlist_fff$CAZ[i] <- paste(paste(z$Var1," x",z$Freq,sep = ""),collapse = "; ")
}

cazlist_f4 <- cazlist_fff
cazlist_f4$CAZf <- cazlist_f4$CAZ

cazlist_f4 <- cazlist_f4 %>%
  dplyr::select(Gene.ID, everything())%>%
  mutate(CAZf=str_split(CAZf, "; ")) %>% 
  unnest_wider(where(is.list), names_sep = "")

for (i in 25:32) {
  cazlist_f4[,i] <- lapply(cazlist_f4[,i], function(x) replace(x, grepl("x1", x), NA))
}

cazlist_f4 <- cazlist_f4 %>%
  unite(CAZ_name, starts_with('CAZf'), 
        na.rm = TRUE, remove = TRUE, sep = '; ')

cazlist_f4$CAZ_anno <- cazlist_f4$CAZ_name
cazlist_f4$CAZ_anno <- gsub(" x2","",cazlist_f4$CAZ_anno)
cazlist_f4$CAZ_anno <- gsub(" x3","",cazlist_f4$CAZ_anno)

cazlist_f5 <- cazlist_f4
cazlist_f5 <- cazlist_f5[!cazlist_f5$CAZ_anno == "",]
cazlist_f5 <- cazlist_f5 %>%
  rename("ProteinID" = "Gene.ID")

write.csv(cazlist_f5, "/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/dbcan/cazlist_08152025.csv")
caz_anno <- read.csv("/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/dbcan/cazlist_08152025.csv")
caz_anno <- caz_anno[,c("ProteinID","CAZ_anno")]

vary_43caz <- merge(vary_43,caz_anno, by = "ProteinID", all.x = TRUE)
write_csv(vary_43caz, "/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/Tables/s3_vary43_s3cavss3_anno_cazanno_08152025.csv")
vary_43caz <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/Tables/s3_vary43_s3cavss3_anno_cazanno_08152025.csv")

vary_43cazonly <- vary_43caz[is.na(vary_43caz$CAZ_anno) == FALSE,] 
write_csv(vary_43cazonly, "/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/Tables/s3_vary43_s3cavss3_anno_cazonly_08152025.csv")


# look at the significantly upregulated (and maybe downregulated) genes in the vary_4c3caz and vary_43cazonly datasets
vary_43caz_sig <- vary_43caz[vary_43caz$S3CA_S3_padj < 0.05,]
vary_43caz_sig <- vary_43caz_sig[order(vary_43caz_sig$S3CA_S3_padj,decreasing=FALSE),]

vary_43caz_sig_poslfc <- vary_43caz_sig[vary_43caz_sig$S3CA_S3_log2FoldChange > 0,]
vary_43caz_sig_poslfc <- vary_43caz_sig_poslfc[order(vary_43caz_sig_poslfc$S3CA_S3_padj,decreasing=FALSE),]

vary_43caz_sig_neglfc <- vary_43caz_sig[vary_43caz_sig$S3CA_S3_log2FoldChange < 0,]
vary_43caz_sig_neglfc <- vary_43caz_sig_neglfc[order(vary_43caz_sig_neglfc$S3CA_S3_padj,decreasing=FALSE),]

length(vary_43caz$CAZ_anno[is.na(vary_43caz$CAZ_anno) == FALSE])
length(vary_43caz_sig$CAZ_anno[is.na(vary_43caz_sig$CAZ_anno) == FALSE])
length(vary_43caz_sig_poslfc$CAZ_anno[is.na(vary_43caz_sig_poslfc$CAZ_anno) == FALSE])
length(vary_43caz_sig_neglfc$CAZ_anno[is.na(vary_43caz_sig_neglfc$CAZ_anno) == FALSE])

unique(vary_43caz_sig_poslfc$CAZ_anno)


# Make heatmap 
# newnames <- lapply(
#   colnames(vary_43[,c(2,4,6)]),
#   function(x) bquote(.(x)))
# newnames
# head(vary_43)
# nrow(vary_43)
# rownum_from <- 1
# rownum_to <- 100
# 
# vary_44 <- vary_43[c(rownum_from:rownum_to),]
# 
# ann_cols <- vary_44[,12]
# #install.packages("RColorBrewer")
# library("RColorBrewer")
# display.brewer.pal(n = 12, name = 'Set3')
# brewer.pal(n=12, name="Set3")
# display.brewer.pal(n = 8, name = 'Set2')
# brewer.pal(n=8, name="Set2")
# display.brewer.pal(n = 9, name = 'Pastel1')
# brewer.pal(n=9, name="Pastel1")
# display.brewer.pal(n = 8, name = 'Pastel2')
# brewer.pal(n=8, name="Pastel2")
# display.brewer.pal(n = 8, name = 'Dark2')
# brewer.pal(n=8, name="Dark2")
# display.brewer.pal(n = 12, name = 'Paired')
# brewer.pal(n=12, name="Paired")

# #install.packages("pheatmap")
# library("pheatmap")
# library(grid)
# draw_colnames_45 <- function (coln, gaps, ...) {
#   coord = pheatmap:::find_coordinates(length(coln), gaps)
#   x = coord$coord - 0.5 * coord$size
#   res = textGrob(coln, x = x, y = unit(1, "ns3") - unit(3,"bigpts"), vjust = 0.8, hjust = .5, rot = 360, gp = gpar(...))
#   return(res)}
# assignInNamespace(x="draw_colnames", value="draw_colnames_45",
#                   ns=asNamespace("pheatmap"))
# 
# my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)
# 
# s3_regulation_heatmap<-pheatmap(data.matrix(vary_44[,c(2,4,6)]), col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
#                                 fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(vary_44[,12]), border_color=NA,breaks=seq(-4, 6, length.out=300))
# 



# Make volcano plot instead - 11/20/2025
# making volcano plots instead of heatmmaps to show differentially expressed genes (and maybe just CAZymes)
# taken from https://biostatsquid.com/volcano-plots-r-tutorial/

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# Import annotated vary_43 (vary_43caz or vary_43cazonly)

df <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/Tables/s3_vary43_s3cavss3_anno_cazanno_08152025.csv")

ggplot(data = df, aes(x = S3CA_S3_log2FoldChange, y = -log10(S3CA_S3_padj))) +
  geom_point()

# Add threshold lines
ggplot(data = df, aes(x = S3CA_S3_log2FoldChange, y = -log10(S3CA_S3_padj))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point() 

# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

ggplot(data = df, aes(x = S3CA_S3_log2FoldChange, y = -log10(S3CA_S3_padj))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point() 

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$S3CA_S3_log2FoldChange > 1 & df$S3CA_S3_padj < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$S3CA_S3_log2FoldChange < -1 & df$S3CA_S3_padj < 0.05] <- "DOWN"

# Explore a bit
head(df[order(df$S3CA_S3_padj) & df$diffexpressed == 'DOWN', ])

ggplot(data = df, aes(x = S3CA_S3_log2FoldChange, y = -log10(S3CA_S3_padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  labs(color='Differential Expression') 


# Edit axis labels and limits
ggplot(data = df, aes(x = S3CA_S3_log2FoldChange, y = -log10(S3CA_S3_padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 50), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) # to customise the breaks in the x axis

# Note. with coord_cartesian() even if we have genes with p-values or log2FC ourside our limits, they will still be plotted.

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
dlup <- df[df$diffexpressed == "UP",]
df$delabel <- ifelse(df$GeneID %in% as.matrix(head(dlup[order(dlup$S3CA_S3_padj), "GeneID"], 10)), df$GeneID, NA)

dldown <- df[df$diffexpressed == "DOWN",]
df$delabel <- ifelse(df$GeneID %in% as.matrix(head(dldown[order(dldown$S3CA_S3_padj), "GeneID"], 10)), df$GeneID, df$delabel)

ggplot(data = df, aes(x = S3CA_S3_log2FoldChange, y = -log10(S3CA_S3_padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 40), xlim = c(-4, 7)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-4, 7, 2)) + # to customise the breaks in the x axis
  #geom_point() + geom_text(show.legend = FALSE, nudge_x = .2, nudge_y = 1)
  geom_label_repel(size = 4, show.legend = FALSE, min.segment.length = unit(0, 'lines')) +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.title = element_blank(), legend.position = "bottom")


## 11/25/2025 make heatmap of genes most of interest?
vary_43 <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/Tables/s3_vary43_s3cavss3_anno_cazanno_08152025.csv")
vary_43filt <- vary_43[abs(vary_43$S3CA_S3_log2FoldChange) > 3 & vary_43$S3CA_S3_padj < 0.01,]
# vary_43filtep <- vary_43filt[!grepl("expressed protein", vary_43filt$Protein_name),] # remove "expressed protein" genes waittttt but this removes some CAZymes so maybe ignore
# vary_43filtephp <- vary_43filtep[!grepl("hypothetical protein", vary_43filtep$Protein_name),] # remove "hypothetical protein" genes waittttt but this removes some CAZymes so maybe ignore
alp <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/Alphafold/Alphafold_genes_05082026.csv")
alp[] <- lapply(alp, function(x) {
  if (is.character(x)) gsub(",(?! )", ", ", x, perl = TRUE)
  else x
})

vary_43filta <- merge(vary_43filt, alp, all.x = TRUE)

pids <- c(259404, 987317, 997245, 1036142, 1048004, 1056260, 1289410, 1293984, 1298002, 1433957, 1653640, 1675787, 1725095, 1751806) # proteinIDs from alphafold that are CAZyme/cellulose binding related
vary_44 <- vary_43filta[is.na(vary_43filta$CAZ_anno) == FALSE | vary_43filta$ProteinID %in% pids,]

vary_44 <- vary_44 %>%
  mutate(Namealpha = ifelse(!is.na(Foldseekanno), 
                          paste(Protein_name, " (Foldseek: ", Foldseekanno, ")", sep = ""), 
                          Protein_name))

vary_44 <- vary_44 %>%
  mutate(Namefin = ifelse(!is.na(CAZ_anno), 
                          paste(Namealpha, CAZ_anno, sep = "; CAZyme family: "), 
                          Namealpha))
vary_44$Namefin <- gsub("\t", " ", vary_44$Namefin)
vary_44$Namefin <- gsub("\\n", "\n", vary_44$Namefin, fixed = TRUE)  # if stored as \n
vary_44$Namefin <- sub("^([a-z])", "\\U\\1", vary_44$Namefin, perl = TRUE)

# vary_43filt$Final_name <- NA
# for (i in 1:nrow(vary_43filt)) {
# if (is.na(vary_43filt$CAZ_anno[i]) == FALSE) {
#   vary_43filt$Final_name[i] <- paste(vary_43filt$Name[i],", ", "CAZyme type: ",vary_43filt$CAZ_anno[i], sep = "")
# } else {
#   vary_43filt$Final_name[i] <- vary_43filt$Name[i]
# }
# }
# 
# vary_43filtep <- vary_43filt[!(grepl("expressed protein", vary_43filt$Protein_name) == TRUE & is.na(vary_43filt$CAZ_anno) == TRUE),]
# vary_43filtephp <- vary_43filtep[!(grepl("hypothetical protein", vary_43filtep$Protein_name) == TRUE & is.na(vary_43filtep$CAZ_anno) == TRUE),]

# Make heatmap
newnames <- "CA + S3 vs.\nS3"
newnames <- "CA + S3 vs. S3"
# vary_44 <- as.data.frame(vary_43filtephp)

# ann_cols <- vary_44$Final_name
ann_cols <- vary_44$Namefin

library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "bigpts") - unit(3,"bigpts"), vjust = 0, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)


# make_bold_names <- function(mat, rc_fun, rc_names) {
#   bold_names <- rc_fun(mat)
#   ids <- rc_names %>% match(rc_fun(mat))
#   ids %>%
#     walk(
#       function(i)
#         bold_names[i] <<-
#         bquote(bold(.(rc_fun(mat)[i]))) %>%
#         as.expression()
#     )
#   bold_names
# }
# tobold <-vary_44$Final_name[grepl("CAZyme", vary_44$Final_name)]
# tobold

# row.names(vary_44) <- vary_44$Final_name
# s3_regulation_heatmap<-pheatmap(as.data.frame(vary_44$S3CA_S3_log2FoldChange), col=my_palette, cellwidth=50,cellheight=11, cluster_rows=TRUE, cluster_cols=F, 
                                # show_rownames=T, fontsize_row=8, fontsize_col=12, labels_col = as.expression(newnames), labels_row = make_bold_names(vary_44, rownames, tobold), border_color='white',breaks=seq(-5, 7, length.out=300), angle_col = 0)

s3_regulation_heatmap<-pheatmap(as.data.frame(vary_44$S3CA_S3_log2FoldChange), col=my_palette, cellwidth=50,cellheight=15, cluster_rows=TRUE, cluster_cols=F, 
                                show_rownames=T, fontsize_row=10, fontsize_col=12, labels_col = as.expression(newnames), labels_row = vary_44$Namefin, border_color='white',breaks=seq(-5, 7, length.out=300), angle_col = 0)


hm1 <- s3_regulation_heatmap

# save_pheatmap_pdf <- function(x, filename, width=8, height=11) {
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width=width, height=height)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }
# save_pheatmap_pdf(s3_regulation_heatmap, "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Tfigures/S3CAvsS3_vary_43filtephp_heatmap_12042025.pdf")


# 12/08/2025 making tiff version cus I think msystems prefers that?
s3_regulation_heatmap<-pheatmap(as.data.frame(vary_44$S3CA_S3_log2FoldChange), col=my_palette, cellwidth=50,cellheight=9, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=10, labels_col = as.expression(newnames), labels_row = make_bold_names(vary_44, rownames, tobold), border_color='white',breaks=seq(-5, 7, length.out=300), angle_col = 0)
