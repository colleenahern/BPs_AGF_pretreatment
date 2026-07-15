# Colleen Ahern
# 08/12/2025
# Complete DESeq analysis on S3CA count matrix from featureCounts
# following tutorial: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library("tximport")
library("readr")
library("tximportData")
library(tidyverse)
library(pasilla)

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
colnames(s3countdata) %in% s3sampleinfo$libraryName
table(colnames(s3countdata)==s3sampleinfo$libraryName)
# s3countdata <- s3countdata[,as.character(s3sampleinfo$libraryName)]
# table(colnames(s3countdata)==s3sampleinfo$libraryName) # Now they match
colnames(s3countdata) <- s3sampleinfo$sampleName

# my S3CAAS samples suck so remove those, also remove S3_5 cus it's an outlier
s3countdata <- s3countdata[,c(1:4,6:10)]
s3sampleinfo <- s3sampleinfo[c(1:4,6:10),]

rownames(s3sampleinfo) <- s3sampleinfo$sampleName
s3coldata <- s3sampleinfo
s3cts <- s3countdata

all(rownames(s3coldata) == colnames(s3cts))
## TPM cutoff, not sure if I want to use this rn, just follow the luminal mice tutorial for now (07/24/2025)
# install.packages("remotes")
# remotes::install_github("davidrequena/drfun")
# library(DRnaSeq)
# all(rownames(s3countdata) == s3seqdata$Geneid)
# gl <- s3seqdata$Length
# s3TPM <- as.data.frame(tpm(s3countdata, gl))

# With the count matrix and the sample info, construct a DESeqDataSet:
library("DESeq2")
s3dds <- DESeqDataSetFromMatrix(countData = s3cts,
                                colData = s3coldata,
                                design = ~ groupName)
s3dds

# Pre-filter

keep <- rowSums(counts(s3dds)) > 0
s3dds <- s3dds[keep,]   # reduced number of genes from 27677 to 25498

# Set the factor levels 
s3dds$groupName <- factor(s3dds$groupName, levels = c("S3","S3CA"))

# Differential expression analysis: compare each substrate to glucose

s3dds <- DESeq(s3dds)
s3res <- results(s3dds)
s3res
s3res <- results(s3dds, contrast=c("groupName","S3CA","S3"))
s3res

# S3CA vs S3
s3res_S3CA_S3 <- results(s3dds, contrast=c("groupName","S3CA","S3"))
s3res_S3CA_S3
s3resmat_S3CA_S3 <- as.matrix(s3res_S3CA_S3)
head(s3resmat_S3CA_S3)
colnames(s3resmat_S3CA_S3) <- paste("S3CA_S3", colnames(s3resmat_S3CA_S3), sep = "_")
head(s3resmat_S3CA_S3)
s3resmat_S3CA_S3 <- as.data.frame(s3resmat_S3CA_S3)
head(s3resmat_S3CA_S3)
sum(s3resmat_S3CA_S3$S3CA_S3_padj < 0.05, na.rm=TRUE)

s3_DESeq <- s3resmat_S3CA_S3
# s3_DESeq <- s3_DESeq[, c(2,6)]
s3_DESeq <- as.data.frame(s3_DESeq)
s3_DESeq$S3CA_S3_Sig <- ifelse(s3_DESeq$S3CA_S3_padj < 0.05,"TRUE","FALSE")

s3_DESeq$GeneID <- rownames(s3_DESeq)
# s3_DESeq <- s3_DESeq[,c(4,1:3)]

library(readr)
write_xlsx(s3_DESeq, "/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/s3_DESeq_08132025.xlsx")
write_csv(s3_DESeq, "/Users/colleenahern/Documents/Magda_BPs_experiment/RNAseq_S3CA/DESeq2_analysis/s3_DESeq_ALL_08222025.csv")

