# Colleen Ahern
# 08/11/2025
# MDS plot from S3CA featureCounts count matrix

library(edgeR)
library(limma)

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
all(colnames(s3countdata)==s3sampleinfo$libraryName)
stopifnot(all(colnames(s3countdata) == s3sampleinfo$libraryName))
colnames(s3countdata) <- s3sampleinfo$sampleName

# Remove S3CAAS samples
keep_samples <- s3sampleinfo$sampleName[!grepl("S3CAAS", s3sampleinfo$sampleName)]
s3countdata <- s3countdata[, keep_samples]
s3sampleinfo <- s3sampleinfo[s3sampleinfo$sampleName %in% keep_samples, ]

s3y <- DGEList(s3countdata)
s3y
names(s3y)
s3y$samples

group <- paste(s3sampleinfo$groupName)
group

# Convert to factor
group <- factor(group)
group

# Add the group information into the DGEList
s3y$samples$group <- group
s3y$samples

# The tutorial uses CPM (only normalizes for library size) as a cutoff
s3CPM <- cpm(s3countdata)
head(s3CPM)

# Which values in myCPM are greater than 0.5?
thresh <- s3CPM > 0.15
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)
# Summary of how many TRUEs there are in each row
table(rowSums(thresh))

# We would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
summary(keep)

# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(s3CPM[,1],s3countdata[,1])

# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(s3CPM[,1],s3countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.15) # 0.15 CPM corresponds to 10-15 counts, not 0.5, so use that

# Now that we’ve checked our filtering method we will filter the DGEList object.
s3y <- s3y[keep, keep.lib.sizes=FALSE]

#First, we can check how many reads we have for each sample in the y.
s3y$samples$lib.size

# We can also plot the library sizes as a barplot to see whether there are any major discrepancies between the samples more easily.
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(s3y$samples$lib.size,names=colnames(s3y),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# we can also adjust the labelling if we want
barplot(s3y$samples$lib.size/1e06, names=colnames(s3y), las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")

# Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts. 
# Next we’ll use box plots to check the distribution of the read counts on the log2 scale. We can use the CPM function to get log2 
# counts per million, which are corrected for the different library sizes. The CPM function also adds a small offset to avoid taking 
# log of zero.

# Get log2 counts per million
logcounts <- cpm(s3y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# Apply normalisation to DGEList object
s3y <- calcNormFactors(s3y)
s3y$samples

# Get log2 counts per million
logcounts <- cpm(s3y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (normalised)")

# plot the MDS
stopifnot(identical(levels(factor(s3sampleinfo$groupName)), c("S3", "S3CA")))
col.cell <- c("blue", "deeppink")[factor(s3sampleinfo$groupName)]
data.frame(s3sampleinfo$groupName,col.cell)
par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
plotMDS(s3y)
plotMDS(s3y,col=col.cell, pch=16)
legend("right", inset = c(-.3, 0), fill=c("blue", "deeppink"), legend=c("S3", "S3 + CA"), title = expression(bold("Substrate")), bty = "n") 


# Remove S3_5 from DESeq2 analysis because it is an outlier

