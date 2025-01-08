install.packages(c("edgeR", "limma", "Glimma", "gplots", "RColorBrewer"))
installed.packages("limma")
library(dplyr)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)

# Install matrixStats package
install.packages("matrixStats")

# Load the library
library(matrixStats)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("org.Mm.eg.db", "GO.db"))
BiocManager::install("limma")
BiocManager::install("Glimma")

# Load data
seqdata <- read.delim("GSE60450_Lactation-GenewiseCounts.txt", row.names = "EntrezGeneID")
seqdata <- seqdata[, -c(1)]  # Remove first two columns

col_names_seqdata <- colnames(seqdata)
col_names_seqdata= substr(col_names_seqdata, 1, 7)
colnames(seqdata) = col_names_seqdata

sampleinfo <- read.delim("SampleInfo.txt")

names <- c(sampleinfo[, 2])

seqdata <- seqdata[, names]

##############
#3#
############
myCPM = cpm(seqdata)
boxplot(myCPM, las=2, main="Distribution of Gene Expression (CPM)")

# Find entries > 0.5
above_threshold <- myCPM > 0.5

table(above_threshold)

genes_in_all_samples <- rowSums(myCPM > 0.5) == 12
table(genes_in_all_samples)

above_threshold_count <- rowSums(myCPM > 0.5)
counts.keep <- seqdata[above_threshold_count > 2, ]

dim(counts.keep)

DGElist_obj = DGEList(counts.keep)
str(DGElist_obj)

lib_sizes <- colSums(seqdata)
# Create a PDF file in the results folder
pdf("library_sizes_barplot.pdf", width=10, height=6)

# Create the barplot
barplot(lib_sizes, 
        main="Library Sizes", 
        ylab="Total Read Count", 
        xlab="Samples", 
        names.arg=colnames(seqdata),  # Use sample names on x-axis
        las=2,  # Rotate x-axis labels for better readability
        col="blue")

# Close the PDF device
dev.off()

y <- DGEList(counts=seqdata)

# Add log-transformed CPM data
y$logcpm <- cpm(y, log=TRUE)

boxplot(y$logcpm, 
        main="Log-Transformed Counts Across Samples", 
        xlab="Samples", 
        ylab="Log2 CPM", 
        las=2,  # Rotate x-axis labels
        col="lightblue")

# Add a line for median expression level
abline(h=median(y$logcpm), col="red", lty=2)


plotMDS(y$logcpm, 
        col = as.numeric(factor(sampleinfo$CellType)), 
        pch = 16)
# Add a legend
legend("topright", 
       legend = levels(factor(sampleinfo$CellType)), 
       col = 1:length(unique(sampleinfo$CellType)), 
       pch = 16)

plotMDS(y$logcpm, 
        col = as.numeric(factor(sampleinfo$Status)), 
        pch = 16)
# Add a legend
legend("topright", 
       legend = levels(factor(sampleinfo$Status)), 
       col = 1:length(unique(sampleinfo$Status)), 
       pch = 16)

sampleinfo_cor <- read.delim("SampleInfo_Corrected.txt")

plotMDS(y$logcpm, 
        col = as.numeric(factor(sampleinfo_cor$CellType)), 
        pch = 16)
# Add a legend
legend("topright", 
       legend = levels(factor(sampleinfo_cor$CellType)), 
       col = 1:length(unique(sampleinfo_cor$CellType)), 
       pch = 16)

plotMDS(y$logcpm, 
        col = as.numeric(factor(sampleinfo_cor$Status)), 
        pch = 16)
# Add a legend
legend("topright", 
       legend = levels(factor(sampleinfo_cor$Status)), 
       col = 1:length(unique(sampleinfo_cor$Status)), 
       pch = 16)

sel = order(apply(y$logcpm, 1, var), decreasing=TRUE)[1:500]

gene_variance <- rowVars(y$logcpm)
# Order genes by variance in descending order
var_order <- order(gene_variance, decreasing = TRUE)
# Select the top 500 most variable genes
top_500_genes <- var_order[1:500]
# Subset the original data with these top 500 genes
highly_variable_genes <- y$logcpm[top_500_genes, ]

heatmap(highly_variable_genes)

heatmap.2(highly_variable_genes)
# Select color palette
heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(256)
library(RColorBrewer)
pdf("heatmap_top500_variable_genes.pdf", width=10, height=12)
# Create heatmap with multiple customizations
heatmap.2(highly_variable_genes, 
          col = heatmap_colors,           # Color palette
          scale = "row",                  # Normalize by row (gene)
          trace = "none",                 # Remove trace line
          density.info = "none",          # Remove density info
          key = TRUE,                     # Show color key
          keysize = 1.5,                  # Adjust key size
          cexCol = 0.7,                   # Adjust column label size
          cexRow = 0.2,                   # Adjust row label size
          margins = c(10, 5),             # Adjust margins
          distfun = function(x) dist(x, method = "manhattan"),  # Manhattan distance
          main = "Top 500 Most Variable Genes")
dev.off()
