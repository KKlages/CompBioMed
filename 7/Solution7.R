# Install necessary packages and load libraries
install.packages(c("edgeR", "limma", "Glimma", "gplots", "RColorBrewer", "matrixStats"))
library(dplyr)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(matrixStats)

# Install Bioconductor packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("org.Mm.eg.db", "GO.db", "annotate"))
library(annotate)
# Load count matrix and sample info
y <- readRDS("data_mat.rds")
sampleinfo <- readRDS("sample_info.rds")

# Calculate logCPM
logCPM <- cpm(y, log = TRUE)

# Print dimensions to verify data load
cat("Dimensions of count matrix:", dim(y), "\n")
cat("Dimensions of logCPM matrix:", dim(logCPM), "\n")

# Normalize data
y <- calcNormFactors(y)

# Sort samples by normalization factors
sorted_samples <- y$samples[order(y$samples$norm.factors, decreasing = TRUE), ]

# Print top two samples with largest norm.factors
print(sorted_samples[1:2, ])

top_samples <- sorted_samples[1:2, ]
sample_names <- rownames(top_samples)

# Remove rows with NaN or infinite values
is_valid <- function(x) {
  !is.na(x) & is.finite(x)
}

# Ensure proper data frame for subsets
logCPM_valid <- logCPM[is_valid(logCPM[,1]) & is_valid(logCPM[,2]), ]

# Create manual MA plot for unadjusted log CPM
A <- rowMeans(logCPM_valid)
M <- logCPM_valid[,1] - logCPM_valid[,2]

plot(A, M, 
     main = "MA Plot of Log CPM (Unadjusted)",
     xlab = "A: Average log CPM",
     ylab = "M: Log Fold Change",
     pch = 20)
abline(h = 0, col = "red")

# Normalize logCPM for valid subset (if needed)
logCPM_norm_valid <- logCPM[is_valid(logCPM[,1]) & is_valid(logCPM[,2]), ]

# Create manual MA plot for normalized log CPM
A_norm <- rowMeans(logCPM_norm_valid)
M_norm <- logCPM_norm_valid[,1] - logCPM_norm_valid[,2]

plot(A_norm, M_norm, 
     main = "MA Plot of Normalized Log CPM",
     xlab = "A: Average log CPM",
     ylab = "M: Log Fold Change",
     pch = 20)
abline(h = 0, col = "red")

# Create combination variable
combinations <- paste(sampleinfo$CellType, sampleinfo$Status, sep=".")

# Create factor with correct levels
f <- factor(combinations, levels=unique(combinations))

# Create design matrix
design <- model.matrix(~ 0 + f)

# Rename columns to match your groups
colnames(design) <- levels(f)

# Verify dimensions match
cat("Columns in design matrix:", ncol(design), "\n")
cat("Columns in count data:", ncol(y), "\n")

# Proceed with voom normalization
v <- voom(y, design, plot = TRUE)


# Create boxplots side by side
par(mfrow=c(1,2))

# Boxplot of log-transformed CPMs
boxplot(logCPM, 
        main = "Log-transformed CPMs", 
        ylab = "Log Expression",
        las = 2)

# Boxplot of voom-normalized data
boxplot(v$E, 
        main = "Voom Normalized Data", 
        ylab = "Normalized Expression",
        las = 2)

# Reset plot parameters
par(mfrow=c(1,1))

fit <- lmFit(v, design)

# Create contrast matrix
# Specify the exact comparison you want
cont.matrix <- makeContrasts(
  BasalPregnantVsLactate = basal.pregnant - basal.lactate, 
  levels = design
)

# Fit the contrasts
fit.contrasts <- contrasts.fit(fit, cont.matrix)

# Apply empirical Bayes moderation to get p-values
fit.eb <- eBayes(fit.contrasts)

# Look at top differentially expressed genes
top_genes <- topTable(fit.eb, number = 10)
print(top_genes)

# Optional: Count significantly differentially expressed genes
summary(decideTests(fit.eb))

# Apply empirical Bayes moderation
fit.eb <- eBayes(fit.contrasts)

# Check results using topTable()
top_genes <- topTable(fit.eb)
print(top_genes)

# Use decideTests() to identify significantly regulated genes
de_results <- decideTests(fit.eb)

# Count significantly up and down regulated genes
sig_genes_summary <- summary(de_results)
print("Significantly regulated genes:")
print(sig_genes_summary)






