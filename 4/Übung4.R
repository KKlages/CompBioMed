install.packages('R.utils')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("maftools")

library(maftools)

laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
samplesummary =getSampleSummary(laml)
genesummary=getGeneSummary(laml)
clinicalsummary =getClinicalData(laml)

plotmafSummary(maf = laml)

oncoplot(
  maf=laml,
  top = 5
)

laml.titv = titv(maf = laml)
plotTiTv(laml.titv)

lollipopPlot(maf = laml, gene = 'FLT3')

######################################################

install.packages("survival")
install.packages("ggfortify")
library(survival)
library(ggfortify)


# Create survival object for all patients
surv_obj <- Surv(time = clinicalsummary$days_to_last_followup, 
                 event = clinicalsummary$Overall_Survival_Status)

# Fit survival curve
fit <- survfit(surv_obj ~ 1)  # The "~ 1" means we're not grouping by any variable

# Create the plot
plot(fit, 
     xlab = "Time (days)", 
     ylab = "Overall Survival Probability",
     main = "Overall Survival in AML Patients",
     conf.int = FALSE)



surv_fit_fab <- survfit(Surv(days_to_last_followup, Overall_Survival_Status) ~ FAB_classification, 
                        data = clinicalsummary)


# Create plot with autoplot
autoplot(surv_fit_fab) +
  labs(x = "Time (days)",
       y = "Overall Survival Probability",
       title = "Overall Survival by FAB Classification",
       color = "FAB Class") +
  theme_minimal()

# Remove rows with -Inf in any of the columns
clean_clinical <- clinicalsummary[(clinicalsummary$days_to_last_followup) != -Inf ]

# Verify that there are no more -Inf values in the data
sapply(clean_clinical, function(x) sum(is.infinite(x)))  # Should return 0 for all columns


test1 <- list(time=clean_clinical$days_to_last_followup, 
              status=clean_clinical$Overall_Survival_Status,
              x = clean_clinical$FAB_classification)

cox = coxph(Surv(time, status) ~ x , test1) 
cox_summary <- summary(cox)

a = mafSurvival(maf = laml, genes = 'RUNX1', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)
############################################################################
# Survival Analysis on Gene Mutations
############################################################################

# Top 20 geneset survival analysis
prog_geneset <- survGroup(
  maf = laml,
  top = 20,
  geneSetSize = 1,
  time = "days_to_last_followup",
  Status = "Overall_Survival_Status",
  verbose = FALSE
)

print(prog_geneset)

# Extract p-values and visualize distribution
P_values <- prog_geneset$P_value
#dev.off()  # Close any open graphics devices
#windows()  # Open a new Windows graphics device

# Histogram of P-values
hist(
  P_values, 
  breaks = 10,
  main = "Distribution of P-values for Gene Mutation Survival Analysis",
  xlab = "P-value",
  ylab = "Frequency",
  col = "lightblue",
  border = "black"
)
abline(v = 0.05, col = "red", lty = 2, lwd = 2)

# Add legend for the threshold line
legend(
  "topright", 
  legend = "p=0.05 threshold",
  col = "red", 
  lty = 2, 
  lwd = 2
)

# Analyze survival for a specific gene
a <- mafSurvival(
  maf = laml, 
  genes = 'TP53', 
  time = 'days_to_last_followup', 
  Status = 'Overall_Survival_Status', 
  isTCGA = TRUE
)

# Perform p-value adjustments
bonf_adj <- p.adjust(P_values, method = "bonferroni")
bh_adj <- p.adjust(P_values, method = "BH")

# Create a comparison table
results_table <- data.frame(
  Gene = prog_geneset$Gene_combination,
  Original_P = P_values,
  Bonferroni_adj = bonf_adj,
  BH_adj = bh_adj
)

# Sort the table by original p-value
results_table <- results_table[order(results_table$Original_P),]
#Bonferroni more conservative
############################################################################
# SARS-CoV Sequence Retrieval and ORF Identification
############################################################################

# Install and load required packages
if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
library(rentrez)

# Search and fetch SARS-CoV sequence
sars_search <- entrez_search(db = "nucleotide", term = "AY274119.3")
sars_seq <- entrez_fetch(db = "nucleotide", id = sars_search$ids, rettype = "fasta", retmode = "text")

# Extract sequence
seq_lines <- strsplit(sars_seq, "\n")[[1]]
seq_lines <- seq_lines[!grepl("^>", seq_lines)]
seq_string <- paste(seq_lines, collapse = "")

# Identify ORFs
if (!requireNamespace("ORFik", quietly = TRUE)) BiocManager::install("ORFik")
library(ORFik)

findORFs(seq_string, longestORF = TRUE, minimumLength = 100)

# Extract ORF sequence
ORF_seq <- substr(seq_string, start = 21492, stop = 25259)

# Convert to amino acid sequence
dna <- DNAString(x = ORF_seq)
cov_aa <- translate(dna)

############################################################################
# Multiple Sequence Alignment of SARS-CoV and Related Viruses
############################################################################
viruses <- readRDS("F:/CoBi/BioMed/Tutorial/4/viruses.rds")
human_cov_seq <- viruses[[3]] 
human_cov <- AAString(human_cov_seq)
length(cov_aa)
length(human_cov)

install.packages("pwalign")
library(pwalign)

alignment <- pwalign::pairwiseAlignment(pattern = cov_aa, subject = human_cov, substitutionMatrix = "BLOSUM62")
alignment
#############################################################################
BiocManager::install(c("Biostrings", "muscle", "msa"))
library(Biostrings)
library(muscle)
library(msa)
viruses_seqs <- AAStringSet(c(lapply(viruses, AAString), cov_aa))
alignment_seqs<- msa(viruses_seqs)
msaPrettyPrint(alignment_seqs)
#############
install.packages("phangorn")
library(phangorn)
phyDat_obj <- as.phyDat(alignment_seqs)
dis_mat <- dist.ml(phyDat_obj)
nj_tree<- NJ(dis_mat)
plot(nj_tree)
###
#Palm civet closest related species