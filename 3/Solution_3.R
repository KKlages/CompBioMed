install.packages('R.utils')
options(device = "windows")  # For Windows

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

a = mafSurvival(maf = laml, genes = 'FLT3', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)
