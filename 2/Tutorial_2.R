
library(snpStats)
library(hexbin)
library(viridis)  
library(RColorBrewer)
library(lattice)
load("snps_10.dat")
show(snps.10)


sample.qc <- row.summary(snps.10)
use <- sample.qc$Heterozygosity > 0
snps.10 <- snps.10[use, ]
subject.support <- subject.support[use, ]
if.case <- subject.support$cc == 1
if.control <- subject.support$cc == 0


tests <- single.snp.tests(cc, data = subject.support, snp.data = snps.10)
summary(tests)


snpsum <- col.summary(snps.10)
use <- snpsum$MAF > 0.01 & snpsum$z.HWE^2 < 200
tests <- tests[use]
position <- snp.support[use, "position"]
p1 <- p.value(tests, df = 1)


hb <- hexbin(position, -log10(p1), xbin = 50)



hexbinplot(-log10(p1) ~ position,
           main = "Manhattan Plot",
           xlab = "Chromosome Position",
           ylab = "-log10(p-value)",
           colramp = viridis::viridis,
           panel = function(x, y, ...) {
             panel.hexbinplot(x, y, ...)
             panel.abline(h = 2.995732, col = "red", lty = 2, lwd = 2)
           }
)



snps_subset <- snps.10[, 1:1000]
ld_subset <- ld(snps_subset, depth =10 ,stats = c("D.prime", "R.squared"))


R_squared_subset <- as.matrix(ld_subset$R.squared)

levelplot(R_squared_subset, col.regions = viridis::viridis(100),
          xlab = "SNP Index (First 1000)", ylab = "SNP Index (First 1000)",
          main = expression(R^2 ~ " LD Heatmap (First 1000 SNPs)"))


options(repr.plot.width = 25, repr.plot.height = 8)  # Adjust dimensions as needed

manhattan_plot <- function(stratum_name) {
  indices <- subject.support$stratum == stratum_name
  snp_data <- snps.10[indices, ]
  subject_data <- subject.support[indices, ]
  
  tests <- single.snp.tests(cc, data = subject_data, snp.data = snp_data)
  p_values <- p.value(tests, df = 1)
  position <- snp.support$position  
  

  hexbinplot(-log10(p_values) ~ position,
             main = paste("Manhattan Plot -", stratum_name),
             xlab = "Chromosome Position",
             ylab = "-log10(p-value)",
             colramp = viridis::viridis,
             panel = function(x, y, ...) {
               panel.hexbinplot(x, y, ...)
               panel.abline(h = 2.995732, col = "red", lty = 2, lwd = 2)
             }
  )
}


manhattan_plot("CEU")
manhattan_plot("JPT+CHB")




indices <- subject.support$stratum == "CEU"  
snp_data <- snps.10[indices, ]



snps_subset <- snp_data[, 1:1000]
ld_subset <- ld(snps_subset, depth = 10, stats = c("D.prime", "R.squared"))
R_squared_subset <- as.matrix(ld_subset$R.squared)



levelplot(R_squared_subset, col.regions = viridis::viridis(100),
          xlab = "SNP Index (First 1000)", ylab = "SNP Index (First 1000)",
          main = expression(R^2 ~ " LD Heatmap CEU(First 1000 SNPs)"),
          colorkey = TRUE)


indices2 <- subject.support$stratum == "JPT+CHB"  
snp_data2 <- snps.10[indices2, ]



snps_subset2 <- snp_data2[, 1:1000]
ld_subset2 <- ld(snps_subset2, depth = 10, stats = c("D.prime", "R.squared"))
R_squared_subset2 <- as.matrix(ld_subset2$R.squared)



levelplot(R_squared_subset2, col.regions = viridis::viridis(100),
          xlab = "SNP Index (First 1000)", ylab = "SNP Index (First 1000)",
          main = expression(R^2 ~ " LD Heatmap JPT+CHB(First 1000 SNPs)"),
          colorkey = TRUE)

