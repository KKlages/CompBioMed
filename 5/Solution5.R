# Install and load required libraries
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(c("msa", "Biostrings", "ape", "phangorn"))

library(msa)
library(Biostrings)
library(ape)
library(phangorn)

# Read FASTA sequences
fasta_sequences <- readDNAStringSet("sars_spike.fasta")

# Multiple Sequence Alignment
alignment_muscle <- msaMuscle(fasta_sequences)

# Convert to phyDat object
phyDat_obj <- as.phyDat(alignment_muscle)

# Calculate distance matrix
dist_matrix <- dist.ml(phyDat_obj)

# Construct trees using different methods
tree_nj <- NJ(dist_matrix)      # Neighbor Joining
tree_upgma <- upgma(dist_matrix)  # UPGMA
tree_wpgma <- wpgma(dist_matrix)  # WPGMA

# Plot trees
par(mfrow=c(1,3))
plot(tree_nj, main="Neighbor Joining")
plot(tree_upgma, main="UPGMA")
plot(tree_wpgma, main="WPGMA")
###############################################
dist_matrix_JC69 <- dist.ml(phyDat_obj, model = "JC69")
dist_matrix_F81 <- dist.ml(phyDat_obj, model = "F81")

tree_nj_JC69 <- NJ(dist_matrix_JC69)      # Neighbor Joining
tree_upgma_JC69 <- upgma(dist_matrix_JC69)  # UPGMA
tree_wpgma_JC69 <- wpgma(dist_matrix_JC69)  # WPGMA

tree_nj_F81 <- NJ(dist_matrix_F81)      # Neighbor Joining
tree_upgma_F81 <- upgma(dist_matrix_F81)  # UPGMA
tree_wpgma_F81 <- wpgma(dist_matrix_F81)  # WPGMA

# Plot trees
par(mfrow=c(1,3))
plot(tree_nj_JC69, main="Neighbor Joining")
plot(tree_upgma_JC69, main="UPGMA")
plot(tree_wpgma_JC69, main="WPGMA")

# Plot trees
par(mfrow=c(1,3))
plot(tree_nj_F81, main="Neighbor Joining")
plot(tree_upgma_F81, main="UPGMA")
plot(tree_wpgma_F81, main="WPGMA")
#####
identical(dist_matrix_JC69, dist_matrix_F81)
#####

DNAbin_file <- as.DNAbin(fasta_sequences)

dist_matrix_K80 <- dist.dna(DNAbin_file, model = "K80")
dist_matrix_TN93 <- dist.dna(DNAbin_file, model = "TN93")
dist_matrix_GG95 <- dist.dna(DNAbin_file, model = "GG95")

tree_nj_K80 <- NJ(dist_matrix_K80)      # Neighbor Joining
tree_upgma_K80 <- upgma(dist_matrix_K80)  # UPGMA
tree_wpgma_K80 <- wpgma(dist_matrix_K80)  # WPGMA

tree_nj_TN93 <- NJ(dist_matrix_TN93)      # Neighbor Joining
tree_upgma_TN93 <- upgma(dist_matrix_TN93)  # UPGMA
tree_wpgma_TN93 <- wpgma(dist_matrix_TN93)  # WPGMA

tree_nj_GG95 <- NJ(dist_matrix_GG95)      # Neighbor Joining
tree_upgma_GG95 <- upgma(dist_matrix_GG95)  # UPGMA
tree_wpgma_GG95 <- wpgma(dist_matrix_GG95)  # WPGMA

# Plot trees
par(mfrow=c(1,3))
plot(tree_nj_K80, main="Neighbor Joining")
plot(tree_upgma_K80, main="UPGMA")
plot(tree_wpgma_K80, main="WPGMA")

# Plot trees
par(mfrow=c(1,3))
plot(tree_nj_TN93, main="Neighbor Joining")
plot(tree_upgma_TN93, main="UPGMA")
plot(tree_wpgma_TN93, main="WPGMA")

# Plot trees
par(mfrow=c(1,3))
plot(tree_nj_GG95, main="Neighbor Joining")
plot(tree_upgma_GG95, main="UPGMA")
plot(tree_wpgma_GG95, main="WPGMA")

#########################################
#2
#########################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# data.frame to be filled
wf_df <- data.frame()

# effective population sizes
sizes <- c(50, 100, 1000, 5000)

# starting allele frequencies
starting_p <- c(.01, .1, .5, .8)

# number of generations
n_gen <- 100

# number of replicates per simulation
n_reps <- 50

# run the simulations
for(N in sizes){
  for(p in starting_p){
    p0 <- p
    for(j in 1:n_gen){
      X <- rbinom(n_reps, 2*N, p)
      p <- X / (2*N)
      rows <- data.frame(replicate = 1:n_reps, N = rep(N, n_reps), 
                         gen = rep(j, n_reps), p0 = rep(p0, n_reps), 
                         p = p)
      wf_df <- bind_rows(wf_df, rows)
    }
  }
}

# plot it up!
p <- ggplot(wf_df, aes(x = gen, y = p, group = replicate)) +
  geom_path(alpha = .5) + facet_grid(N ~ p0) + guides(colour=TRUE)
p









