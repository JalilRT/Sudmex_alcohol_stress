#!/usr/bin/env Rscript

# Clear workspace
rm(list = ls())

# Libraries
library(NBR)

# Input/output directories
data_dir <- # Where the "fc_dataset.csv" and "fc_matrices.rds" are stored
out_dir <- # Where the results are going to be located
  
# Phenotypic info
phen <- read.csv(file = file.path(data_dir,"fc_dataset.csv"),
                   stringsAsFactors = T)
# Relevel: setting control group as reference
# Intake default already 'Ctrl'
# Network matrices
net <- readRDS(file.path(data_dir,"fc_matrices.rds"))
net_dim <- dim(net)

# Models
# Modelo 2. Intake-Age interaction.
mod <- "~Age*Intake + Sex  + Batch"
rdm <- "~1|RID"

# Number of cores
n_cores <- # Something reasonable
  
# Number of permutations
n_perm <- 1000

# A priori threshold
ap_thr <- 0.05

# Run NBR
tic <- Sys.time()
fit_nbr_mod2_aov_thr05 <- nbr_lme_aov(net = net,
                                 nnodes = net_dim[1],
                                 idata = phen,
                                 mod = mod,
                                 rdm = rdm,
                                 nperm = n_perm,
                                 thrP = ap_thr,
                                 nudist = T,
                                 cores = n_cores)
toc <- Sys.time()
print(toc-tic)

# It is recommendable to test more than one a priori threshold
# to test the consistency of the results.
# A priori threshold
ap_thr <- 0.01

# Run NBR
tic <- Sys.time()
fit_nbr_mod2_aov_thr01 <- nbr_lme_aov(net = net,
                                 nnodes = net_dim[1],
                                 idata = phen,
                                 mod = mod,
                                 rdm = rdm,
                                 nperm = n_perm,
                                 thrP = ap_thr,
                                 nudist = T,
                                 cores = n_cores)
toc <- Sys.time()
print(toc-tic)

# A priori threshold
ap_thr <- 0.1

# Run NBR
tic <- Sys.time()
fit_nbr_mod2_aov_thr10 <- nbr_lme_aov(net = net,
                                 nnodes = net_dim[1],
                                 idata = phen,
                                 mod = mod,
                                 rdm = rdm,
                                 nperm = n_perm,
                                 thrP = ap_thr,
                                 nudist = T,
                                 cores = n_cores)
toc <- Sys.time()
print(toc-tic)

# Save results
save(fit_nbr_mod2_aov_thr10,
     fit_nbr_mod2_aov_thr05,
     fit_nbr_mod2_aov_thr01,
     file = file.path(out_dir,"fit_nbr_mod2_aov_thrs10_05_01.Rdata"))