#!/usr/bin/env Rscript

# Clear workspace
rm(list = ls())

# Libraries
library(NBR)

# Input/output directories
data_dir <- paste0(getwd(),"/NBR/") # Where the "fc_dataset.csv" and "fc_matrices.rds" are stored
out_dir <- paste0(getwd(),"/NBR/") # Where the results are going to be located

# Phenotypic info
phen <- read.csv(file = file.path(data_dir,"fc_dataset.csv"),
                 stringsAsFactors = T)

# Filter by Group (i.e., remove controls)
idx <- which(phen$Group != "Ctrl")
phen <- phen[idx,]

# Relevel: setting low group as reference
phen$Group <- factor(phen$Group, levels = c("Alc","AlcStr","Str"))
# Intake default already 'Ctrl'
# Network matrices
net <- readRDS(file.path(data_dir,"fc_matrices.rds"))
net <- net[,,idx]
net_dim <- dim(net)

# Models
# Modelo 5. Intake-Age interaction.
mod <- "~Age*Group*Sex + Batch"
rdm <- "~1|RID"

# Number of cores
n_cores <- 60 # Something reasonable

# Number of permutations
n_perm <- 1000

# A priori threshold
ap_thr <- 0.05

# Run NBR
# Alternative greater than
tic <- Sys.time()
fit_nbr_mod5_gt_thr05 <- nbr_lme(net = net,
                   nnodes = net_dim[1],
                   idata = phen,
                   mod = mod,
                   rdm = rdm,
                   alternative = "greater",
                   nperm = n_perm,
                   thrP = ap_thr/2, #divided by two to account for the bi-sided
                   nudist = T,
                   cores = n_cores)
toc <- Sys.time()
print(toc-tic)

# Run NBR
# Alternative lower than
tic <- Sys.time()
fit_nbr_mod5_lt_thr05 <- nbr_lme(net = net,
                                 nnodes = net_dim[1],
                                 idata = phen,
                                 mod = mod,
                                 rdm = rdm,
                                 alternative = "lower",
                                 nperm = n_perm,
                                 thrP = ap_thr/2, #divided by two to account for the bi-sided
                                 nudist = T,
                                 cores = n_cores)
toc <- Sys.time()
print(toc-tic)

# It is recommendable to test more than one a priori threshold
# to test the consistency of the results.
# A priori threshold
ap_thr <- 0.01

# Run NBR
# Alternative greater than
tic <- Sys.time()
fit_nbr_mod5_gt_thr01 <- nbr_lme(net = net,
                                 nnodes = net_dim[1],
                                 idata = phen,
                                 mod = mod,
                                 rdm = rdm,
                                 alternative = "greater",
                                 nperm = n_perm,
                                 thrP = ap_thr/2, #divided by two to account for the bi-sided
                                 nudist = T,
                                 cores = n_cores)
toc <- Sys.time()
print(toc-tic)

# Run NBR
# Alternative lower than
tic <- Sys.time()
fit_nbr_mod5_lt_thr01 <- nbr_lme(net = net,
                                 nnodes = net_dim[1],
                                 idata = phen,
                                 mod = mod,
                                 rdm = rdm,
                                 alternative = "lower",
                                 nperm = n_perm,
                                 thrP = ap_thr/2, #divided by two to account for the bi-sided
                                 nudist = T,
                                 cores = n_cores)
toc <- Sys.time()
print(toc-tic)

# It is recommendable to test more than one a priori threshold
# to test the consistency of the results.
# A priori threshold
ap_thr <- 0.1

# Run NBR
# Alternative greater than
tic <- Sys.time()
fit_nbr_mod5_gt_thr10 <- nbr_lme(net = net,
                                 nnodes = net_dim[1],
                                 idata = phen,
                                 mod = mod,
                                 rdm = rdm,
                                 alternative = "greater",
                                 nperm = n_perm,
                                 thrP = ap_thr/2, #divided by two to account for the bi-sided
                                 nudist = T,
                                 cores = n_cores)
toc <- Sys.time()
print(toc-tic)

# Run NBR
# Alternative lower than
tic <- Sys.time()
fit_nbr_mod5_lt_thr10 <- nbr_lme(net = net,
                                 nnodes = net_dim[1],
                                 idata = phen,
                                 mod = mod,
                                 rdm = rdm,
                                 alternative = "lower",
                                 nperm = n_perm,
                                 thrP = ap_thr/2, #divided by two to account for the bi-sided
                                 nudist = T,
                                 cores = n_cores)
toc <- Sys.time()
print(toc-tic)

# Save results
save(fit_nbr_mod5_gt_thr10,
     fit_nbr_mod5_lt_thr10,
     fit_nbr_mod5_gt_thr05,
     fit_nbr_mod5_lt_thr05,
     fit_nbr_mod5_gt_thr01,
     fit_nbr_mod5_lt_thr01,
     file = file.path(out_dir,"fit_nbr_mod5_thrs10_05_01.Rdata"))
