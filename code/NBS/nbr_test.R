#!/usr/bin/env Rscript

# Clear workspace
rm(list = ls())
# Save original parameters
op <- par()
# Set working directory
setwd("Desktop/Zeus/Sudomex-alcohol-rat")

# Libraries
library(nlme)
library(NBR)

# Input/output directories
data_dir <- file.path(getwd(),"Data")
out_dir <- file.path(getwd(),"Results")

# Phenotypic info
phen <- read.csv(file = file.path(data_dir,"fc_dataset.csv"),
                 stringsAsFactors = T)
phen$Group <- factor(phen$Group, levels = c("Ctrl","Alc"))
# Network matrices
net <- readRDS(file.path(data_dir,"fc_matrices.rds"))
net_dim <- dim(net)

# Models
# Modelo 1. Group-Age interaction.
mod <- "~Age*Group + Sex + Batch"
rdm <- "~1|RID"

# Estimate time in one edge
phen$ytest <- net[1,2,]
tic <- Sys.time()
fit_test <- lme(fixed = as.formula(paste("ytest",mod)),
                data = phen,
                random = as.formula(rdm))
toc <- Sys.time()
# One edge
print(toc-tic)
# Whole-network
edge_n <- (net_dim[1]*(net_dim[2]-1)*0.5)
print((toc-tic)*edge_n)

# Number of permutations
n_perm <- 1

# A priori threshold
ap_thr <- 0.05

# Run NBR
# Alternative greater than
tic <- Sys.time()
fit_nbr_mod1_gt_thr05 <- nbr_lme(net = net,
                   nnodes = net_dim[1],
                   idata = phen,
                   mod = mod,
                   rdm = rdm,
                   alternative = "greater",
                   nperm = n_perm,
                   thrP = ap_thr/2, #divided by two to account for the bi-sided
                   nudist = T,
                   cores = 4)
toc <- Sys.time()
print(toc-tic)

# Run NBR
# Alternative lower than
tic <- Sys.time()
fit_nbr_mod1_lt_thr05 <- nbr_lme(net = net,
                                 nnodes = net_dim[1],
                                 idata = phen,
                                 mod = mod,
                                 rdm = rdm,
                                 alternative = "lower",
                                 nperm = n_perm,
                                 thrP = ap_thr/2, #divided by two to account for the bi-sided
                                 nudist = T,
                                 cores = 4)
toc <- Sys.time()
print(toc-tic)

# It is recommendable to test more than one a priori threshold
# to test the consistency of the results.
# A priori threshold
ap_thr <- 0.01

# Run NBR
# Alternative greater than
tic <- Sys.time()
fit_nbr_mod1_gt_thr01 <- nbr_lme(net = net,
                                 nnodes = net_dim[1],
                                 idata = phen,
                                 mod = mod,
                                 rdm = rdm,
                                 alternative = "greater",
                                 nperm = n_perm,
                                 thrP = ap_thr/2, #divided by two to account for the bi-sided
                                 nudist = T,
                                 cores = 4)
toc <- Sys.time()
print(toc-tic)

# Run NBR
# Alternative lower than
tic <- Sys.time()
fit_nbr_mod1_lt_thr01 <- nbr_lme(net = net,
                                 nnodes = net_dim[1],
                                 idata = phen,
                                 mod = mod,
                                 rdm = rdm,
                                 alternative = "lower",
                                 nperm = n_perm,
                                 thrP = ap_thr/2, #divided by two to account for the bi-sided
                                 nudist = T,
                                 cores = 4)
toc <- Sys.time()
print(toc-tic)

# Save results
save(fit_nbr_mod1_gt_thr05,
     fit_nbr_mod1_lt_thr05,
     fit_nbr_mod1_gt_thr01,
     fit_nbr_mod1_lt_thr01,
     file = file.path(out_dir,"fit_nbr_mod1_thrs05_01.Rdata"))
