# This file contains the code that was used for simulating the four different
# synthetic scRNA-seq scenarios and for batch correction of Scenarios 1 and 4.
# Specifically, it was used for producing the plots in Figures 3.1-3.6 of my thesis.

library(scater)   # SingleCellExperiment class and accompanying functions
library(splatter)  # for scRNA-seq data simulation
library(scran)  # for HVG selection
library(batchelor)  # for fastMNN batch correction
library(harmony)  # for Harmony batch correction
library(Seurat)  # for Seurat batch correction
library(kBET)  # for kBET batch mixing metric
library(lisi)  # for LISI batch mixing metric
library(FNN)  # fast NN search for local KL divergence metric
library(kSamples)  # for the Anderson-Darling test 
library(tidyverse)  # for dataframe manipulations and ggplot
library(egg)  # for combining ggplots
source("../helpers.R")  # contains some handmade functions

# Specify output directory for plots
output_dir <- "../plots_report/sim/"

#### Simulation scenarios ####
## Scenario 1: 1 group, 2 batches
## Varying batch.facLoc and batch.facScale parameters (0.1, 0.07, and 0.04)
# param 0.1
params_0.1 <- newSplatParams(
  nGenes = 15000,
  batchCells = c(300, 300)
  # default batch.facLoc = 0.1 and batch.facScale = 0.1
)
sim_0.1 <- 
  sim_prep(params_0.1) %>%
  compute_svd()

# param 0.07
params_0.07 <- newSplatParams(
  nGenes = 15000,
  batchCells = c(300, 300),
  batch.facLoc = 0.07,
  batch.facScale = 0.07
)
sim_0.07 <- 
  sim_prep(params_0.07) %>%
  compute_svd()

# param 0.04
params_0.04 <- newSplatParams(
  nGenes = 15000,
  batchCells = c(300, 300),
  batch.facLoc = 0.04,
  batch.facScale = 0.04,
  seed = 1
)
sim_0.04 <- 
  sim_prep(params_0.04) %>%
  compute_svd()


# Plots for Figure 3.1
sce_list <- list(sim_0.1 = sim_0.1, sim_0.07 = sim_0.07, sim_0.04 = sim_0.04)
singvals_params <- find_plot_params(sce_list, sv_test = FALSE)
singvect_params <- find_plot_params(sce_list, plot_type = "singvect", vectors = 1:2)

x <- singvals_params$xlims[2]
y <- singvals_params$ylims[2]
for (k in c(0.1, 0.07, 0.04)) {
  p_svals <-
    plot_hist(sce_list[[paste0("sim_", k)]], pparams = singvals_params) +
    annotate(
      "text", x = x, y = y,
      label = paste("list(mu,beta)", "==", k), parse = TRUE,
      hjust = "inward", vjust = "inward"
    )
  filepath <- paste0("../plots_report/sim/singvals_1g2b_", k, ".pdf")
  
  p_svect <- 
    plot_vector(sce_list[[paste0("sim_", k)]], vectors = 1:2, pparams = singvect_params) +
    scale_y_continuous(breaks = c(-0.1, 0, 0.1))

  p <- ggarrange(p_svals, p_svect, ncol = 1)
  filepath <- paste0("../plots_report/sim/1g2b_", k, ".pdf")
  ggsave(filepath, p, width = 60, height = 120, units = "mm", device = "pdf")
}


## Scenario 2: 1 group, multiple batches #
# 1 group, 3 batches
params_3b <- newSplatParams(
  nGenes = 15000,
  batchCells = rep(200, 3),
)
sim_3b <- 
  sim_prep(params_3b) %>% # includes simple library size normalization
  compute_svd()

# 1 group, 4 batches
params_4b <- newSplatParams(
  nGenes = 15000,
  batchCells = rep(150, 4),
)
sim_4b <- 
  sim_prep(params_4b) %>%
  compute_svd()

# 1 group, 5 batches
params_5b <- newSplatParams(
  nGenes = 15000,
  batchCells = rep(120, 5),
#  batch.facLoc = 0.2,
#  batch.facScale = 0.2,
#  seed = 1  # to get more distinct sv's
)
sim_5b <- 
  sim_prep(params_5b) %>%
  compute_svd()

# Plots for Figure 3.2
sce_list <- list(sim_3b = sim_3b, sim_4b = sim_4b, sim_5b = sim_5b)
singvals_params <- find_plot_params(sce_list, sv_test = FALSE)
singvect_params <- find_plot_params(sce_list, plot_type = "singvect", vectors = c(3,1))

for (k in 3:5) {
  p_singvals <-
    plot_hist(sce_list[[paste0("sim_", k, "b")]], pparams = singvals_params) +
    scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9, 1.2))

  p_singvect <- 
    plot_vector(sce_list[[paste0("sim_", k, "b")]], vectors = c(k-1, 1), pparams = singvect_params) +
    scale_y_continuous(breaks = c(-0.07, 0, 0.07))
  
  p <- ggarrange(p_singvals, p_singvect, ncol = 1)
  filepath <- paste0("../plots_report/sim/1g", k, "b.pdf")
  ggsave(filepath, p, width = 60, height = 120, units = "mm", device = "pdf")
}


## Scenario 3: 2 cell types, 2 batches
sim_2g2b <- list()
batch_effect <- c(0.12, 0.15, 0.18)
for (k in batch_effect) {
  params <- newSplatParams(
    nGenes = 15000,
    batchCells = c(300, 300),
    group.prob = c(0.5, 0.5),
    batch.facLoc = k,
    batch.facScale = k
  )
  sim_2g2b[[paste0("sim_", k)]] <- 
    sim_prep(params) %>%
    compute_svd()
}

# Plots for Figure 3.3
singvals_params <- find_plot_params(sim_2g2b, sv_test = FALSE)
singvect_params <- find_plot_params(sim_2g2b, plot_type = "singvect", vectors = c(1,2))

for (k in batch_effect) {
  p_singvals <-
    plot_hist(sim_2g2b[[paste0("sim_", k)]], pparams = singvals_params)
  
  p_singvect <- 
    plot_vector(sim_2g2b[[paste0("sim_", k)]], vectors = c(1, 2), pparams = singvect_params) +
    scale_y_continuous(breaks = c(-0.05, 0, 0.05))
  
  p <- ggarrange(p_singvals, p_singvect, ncol = 1)
  filepath <- paste0("../plots_report/sim/2g2b_", k, ".pdf")
  ggsave(filepath, p, width = 60, height = 120, units = "mm", device = "pdf")
}


## Scenario 4: Unbalanced data
params_sc4 <- newSplatParams(
  nGenes = 15000,
  batchCells = c(400, 400),
  group.prob = c(0.4, 0.2, 0.2, 0.2),
  batch.facLoc = 0.1,
  batch.facScale = 0.1
)
sc4 <- 
  splatSimulate(params_sc4, method = "groups", verbose = FALSE) %>%
  minimiseSCE(colData.keep = c("Batch", "Group"), verbose = FALSE) %$%
  
  # Remove all cells in Batch 2 from Group 4
  .[, !(.$Batch == "Batch2" & .$Group == "Group4")] %$%
  
  # Remove undetected genes
  .[nexprs(., byrow = TRUE, detection_limit = 1) >= 2, ] %>%
  
  logNormCounts() %>%
  filter_HVGs(percent_topHVGs = 0.10) %>%
  compute_svd()

reducedDim(sc4, "PCA_logcounts") <- calculatePCA(sc4, ntop = nrow(sc4))
plotReducedDim(sc4, "PCA_logcounts", colour_by = "Batch")

# Plots for Figure 3.4
plot_hist(sc4, assayname = "logcounts") +
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9))
ggsave(
  paste0(output_dir, "unbalanced_svhist.pdf"), 
  width = 60, height = 60, units = "mm", device = "pdf"
)
plot_vector(sc4, vectors = 1:2) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05))
ggsave(
  paste0(output_dir, "unbalanced_vector12.pdf"), 
  width = 60, height = 60, units = "mm", device = "pdf"
)
plot_vector(sc4, vectors = 3:4)
ggsave(
  paste0(output_dir, "unbalanced_vector34.pdf"), 
  width = 60, height = 60, units = "mm", device = "pdf"
)


#### Batch correction ####

# Figure 3.5
singvals_params <- find_plot_params(list(sim_0.1))
singvals_params$ylims <-
  sim_0.1 %>%
  harmony_wrapper(d = 2) %>%
  compute_svd("harmony") %>%
  plot_hist("harmony", sv_test = TRUE) %>%
  get_lims("y")

# Different values for d for fastMNN and Harmony and Seurat
set.seed(17)
D <- c(2, 20, seq(50, 600, by = 50))
batch_singvals <- list(
  mnn = numeric(length(D)), 
  harmony = numeric(length(D)),
  seurat = numeric(length(D))
)

for (i in 1:length(D)) {
  sim_0.1 <-
    sim_0.1 %>%
    harmony_wrapper(d = D[i]) %>%
    fastmnn_wrapper2(d = D[i]) %>%
    compute_svd("harmony") %>%
    compute_svd("mnn")
  
  # Seurat doesn't work if d >= #cells in the smallest batch
  if (D[i] < 300) {
    sim_0.1 <-
      seurat_wrapper(sim_0.1, d = D[i]) %>%
      compute_svd("seurat")
  }
  
  # Compute max singval associated with batch effect
  P <- nrow(sim_0.1)
  for (assayname in c("harmony", "mnn", "seurat")) {
    SVD <- metadata(sim_0.1)[[assayname]]$SVD
    idx <- identify_batch_vectors(
      vectors = SVD$v[, -ncol(SVD$v)], batch = sim_0.1$Batch, 
      adjust = TRUE, idx = TRUE, threshold = 0.05
    )
    
    batch_singvals[[assayname]][i] <- SVD$d[idx[1]] / sqrt(P - 1)  # record highest batch effect singular value
  }

  caption_dict <- list("Uncorrected", "Harmony", "fastMNN", "Seurat")
  names(caption_dict) <- names(metadata(sim_0.1))
  
  # Plot singular value histograms (Figure 3.5B)
  annot_x <- max(singvals_params$breaks)
  annot_y <- max(singvals_params$ylims)
  for (assayname in names(metadata(sim_0.1))) {
    if (assayname != "seurat" | D[i] < 300) {  # make sure that if Seurat, only plot if d < 300
      p_singvals <- 
        plot_hist(sim_0.1, assayname, sv_test = TRUE, pparams = singvals_params) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        annotate(
          "text", x = annot_x, y = annot_y,
          label = paste0(caption_dict[[assayname]], "\n", "d = ", D[i]),
          hjust = "inward", vjust = "inward", size = 3
        )
      filepath <- paste0("../plots_report/sim/dim_", assayname, D[i], ".pdf")
      ggsave(filepath, p_singvals, width = 60, height = 60, units = "mm", device = "pdf")
    }
  }
}
# Save batch effect-associated singular values (for Figure 3.5C)
saveRDS(batch_singvals, "batch_singvals.rds")
batch_singvalsreadRDS("output/batch_singvals.rds")

# Plot highest batch singval as a function of $d$ (Figure 3.5C)
N <- ncol(sim_0.1)
P <- nrow(sim_0.1)
c <- min(N/P, P/N)
a <- (1 - sqrt(c))^2
b <- (1 + sqrt(c))^2
uncorrected_sv <- metadata(sim_0.1)$logcounts$SVD$d[1] / sqrt(P - 1)

names(batch_singvals) <- c("fastMNN", "Harmony", "Seurat")
batch_singvals$Seurat[8:14] <- NA
data.frame(D, data.frame(batch_singvals)) %>%
  gather(key = "method", value = "singval", -D) %>%
  ggplot(aes(D, singval, color = method)) +
  geom_hline(yintercept = sqrt(c(a,b )), linetype = "dashed", alpha = 0.5) +
  geom_point() + geom_line() +
  geom_hline(yintercept = uncorrected_sv, alpha = 0.5) +
  scale_color_manual(values = c("#A02C2C", "#2C2CA0", "#2CA02C")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    #legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
    legend.position = "none",
    legend.box.background = element_rect(colour = "black")
  ) +
  xlab("d") + ylab(expression(sigma / sqrt(P - 1)))
filepath <- paste0("../plots_report/sim/dim2.pdf")
ggsave(filepath, width = 60, height = 60, units = "mm", device = "pdf")

# Note: Figure 3.5A is the same plot as Figure 3.1A, upper plot.

####################

# Scenario 4 batch correction (Figure 3.6)
set.seed(100)
D <- seq(0, 100, 10)
mixing_sc4 <- NULL # Save mixing metrics
SVDs <- list() # Save all SVDs
for (d in D) {
  if (d > 0) {
    sc4 <- 
      sc4 %>%
      fastmnn_wrapper2(d = d) %>%
      harmony_wrapper(d = d) %>%
      seurat_wrapper(d = d) %>%
      compute_svd("harmony") %>%
      compute_svd("mnn") %>%
      compute_svd("seurat")
    
    # Save SVDs
    for (assayname in c("harmony", "mnn", "seurat")) {
      SVDs[[paste0("d", d)]][[assayname]] <- metadata(sc4)[[assayname]]$SVD
    }
  }
  
  # Compute mixing with 
  for (assayname in c("mnn", "harmony", "seurat")) {
    if (d == 0) {
      method <- "logcounts"
    } else {
      method <- assayname
    }
    mixing_sc4 <- rbind(
      mixing_sc4, 
      data.frame(d = d, method = assayname, compute_metrics(sc4, method, k = 75))
    )
  }
}
# Save computed metrics
#saveRDS(mixing_sc4, "output/mixing_sc4.rds")
#saveRDS(SVDs, "output/SVDs_sc4.rds")

# (Optional) load previously saved computations
mixing_sc4 <- readRDS("output/mixing_sc4.rds")
SVDs <- readRDS("output/SVDs_sc4.rds")

# kBET boxplots (Figure 3.6A)
p_kbet <- 
  mixing_sc4 %>%
  filter(metric == "kBET") %>%
  mutate(d = as.factor(d)) %>%
  mutate(method = factor(
    method, 
    levels = c("mnn", "harmony", "seurat"),
    labels = c("fastMNN", "Harmony", "Seurat")
  )) %>%
  ggplot(aes(d, value, fill = method)) +
  facet_grid(~method) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_fill_manual(values = c("#A02C2C", "#2C2CA0", "#2CA02C")) +
  scale_x_discrete(labels = ifelse(0:10 %% 2, "", as.character((0:10)*10))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    legend.position = "none"
  ) +
  ylab("kBET")

# LISI boxplots (Figure 3.6B)
p_lisi <-
  mixing_sc4 %>%
  filter(metric == "LISI") %>%
  mutate(d = as.factor(d)) %>%
  mutate(method = factor(
    method, 
    levels = c("mnn", "harmony", "seurat"),
    labels = c("fastMNN", "Harmony", "Seurat")
  )) %>%
  ggplot(aes(d, value, fill = method)) +
  facet_grid(~method) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_fill_manual(values = c("#A02C2C", "#2C2CA0", "#2CA02C")) +
  scale_x_discrete(labels = ifelse(0:10 %% 2, "", as.character((0:10)*10))) +
  scale_y_continuous(breaks = c(1, 1.5, 2)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    legend.position = "none"
  ) +
  ylab("LISI")

# local KL boxplots (Figure 3.6C)
p_kl <-
  mixing_sc4 %>%
  filter(metric == "KL") %>%
  mutate(d = as.factor(d)) %>%
  mutate(method = factor(
    method, 
    levels = c("mnn", "harmony", "seurat"),
    labels = c("fastMNN", "Harmony", "Seurat")
  )) %>%
  ggplot(aes(d, value, fill = method)) +
  facet_grid(~method) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_fill_manual(values = c("#A02C2C", "#2C2CA0", "#2CA02C")) +
  scale_x_discrete(labels = ifelse(0:10 %% 2, "", as.character((0:10)*10))) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(), 
    legend.position = "none"
  ) +
  ylab("KL")

# Compute SVDs for d = 0
SVDs$d0 <- replicate(3, metadata(sc4)$logcounts$SVD, simplify = FALSE)
names(SVDs$d0) <- c("harmony", "mnn", "seurat")

# Extract and rescale singvals, and identify BSVs
N <- ncol(sc4)
P <- nrow(sc4)
singvals <- lapply(seq(0, 100, 10), function(i) {
  lapply(c("mnn", "harmony", "seurat"), function(method) {
    SVD <- SVDs[[paste0("d", i)]][[method]]
    singval <- SVD$d / sqrt(max(N,P) - 1)
    p.val <- identify_batch_vectors(SVD$v, batch = sc4$Batch, method = "AD", idx = FALSE)
    data.frame(d = i, method, singval, p.val)
  }) %>%
    do.call(rbind, .)
}) %>% 
  do.call(rbind, .) %>%
  mutate(method = factor(
    method, 
    levels = c("mnn", "harmony", "seurat"),
    labels = c("fastMNN", "Harmony", "Seurat")
  )) %>%
  mutate(BSV = p.val < 0.05)

# BSV plot (Figure 3.6D)
N <- ncol(sc4)
P <- nrow(sc4)
c <- min(N/P, P/N)
bulk_bounds <- c(1 - sqrt(c), 1 + sqrt(c))
p_svals <-
  singvals %>%
  arrange(BSV) %>%  # Make sure BSVs are on top of non_BSVs
  ggplot(aes(x = d, y = singval, color = BSV)) +
  geom_hline(yintercept = bulk_bounds, linetype = "dashed", alpha = 0.5) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_x_continuous(
    breaks = seq(0, 100, 10),
    labels = ifelse(0:10 %% 2, "", as.character((0:10)*10))
  ) +
  facet_grid(~method) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.99, 0.99), legend.justification = c(1, 1)
  ) +
  ylab(expression(sigma / sqrt(P - 1)))

# Combine plots (Figure 3.6)
sc4_bc <- ggarrange(p_kbet, p_lisi, p_kl, p_svals, ncol = 1)
filepath <- paste0(output_dir, "sc4_bc.pdf")
ggsave(filepath, sc4_bc, width = 180, height = 240, units = "mm", device = "pdf")


