# This file contains the analysis of the dendritic cell data.
# It was used for producing the plots in Figures 3.6, 3.7, and 3.8.

library(scater)  # contains the SingleCellExperiment class and accompanying functions
library(scran)  # contains the filterHVGs() function
library(batchelor)  # for fastMNN batch correction
library(harmony)  # for Harmony batch correction
library(Seurat)  # for Seurat batch correction
library(tidyverse)  # for ggplots and dataframe manipulations
library(egg)  # for arranging some ggplots together
library(kSamples)  # for the Anderson-Darling test
library(kBET)  # for kBET batch mixing metric
library(lisi)  # for LISI batch mixing metric
library(FNN)  # fast NN search for local KL divergence mixing metric
source("../helpers.R")  # contains a variety of handmade functions

# Specify output directory for plots
output_dir <- "../plots_report/dendrite/" 

# Load data and filter out low quality cells (low quality cells were determined in preprocess.R)
sce_full <- readRDS("output/dendrite.rds")
sce_full <- sce_full[rowData(sce_full)$QC, colData(sce_full)$QC]

# Normalize
sce_full <- logNormCounts(sce_full)  # almost not necessary bc already CPM normalized

# Select 10% top HVGs
sce <- filter_HVGs(sce_full)

# Compute PCA
reducedDim(sce, "PCA_logcounts") <- calculatePCA(sce, ntop = nrow(sce))

# UMAP of uncorrected data (for Figure 3.7)
set.seed(731)
reducedDim(sce, "UMAP_logcounts") <- calculateUMAP(sce, dimred = "PCA_logcounts", n_dimred = 50)
plot_data(sce, "logcounts", type = "UMAP", size = 0.7)
filepath <- paste0(output_dir, "UMAP_logcounts.pdf")
ggsave(filepath, width = 60, height = 60, units = "mm", device = "pdf")


# Batch mixing metrics for different values of d (Figure 3.8)
set.seed(100)
D <- seq(0, 100, 10)
mixing <- NULL # Save mixing metrics
SVDs <- list() # Save all SVDs
for (d in D) {
  if (d > 0) {
    sce <- 
      sce %>%
      fastmnn_wrapper2(d = d) %>%
      harmony_wrapper(d = d) %>%
      seurat_wrapper(d = d) %>%
      compute_svd("harmony") %>%
      compute_svd("mnn") %>%
      compute_svd("seurat")
    
    # Save SVDs
    for (assayname in c("harmony", "mnn", "seurat")) {
      SVDs[[paste0("d", d)]][[assayname]] <- metadata(sce)[[assayname]]$SVD
    }
  }
  
  # Compute batch mixing
  for (assayname in c("mnn", "harmony", "seurat")) {
    if (d == 0) {
      method <- "logcounts"
    } else {
      method <- assayname
    }
    mixing <- rbind(
      mixing, 
      data.frame(d = d, method = assayname, compute_metrics(sce, method, k = 90))
    )
  }
}

# Save rds objects to disk for time efficiency
saveRDS(mixing, "output/mixing_DC.rds")
saveRDS(SVDs, "output/SVDs_DC.rds")

# Load previously saved metrics (optional)
#mixing <- readRDS("output/mixing_DC.rds")
#SVDs <- readRDS("output/SVDs_DC.rds")

# kBET batch mixing boxplots (Figure 3.8A)
p_kbet <- 
  mixing %>%
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

# LISI batch mixing boxplots (Figure 3.8B)
p_lisi <-
  mixing %>%
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

# local KL divergence batch mixing boxplots (Figure 3.8C)
p_kl <-
  mixing %>%
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
sce <- compute_svd(sce)
SVDs$d0 <- replicate(3, metadata(sce)$logcounts$SVD, simplify = FALSE)
names(SVDs$d0) <- c("harmony", "mnn", "seurat")

# Extract and rescale singular values, and identify BSVs
N <- ncol(sce)
P <- nrow(sce)
singvals <- lapply(seq(0, 100, 10), function(i) {
  lapply(c("mnn", "harmony", "seurat"), function(method) {
    SVD <- SVDs[[paste0("d", i)]][[method]]
    singval <- SVD$d / sqrt(max(N,P) - 1)
    p.val <- identify_batch_vectors(SVD$v, batch = sce$Batch, method = "AD", idx = FALSE)
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

# BSV plots (Figure 3.8D)
N <- ncol(sce)
P <- nrow(sce)
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

# Combine plots together (Figure 3.8)
DC_bc <- ggarrange(p_kbet, p_lisi, p_kl, p_svals, ncol = 1)
filepath <- paste0(output_dir, "DC_bc.pdf")
ggsave(filepath, DC_bc, width = 180, height = 240, units = "mm", device = "pdf")

# Batch correction with selected d values
D <- c(0, 50, 30, 10)
sce <-
  sce %>%
  fastmnn_wrapper2(d = D[2]) %>%
  harmony_wrapper(d = D[3]) %>%
  seurat_wrapper(d = D[4]) %>%
  compute_svd("logcounts") %>%
  compute_svd("mnn") %>%
  compute_svd("harmony") %>%
  compute_svd("seurat")

# Compute and plot UMAPs after batch correction (not included in thesis)
set.seed(113)
for (method in c("logcounts", "mnn", "harmony", "seurat")) {
  reducedDim(sce, paste0("UMAP_", method)) <- calculateUMAP(sce, dimred = paste0("PCA_", method))
  p_umap <- plot_data(sce, method, type = "UMAP", size = 0.7, legend = FALSE)
  filepath <- paste0(output_dir, "sce_UMAP_", method, ".pdf")
  ggsave(filepath, plot = p_umap, width = 60, height = 60, units = "mm", device = "pdf")
}

# Make sv histograms (Figure 3.9, upper row)
singvals_params <- find_plot_params(list(sce), assayname = "mnn", sv_test = "AD")
caption_dict <- list(
  logcounts = "Uncorrected", 
  mnn = "fastMNN", 
  harmony = "Harmony", 
  seurat= "Seurat"
)
annot_x <- max(singvals_params$breaks)
annot_y <- max(singvals_params$ylims)
color_list <- c(
  logcounts = "#A0A0A0", 
  mnn = "#A02C2C", 
  harmony = "#2C2CA0", 
  seurat = "#2CA02C"
)
i <- 0
for (assayname in names(metadata(sce))) {
  i <- i + 1
  p_singvals <- 
    plot_hist(
      sce, assayname, sv_test = "AD", pparams = singvals_params, 
      fill = color_list[[assayname]]
    ) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    annotate(
      "text", x = annot_x, y = annot_y,
      label = paste0(caption_dict[[assayname]], "\n", "d = ", D[i]),
      hjust = "inward", vjust = "inward", size = 3
    )
  print(p_singvals)
  filepath <- paste0(output_dir, "singvals_", assayname, ".pdf")
  ggsave(filepath, p_singvals, width = 60, height = 60, units = "mm", device = "pdf")
}

# Plot histograms of 4th singular vector elements (Figure 3.9, lower row)
for (assayname in names(metadata(sce))) {
  plot_vector(sce, assayname, vectors = 4, bins = 50) + 
    scale_y_continuous(limits = c(0, 61)) +
    theme(legend.position = "none")
  filepath <- paste0(output_dir, "v4_", assayname, ".pdf")
  ggsave(filepath, width = 60, height = 60, units = "mm", device = "pdf")
}