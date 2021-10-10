# This file contains the analysis of the lymphocyte data.
# It was used for producing the plots in Figures 3.10-3.14.

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
library(FNN)  # fast NN search for local KL divergence metric
source("../helpers.R")  # contains a variety of handmade functions

# Specify output directory for plots
output_dir <- "../plots_report/lympho/"

# Read data
S1_full <- readRDS("output/Subject1.rds")
S2_full <- readRDS("output/Subject2.rds")

# Select top 20% HVGs
S1 <- filter_HVGs(S1_full, percent_topHVGs = 0.2)
S2 <- filter_HVGs(S2_full, percent_topHVGs = 0.2)

# Compute PCA
reducedDim(S1, "PCA_logcounts") <- calculatePCA(S1, ntop = nrow(S1))
reducedDim(S2, "PCA_logcounts") <- calculatePCA(S2, ntop = nrow(S2))


# UMAPs of uncorrected data (for Figure 3.10)
set.seed(731)
reducedDim(S1, "UMAP_logcounts") <- calculateUMAP(S1, dimred = "PCA_logcounts", n_dimred = 50)
plot_data(S1, "logcounts", type = "UMAP", size = 0.7, legend = FALSE)
filepath <- paste0(output_dir, "S1_UMAP_logcounts.pdf")
ggsave(filepath, width = 60, height = 60, units = "mm", device = "pdf")

reducedDim(S2, "UMAP_logcounts") <- calculateUMAP(S2, dimred = "PCA_logcounts", n_dimred = 50)
plot_data(S2, "logcounts", type = "UMAP", size = 0.7, legend = FALSE)
filepath <- paste0(output_dir, "S2_UMAP_logcounts_v2.pdf")
ggsave(filepath, width = 60, height = 60, units = "mm", device = "pdf")

#### Patient 1
# Compute batch mixing metrics for Patient 1
set.seed(100)
D <- seq(0, 100, 10)
mixing_S1 <- NULL # Save mixing metrics
SVDs_S1 <- list() # Save all SVDs
for (d in D) {
  if (d > 0) {
    S1 <- 
      S1 %>%
      fastmnn_wrapper2(d = d) %>%
      harmony_wrapper(d = d) %>%
      seurat_wrapper(d = d) %>%
      compute_svd("harmony") %>%
      compute_svd("mnn") %>%
      compute_svd("seurat")
    
    # Save SVDs
    for (assayname in c("harmony", "mnn", "seurat")) {
      SVDs_S1[[paste0("d", d)]][[assayname]] <- metadata(S1)[[assayname]]$SVD
    }
  }
  
  # Compute mixing with 
  for (assayname in c("mnn", "harmony", "seurat")) {
    if (d == 0) {
      method <- "logcounts"
    } else {
      method <- assayname
    }
    mixing_S1 <- rbind(
      mixing_S1, 
      data.frame(d = d, method = assayname, compute_metrics(S1, method, k = 90))
    )
  }
}

# Save mixing metrics
saveRDS(mixing_S1, "output/mixing_S1.rds")
saveRDS(SVDs_S1, "output/SVDs_S1.rds")

# Load previously saved metrics (optional)
#mixing <- readRDS("output/mixing_DC.rds")
#SVDs <- readRDS("output/SVDs_DC.rds")

# Patient 1 kBET batch mixing boxplots (Figure 3.11A)
p_kbet_S1 <- 
  mixing_S1 %>%
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

# Patient 1 LISI batch mixing boxplots (Figure 3.11B)
p_lisi_S1 <-
  mixing_S1 %>%
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

# Patient 1 LISI batch mixing boxplots (Figure 3.11C)
p_kl_S1 <-
  mixing_S1 %>%
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
S1 <- compute_svd(S1)
SVDs_S1$d0 <- replicate(3, metadata(S1)$logcounts$SVD, simplify = FALSE)
names(SVDs_S1$d0) <- c("harmony", "mnn", "seurat")

# Extract and rescale singular values, and identify BSVs
N <- ncol(S1)
P <- nrow(S1)
singvals_S1 <- lapply(seq(0, 100, 10), function(i) {
  lapply(c("mnn", "harmony", "seurat"), function(method) {
    SVD <- SVDs_S1[[paste0("d", i)]][[method]]
    singval <- SVD$d / sqrt(max(N,P) - 1)
    p.val <- identify_batch_vectors(SVD$v, batch = S1$Batch, method = "AD", idx = FALSE)
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

# BSV plots (Figure 3.11D)
c <- min(N/P, P/N)
bulk_bounds <- c(1 - sqrt(c), 1 + sqrt(c))
p_svals_S1 <-
  singvals_S1 %>%
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
    legend.position = "none"
  ) +
  ylab(expression(sigma / sqrt(P - 1)))

# Combine plots together (Figure 3.11)
S1_bc <- ggarrange(p_kbet_S1, p_lisi_S1, p_kl_S1, p_svals_S1, ncol = 1)
filepath <- paste0(output_dir, "S1_bc.pdf")
ggsave(filepath, S1_bc, width = 180, height = 240, units = "mm", device = "pdf")


# Make sv histograms (Figure 3.12, upper row)
D <- c(logcounts = 0, mnn = 50, harmony = 50, seurat = 40)
set.seed(234)
S1 <-
  S1 %>%
  fastmnn_wrapper2(d = D[2]) %>%
  harmony_wrapper(d = D[3]) %>%
  seurat_wrapper(d = D[4]) %>%
  compute_svd("logcounts") %>%
  compute_svd("mnn") %>%
  compute_svd("harmony") %>%
  compute_svd("seurat")

# Compute UMAPs
set.seed(112)
reducedDim(S1, "UMAP_mnn") <- calculateUMAP(S1, dimred = "PCA_mnn")
reducedDim(S1, "UMAP_harmony") <- calculateUMAP(S1, dimred = "PCA_harmony")
reducedDim(S1, "UMAP_seurat") <- calculateUMAP(S1, dimred = "PCA_seurat")

# Compute and plot UMAPs after batch correction (not included in thesis)
for (method in c("mnn", "harmony", "seurat")) {
  reducedDim(S1, paste0("UMAP_", method)) <- calculateUMAP(S1, dimred = paste0("PCA_", method))
  p_umap <- plot_data(S1, method, type = "UMAP", size = 0.5, legend = FALSE)
  filepath <- paste0(output_dir, "S1_UMAP_", method, ".pdf")
  ggsave(filepath, plot = p_umap, width = 60, height = 60, units = "mm", device = "pdf")
}

# Make sv histograms (Figure 3.12, upper row)
singvals_params <- find_plot_params(list(S1), assayname = "seurat", sv_test = "AD")
caption_dict <- list(
  logcounts = "Uncorrected", 
  mnn = "fastMNN",
  harmony = "Harmony", 
  seurat = "Seurat"
)
annot_x <- max(singvals_params$breaks)
annot_y <- max(singvals_params$ylims)
color_list <- c(
  logcounts = "#A0A0A0", 
  mnn = "#A02C2C", 
  harmony = "#2C2CA0", 
  seurat = "#2CA02C"
)
for (assayname in names(metadata(S1))) {
  p_singvals <- 
    plot_hist(
      S1, assayname, sv_test = "AD", 
      pparams = singvals_params, fill = color_list[[assayname]]
    ) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    annotate(
      "text", x = annot_x, y = annot_y,
      label = paste0(caption_dict[[assayname]], "\n", "d = ", D[[assayname]]),
      hjust = "inward", vjust = "inward", size = 3
    )
  print(p_singvals)
  filepath <- paste0(output_dir, "S1_singvals_", assayname, ".pdf")
  ggsave(filepath, p_singvals, width = 60, height = 60, units = "mm", device = "pdf")
}

# Plot histograms of 1st singular vector elements (Figure 3.12, lower row)
for (assayname in names(metadata(S1))) {
  plot_vector(S1, assayname, vectors = 1, bins = 50) + 
    scale_x_continuous(breaks = c(-0.03, 0, 0.03)) +
    scale_y_continuous(limits = c(0, 88)) +
    theme(legend.position = "none")
  filepath <- paste0(output_dir, "S1_v1_", assayname, ".pdf")
  ggsave(filepath, width = 60, height = 60, units = "mm", device = "pdf")
}

#### Patient 2
# Batch mixing metrics for different values of d (Figure 3.13)
set.seed(100)
D <- seq(0, 100, 10)
mixing_S2 <- NULL # Save mixing metrics
SVDs_S2 <- list() # Save all SVDs
for (d in D) {
  if (d > 0) {
    S2 <- 
      S2 %>%
      fastmnn_wrapper2(d = d) %>%
      harmony_wrapper(d = d) %>%
      seurat_wrapper(d = d) %>%
      compute_svd("harmony") %>%
      compute_svd("mnn") %>%
      compute_svd("seurat")
    
    # Save SVDs
    for (assayname in c("harmony", "mnn", "seurat")) {
      SVDs_S2[[paste0("d", d)]][[assayname]] <- metadata(S2)[[assayname]]$SVD
    }
  }
  
  # Compute mixing metrics
  for (assayname in c("mnn", "harmony", "seurat")) {
    if (d == 0) {
      method <- "logcounts"
    } else {
      method <- assayname
    }
    mixing_S2 <- rbind(
      mixing_S2, 
      data.frame(d = d, method = assayname, compute_metrics(S2, method, k = 90))
    )
  }
}
# Save metrics 
saveRDS(mixing_S2, "output/mixing_S2_v2.rds")
saveRDS(SVDs_S2, "output/SVDs_S2_v2.rds")

# Load previously saved metrics (optional)
#mixing <- readRDS("output/mixing_S2_v2.rds")
#SVDs <- readRDS("output/SVDs_S2_v2.rds")

# Patient 2 kBET batch mixing boxplots (Figure 3.13A)
p_kbet_S2 <- 
  mixing_S2 %>%
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

# Patient 2 LISI batch mixing boxplots (Figure 3.13B)
p_lisi_S2 <-
  mixing_S2 %>%
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

# Patient 2 local KL divergence boxplots (Figure 3.13C)
p_kl_S2 <-
  mixing_S2 %>%
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
S2 <- compute_svd(S2)
SVDs_S2$d0 <- replicate(3, metadata(S2)$logcounts$SVD, simplify = FALSE)
names(SVDs_S2$d0) <- c("harmony", "mnn", "seurat")

# Extract and rescale singular values, and identify BSVs
N <- ncol(S2)
P <- nrow(S2)
singvals_S2 <- lapply(seq(0, 100, 10), function(i) {
  lapply(c("mnn", "harmony", "seurat"), function(method) {
    SVD <- SVDs_S2[[paste0("d", i)]][[method]]
    singval <- SVD$d / sqrt(max(N,P) - 1)
    p.val <- identify_batch_vectors(SVD$v, batch = S2$Batch, method = "AD", idx = FALSE)
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

# BSV plots (Figure 3.13D)
c <- min(N/P, P/N)
bulk_bounds <- c(1 - sqrt(c), 1 + sqrt(c))
p_svals_S2 <-
  singvals_S2 %>%
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
    legend.position = "none"
  ) +
  ylab(expression(sigma / sqrt(P - 1)))

# Combine plots together (Figure 3.13)
S2_bc <- ggarrange(p_kbet_S2, p_lisi_S2, p_kl_S2, p_svals_S2, ncol = 1)
filepath <- paste0(output_dir, "S2_bc_v2.pdf")
ggsave(filepath, S2_bc, width = 180, height = 240, units = "mm", device = "pdf")

# Batch correction with selected d values
set.seed(1234)
D <- c(logcounts = 0, mnn = 50, harmony = 30, seurat = 10)
S2 <-
  S2 %>%
  fastmnn_wrapper2(d = D[2]) %>%
  harmony_wrapper(d = D[3]) %>%
  seurat_wrapper(d = D[4]) %>%
  compute_svd("logcounts") %>%
  compute_svd("mnn") %>%
  compute_svd("harmony") %>%
  compute_svd("seurat")

# Compute and plot UMAPs after batch correction (not included in thesis)
set.seed(712)
for (method in c("mnn", "harmony", "seurat")) {
  reducedDim(S2, paste0("UMAP_", method)) <- calculateUMAP(S2, dimred = paste0("PCA_", method))
  p_umap <- plot_data(S2, method, type = "UMAP", size = 0.5, legend = FALSE)
  filepath <- paste0(output_dir, "S2_UMAP_", method, ".pdf")
  ggsave(filepath, plot = p_umap, width = 60, height = 60, units = "mm", device = "pdf")
}

# Make sv histograms (Figure 3.14, upper row)
singvals_params <- find_plot_params(list(S2), assayname = "seurat", sv_test = "AD")
annot_x <- max(singvals_params$breaks)
annot_y <- max(singvals_params$ylims)
caption_dict <- list(
  logcounts = "Uncorrected", 
  mnn = "fastMNN",
  harmony = "Harmony", 
  seurat = "Seurat"
)
color_list <- c(
  logcounts = "#A0A0A0", 
  mnn = "#A02C2C", 
  harmony = "#2C2CA0", 
  seurat = "#2CA02C"
)
for (assayname in names(metadata(S2))) {
  p_singvals <- 
    plot_hist(S2, assayname, sv_test = "AD", pparams = singvals_params, fill = color_list[[assayname]]) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    annotate(
      "text", x = annot_x, y = annot_y,
      label = paste0(caption_dict[[assayname]], "\n", "d = ", D[[assayname]]),
      hjust = "inward", vjust = "inward", size = 3
    )
  print(p_singvals)
  filepath <- paste0(output_dir, "S2_singvals_", assayname, "_v2.pdf")
  ggsave(filepath, p_singvals, width = 60, height = 60, units = "mm", device = "pdf")
}

# Plot histograms of 1st singular vector elements (Figure 3.14, lower row)
for (assayname in names(metadata(S2))) {
  plot_vector(S2, assayname, vectors = 1, bins = 50) + 
    scale_x_continuous(breaks = c(-0.03, 0, 0.03)) +
    scale_y_continuous(limits = c(0, 100)) +
    theme(legend.position = "none")
  filepath <- paste0(output_dir, "S2_v1_", assayname, "_v2.pdf")
  ggsave(filepath, width = 60, height = 60, units = "mm", device = "pdf")
}

