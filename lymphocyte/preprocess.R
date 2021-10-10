# This file contains the code for preprocessing the lymphocyte data.

library(scater)  # for SingleCellExperiment class and associated functions (e.g. for QC)
library(tidyverse)  # for data manipulations and plotting QC metrics

# Read data
dat <- readRDS("data/K2nB_K3nB.rds")

# Cells per subject per batch
table(dat$cellAnn$p, dat$cellAnn$calledS)
table(dat$cellAnn$status)  # already filtered for singlets

# Create SingleCellExperiment object
batch_barcode <- strsplit(dat$cellAnn$c, split = ":")
barcode <- sapply(batch_barcode, function(x) x[2])
coldat <- dat$cellAnn[, c(1, 3)]
names(coldat) <- c("Batch", "Subject")
sce <- SingleCellExperiment(
  assays = list(counts = dat$gex),
  colData = coldat
)

# Identify mitochondrial genes
mt_genes <- read.csv("mt_genes.csv")
mt_HGNC <- 
  mt_genes$alias %>%
  strsplit(", ") %>%
  sapply(function(x) x[length(x)]) %>%
  sub(pattern = "MT", replacement = "MT-")
is_mt <- rownames(sce) %in% mt_HGNC
feature_type <- rep("endogenous", nrow(sce))
feature_type[is_mt] <- "MT"
sce <- splitAltExps(sce, feature_type)

# Compute quality metrics
colData(sce) <- cbind(colData(sce), perCellQCMetrics(sce))

# Save the whole thing
saveRDS(sce, "output/lymphocyte_full.rds")

# Split over subjects
S1 <- sce[, sce$Subject == "Subject1"]
S2 <- sce[, sce$Subject == "Subject2"]


#### Cell QC ####
# Visualize library size per batch
as.data.frame(colData(sce)) %>%
  ggplot(aes(sum)) +
  geom_histogram(bins = 50, color = "#000000", fill = "#999999") +
  facet_wrap(Subject ~ Batch) +
  scale_x_log10() + 
  geom_vline(
    data = filter(as.data.frame(colData(sce)), Batch == "K2nB", Subject == "Subject2"), 
    aes(xintercept = 1000)
  ) +
  geom_vline(
    data = filter(as.data.frame(colData(sce)), Batch == "K2nB", Subject == "Subject1"), 
    aes(xintercept = 1000)
  ) +
  geom_vline(
    data = filter(as.data.frame(colData(sce)), Batch == "K3nB", Subject == "Subject1"), 
    aes(xintercept = 1250)
  ) +
  theme_classic()

# Visualize number of detected genes per batch
as.data.frame(colData(sce)) %>%
  ggplot(aes(detected)) +
  geom_histogram(bins = 50, color = "#000000", fill = "#999999") +
  facet_wrap(Subject ~ Batch) +
  scale_x_log10() + 
  geom_vline(
    data = filter(as.data.frame(colData(sce)), Batch == "K2nB", Subject == "Subject2"), 
    aes(xintercept = 430)
  ) +
  geom_vline(
    data = filter(as.data.frame(colData(sce)), Batch == "K2nB", Subject == "Subject1"), 
    aes(xintercept = 500)
  ) +
  geom_vline(
    data = filter(as.data.frame(colData(sce)), Batch == "K3nB", Subject == "Subject2"), 
    aes(xintercept = 650)
  ) +
  geom_vline(
    data = filter(as.data.frame(colData(sce)), Batch == "K3nB", Subject == "Subject1"), 
    aes(xintercept = 600)
  ) +
  theme_classic()

# Visualize library size, detected genes, and % mitochondrial expression
as.data.frame(colData(sce)) %>%
  ggplot(aes(sum, detected, color = altexps_MT_percent)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  facet_wrap(Subject ~ Batch) +
  scale_x_log10() + 
  scale_y_log10() +
  geom_vline(
    data = filter(as.data.frame(colData(sce)), Batch == "K2nB", Subject == "Subject1"), 
    aes(xintercept = 1000)
  ) +
  geom_vline(
    data = filter(as.data.frame(colData(sce)), Batch == "K2nB", Subject == "Subject2"), 
    aes(xintercept = 1000)
  ) +
  geom_vline(
    data = filter(as.data.frame(colData(sce)), Batch == "K3nB", Subject == "Subject1"), 
    aes(xintercept = 1250)
  ) +
  geom_hline(
    data = filter(as.data.frame(colData(sce)), Batch == "K2nB", Subject == "Subject2"), 
    aes(yintercept = 650)
  ) +
  geom_hline(
    data = filter(as.data.frame(colData(sce)), Batch == "K2nB", Subject == "Subject1"), 
    aes(yintercept = 500)
  ) +
  geom_hline(
    data = filter(as.data.frame(colData(sce)), Batch == "K3nB", Subject == "Subject2"), 
    aes(yintercept = 650)
  ) +
  geom_hline(
    data = filter(as.data.frame(colData(sce)), Batch == "K3nB", Subject == "Subject1"), 
    aes(yintercept = 600)
  ) +
  theme_classic()

# Define cutoff thresholds
cut_sum_S1 <- c(1000, 1250)
cut_sum_S2 <- c(1000, 0)
cut_detect_S1 <- c(500, 600)
cut_detect_S2 <- c(650, 650)
qc_S1_B1 <- S1$sum >= cut_sum_S1[1] & S1$detected >= cut_detect_S1[1]
qc_S1_B2 <- S1$sum >= cut_sum_S1[2] & S1$detected >= cut_detect_S1[2]
qc_S2_B1 <- S2$sum >= cut_sum_S2[1] & S2$detected >= cut_detect_S2[1]
qc_S2_B2 <- S2$sum >= cut_sum_S2[2] & S2$detected >= cut_detect_S2[2]
S1$QC <- ifelse(S1$Batch == "K2nB", qc_S1_B1, qc_S1_B2)
S2$QC <- ifelse(S2$Batch == "K2nB", qc_S2_B1, qc_S2_B2)

# Visualize library size, detected genes, and % MT expression after QC
as.data.frame(colData(S2[, S2$QC])) %>%
  ggplot(aes(sum, detected, color = altexps_MT_percent)) +
  geom_point(size = 0.8) +
  scale_colour_gradient(low = "gray90", high = "black") +
  facet_wrap(~ Batch) +
  scale_x_log10() + 
  scale_y_log10() +
  theme_classic()

#### Gene QC ####
# Keep genes with a count of >1 in at least 2 cells
rowData(S1)$QC <- nexprs(
  S1[, S1$QC],
  byrow = TRUE,
  detection_limit = 1
) >= 2
rowData(S2)$QC <- nexprs(
  S2[, S2$QC],
  byrow = TRUE,
  detection_limit = 1
) >= 2

# Save all
saveRDS(S1, "output/Subject1_full.rds")
saveRDS(S2, "output/Subject2_full.rds")

# Normalization (library size normalization)
S1 <- S1[rowData(S1)$QC, colData(S1)$QC]
S2 <- S2[rowData(S2)$QC, colData(S2)$QC]
S1 <- logNormCounts(S1)
S2 <- logNormCounts(S2)

# Save separately
saveRDS(S1, "output/Subject1.rds")
saveRDS(S2, "output/Subject2.rds")
