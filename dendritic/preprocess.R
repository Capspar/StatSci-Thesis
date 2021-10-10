# This file contains the code for preprocessing the dendritic cell data.

library(scater)  # for SingleCellExperiment class and associated functions (e.g. for QC)
library(tidyverse)  # for data manipulations

# Read data
# This data file was downloaded from: 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94820
raw <- as.matrix(
  read.delim("data/GSE94820_raw.expMatrix_DCnMono.discovery.set.submission.txt")
)

# Prepare column data from cell identifiers
coldat <- colnames(raw)
coldat <- strsplit(coldat, "_")
coldat <- as.data.frame(lapply(1:3, function(i) sapply(coldat, function(x) x[i])))
names(coldat) <- c("Celltype", "Plate", "Cell_ID")

# Keep plate cells
keep <- startsWith(coldat$Plate, "P")
coldat <- coldat[keep, ]

# Identify batches
coldat$Batch <- ifelse(
  coldat$Plate %in% c("P7", "P8", "P9", "P10"), 
  yes = "Batch1", no = "Batch2"
)

# Row data
rowdat <- rownames(raw)

# Create SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(counts = raw[, keep]),
  colData = coldat,
  rowData = rowdat
)

# Cell QC
# filter out cells with low detected genes
colData(sce) <- cbind(colData(sce), perCellQCMetrics(sce))
hist(log(sce$detected + 1), breaks = 100)
colData(sce)$QC <- log(sce$detected + 1) > 8

# Gene QC
# Keep genes with a count of >1 in at least 2 cells
keep_feature <- nexprs(
  sce[, sce$QC],
  byrow = TRUE,
  detection_limit = 1
) >= 2
rowData(sce)$QC <- keep_feature

saveRDS(sce, "output/dendrite.rds")