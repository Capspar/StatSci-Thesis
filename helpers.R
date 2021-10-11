# This file contains functions for the simulation and analysis of
# single cell data. The functions are imported in the analysis code of the
# simulations, the dendritic cell data, and the lymphocyte data.

# Prepare simulation scenarios
# Includes:
# - Removal of undetected genes
# - HVG selection
# - PCA
sim_prep <- function(params, top_HVGs = 0.10) {
  
  message("Simulating SingleCellExperiment...")
  if (getParam(params, "nGroups") == 1) {  # To get rid of warning message
    method = "single"
  } else {
    method = "groups"
  }
  
  sim <- splatSimulate(params, method = method, verbose = FALSE)
  
  message("Minimizing SingleCellExperiment...")
  sim <- minimiseSCE(sim, colData.keep = c("Batch", "Group"), verbose = FALSE)
  
  message("Removing undetected genes and lognormalizing...")
  sim <- sim[nexprs(sim, byrow = TRUE, detection_limit = 1) >= 2, ]
  sim <- logNormCounts(sim)  # lognormalize
  
  sim <- filter_HVGs(sim, top_HVGs)
  P <- nrow(sim)
  message(paste("Selected", P, "highly variable genes."))
  
  message("Computing PCA...")
  reducedDim(sim, "PCA_logcounts") <- calculatePCA(sim, ntop = P)
  
  sim
}

# Marchenko-Pastur density function
dmp <- function(x, N, P, makenan = FALSE) {
  c <- min(N/P, P/N)
  a <- (1 - sqrt(c))^2
  b <- (1 + sqrt(c))^2
  
  d <- ifelse(
    (x >= sqrt(a)) & (x <= sqrt(b)),
    yes = sqrt((b - x^2) * (x^2 - a)) / (pi * x * c),
    no = 0
  )
  
  d
}

# (Deprecated) Plot singular value histogram + MP fit (base R plotting)
sv_hist <- function(scaled, batch = NULL, threshold = 0.05, breaks = NULL, y = 1, ...) {
  P <- dim(scaled)[1]
  N <- dim(scaled)[2]
  
  SVD <- svd(scaled)
  singvals <- SVD$d / sqrt(max(N, P) - 1)
  
  if (is.null(breaks)) {
    breaks <- 100
  }
  
  hist(singvals, prob = TRUE, breaks = breaks, ...)
  
  if (!is.null(batch)){
    idx <- identify_batch_vectors(SVD$v, batch, threshold = 0.05)
    text(x = singvals[idx], y = y,"*")
  }
  suppressWarnings(curve(dmp(x, N, P, makenan = TRUE), add = TRUE, lwd = 1.5, n = 1001))
}

# Scale rows, then columns (mean 0, var 1)
scale_row_col <- function(X, return.list = FALSE) {
  row_scl <- t(scale(t(X)))
  row_col_scl <- scale(row_scl)
  scaled_mat <- row_col_scl
  
  if (return.list == FALSE) {
    return(scaled_mat)
    
  # (Not used) Optionally, return scaling factors
  } else {
    attr(scaled_mat, "scaled:scale") <- NULL
    attr(scaled_mat, "scaled:center") <- NULL
    list(
      scaled = scaled_mat,
      gene_center = attr(row_scl, "scaled:center"),
      gene_scale = attr(row_scl, "scaled:scale"),
      cell_center = attr(row_col_scl, "scaled:center"),
      cell_scale = attr(row_col_scl, "scaled:scale")
    )
  }
}

# (Not used) Scale back from scale_row_col to original
backscale <- function(scaled_list) {
  out <- t(t(scaled_list$scaled) * scaled_list$cell_scale + scaled_list$cell_center)
  out * scaled_list$gene_scale + scaled_list$gene_center
}

# Perform t-test given matrix of singular vectors (V) and a vector of batch labels
t_test <- function(V, batch) {
  batchnames <- unique(batch)
  N <- tabulate(factor(batch))
  
  splitted <- lapply(batchnames, function(b) V[batch == b, ])
  means <- sapply(splitted, colMeans)
  vars <- sapply(splitted, function(X) apply(X, 2, var))  # variances
  sq_SEs <- t(t(vars)/N)  # squared SEs
  s <- sqrt(rowSums(sq_SEs))
  
  tstat <- (means[, 1] - means[, 2]) / s
  dfs <- (colSums(t(sq_SEs^2) / (N-1)) / s^4)^-1  # Welch-Satterthwaite approximation of df
  
  pvals <- 2*pt(-abs(tstat), dfs)
  
  data.frame(
    t.statistic = tstat,
    df = dfs,
    p.value = pvals
  )
}

# Perform Anderson-Darling test given singular vectors V and a vector of batch labels
# (Currently only for 2 batches)
# This function uses the function ad.pval from the kSamples package
ad_test <- function(V, batch) {
  batchnames <- unique(batch)
  k <- length(batchnames)
  ns <- tabulate(factor(batch))
  N <- sum(ns)
  Z <- apply(V, 2, sort)
  X <- V[batch == batchnames[1], ]
  
  # For each singular vector, compute number of values of X <= values in Z
  stepfuns <- apply(X, 2, function(x) stepfun(sort(x), 0:length(x)))
  M <- sapply(1:ncol(V), function(i) stepfuns[[i]](Z[, i]))
  
  # Compute Anderson-Darling statistic
  AD <- colSums((M[-N, ]*N - ns[1]*(1:(N-1)))^2 / (1:(N-1)*(N - 1:(N-1)))) / (ns[1]*ns[2])
  
  # Compute variance of AD statistic under H0 (from kSamples::ad.test)
  H <- sum(1/ns)
  h <- sum(1/(1:(N - 1)))
  g <- 0
  for (i in 1:(N - 2)) {
    g <- g + (1/(N - i)) * sum(1/((i + 1):(N - 1)))
  }
  coef.a <- (4 * g - 6) * (k - 1) + (10 - 6 * g) * H
  coef.b <- (2 * g - 4) * k^2 + 8 * h * k + (2 * g - 14 * 
                                               h - 4) * H - 8 * h + 4 * g - 6
  coef.c <- (6 * h + 2 * g - 2) * k^2 + (4 * h - 4 * g + 
                                           6) * k + (2 * h - 6) * H + 4 * h
  coef.d <- (2 * h + 6) * k^2 - 4 * h * k
  sig2 <- (coef.a * N ^3 + coef.b * N ^2 + coef.c * N + coef.d)/((N - 
                                                                    1) * (N - 2) * (N - 3))
  sig <- sqrt(sig2)
  
  # Standardize AD statistics using the estimated sd
  AD_std <- (AD - (k - 1)) / sig
  
  # Compute pvals with standardized AD statistics
  pvals <- kSamples::ad.pval(AD_std, k-1, 1)
  
  data.frame(
    AD.standard = AD_std,
    p.value = pvals
  )
}

# Identify batch effect-associated singular values (BSVs)
identify_batch_vectors <- function(vectors, batch, threshold = 0.05, 
  idx = TRUE, adjust = TRUE, method = "t") {

  if (method == "t") {
    pvals <- t_test(vectors, batch)$p.value
  } else if (method == "AD") {
    pvals <- ad_test(vectors, batch)$p.value
  }
  
  if (adjust) {
    pvals <- p.adjust(pvals, method = "BH")  # adjust for multiple testing by controlling the FDR
  }
  # Return indices of BSVs
  if (idx) {
    return(which(pvals < threshold))
  }
  # ... or the p-values
  else {
    pvals
  }
}

# Extract histogram breaks from a ggplot
get_breaks <- function(p) {
  hist_dat <- ggplot_build(p)$data[[1]]
  breaks <- unique(c(hist_dat$xmin[1], hist_dat$xmax))
  
  breaks
}

# Extract x- or y-limits from a ggplot
get_lims <- function(p, axis) layer_scales(p)[[axis]]$get_limits()

# Plot right singular vector elements
# If length(vectors) == 1, plot a histogram
# If length(vectors) == 2, plot a scatterplot
plot_vector <- function(sce, assayname = "logcounts", vectors, 
                        batch = "Batch", group = "Group", 
                        color = "Batch", shape = "Group", 
                        bins = 100, pparams = list(breaks = NULL), position = "identity") {
  V <- metadata(sce)[[assayname]]$SVD$v
  df <- data.frame(V)
  
  df$Batch <- sce[[batch]]
  df$Group <- sce[[group]]
  
  if (is.null(df$Batch)) Batch <- NULL
  if (is.null(df$Group)) Group <- NULL
  if (length(vectors) == 2) {
    p <- 
      scramble(df, rows = TRUE) %>%
      ggplot(aes_string(
        x = paste0("X", vectors[1]), y = paste0("X", vectors[2]),
        color = color, shape = shape
      )) + 
      geom_point(size = 1.2, alpha = 0.7) +
      scale_color_brewer(palette = "Set1") + 
      xlab(bquote(v[.(vectors[1])])) +
      ylab(bquote(v[.(vectors[2])])) +
      coord_cartesian(xlim = pparams$xlims, ylim = pparams$ylims)
      
  } else {
    p <- 
      df %>%
      ggplot(aes_string(x = paste0("X", vectors), fill = color)) +
      geom_histogram(bins = bins, breaks = pparams$breaks, alpha = 0.8, position = position) +
      scale_fill_brewer(palette = "Set1") + 
      xlab(bquote(v[.(vectors[1])]))
    #if (!is.null(xticks)) {
    #  p <- p + scale_x_continuous(breaks = xticks)
    #}
  }
  
  xticks <- round(c(pparams$xlims[1]*3/4, 0, pparams$xlims[2]*3/4), 2)  # choose 3 xticks away from border
  
  p +
    theme_bw() +
    theme(
      panel.grid = element_blank(), 
      legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
      legend.box.background = element_rect(colour = "black")
    )# +
    #guides(col = guide_legend(ncol = length(unique(sce$Batch)))) +
    #scale_x_continuous(breaks = xticks)
}




compute_svd <- function(sce, assayname = "logcounts") {
  SVD <- svd(scale_row_col(assay(sce, assayname)))
  metadata(sce)[[assayname]]$SVD <- SVD
  
  sce
}


# Plot singular value histogram + theoretical MP density (ggplot)
plot_hist <- function(sce, assayname = "logcounts", bins = NULL, pparams = list(), 
                      sv_test = NULL, threshold = 0.05, fill = "#666666") {
  P <- nrow(sce)
  N <- ncol(sce)
  c <- min(N/P, P/N)
  a <- (1 - sqrt(c))^2 # for the min boundary of MP
  b <- (1 + sqrt(c))^2  # for the max boundary of MP
  
  if (N > P) {
    maxNP <- "N"
  } else{
    maxNP <- "P"
  }
  
  # Extract SVD (has to be computed first with compute_svd())
  SVD <- metadata(sce)[[assayname]]$SVD
  if (is.null(SVD)) stop(paste0("No SVD was computed for ", assayname))
  
  # Scale singular values
  singvals <- SVD$d / sqrt(max(N, P)-1)
  
  if (is.null(bins)) bins <- 100
  p <- 
    data.frame(singvals) %>%
    ggplot(aes(x = singvals)) +
    geom_histogram(
      aes(y = ..density..), 
      fill = fill, bins = bins, breaks = pparams$breaks
    ) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    stat_function(
      fun = dmp, 
      args = list(N = N, P = P),
      n = 1001,
      xlim = sqrt(c(a, b)),
      color = 'black',
      size = 0.8
    ) + 
    xlab(bquote(sigma / sqrt(.(maxNP) - 1))) +
    coord_cartesian(ylim = pparams$ylims)
    #scale_y_continuous(breaks = pparams$yticks)  # will be done when saving the plots
  
  # Plot asterisks at BSV locations
  if (!is.null(sv_test)) {
    idx <- identify_batch_vectors(
      vectors = SVD$v[, -ncol(SVD$v)], batch = sce$Batch, 
      adjust = TRUE, idx = TRUE, threshold = threshold,
      method = sv_test
    )
    if (length(idx) > 0) {
      p <-
        p +
        geom_text(
          aes(x = x, y = -0.1, label = "*"), 
          data = data.frame(x = singvals[idx]),
          size = 5
        )
    }
  }
  
  p
}

# Randomly shuffle rows or columns of an array (mostly for plotting)
scramble <- function(array, rows = FALSE) {
  if (rows) {
    array[sample(nrow(array)), ]
  } else {
    array[, sample(ncol(array))]
  }
}

# Compute suitable x- and y-limits for singular value histograms, such that all
# histograms have the same limits.
find_plot_params <- function(sce_list, assayname = "logcounts", plot_type = "singvals", sv_test = FALSE, vectors = 1) {
  n <- length(sce_list)
  best_breaks <- c(0, 0)
  best_xlims <- c(0, 0)
  best_ylims <- c(0, 0)
  for (i in 1:n) {
    if (plot_type == "singvals") {
      p <- plot_hist(sce_list[[i]], assayname, sv_test = sv_test)
    } else if (plot_type == "singvect") {
      p <- plot_vector(sce_list[[i]], assayname, vectors[[i]])
    } else {
      stop("Invalid plot_type")
    }
    
    if (plot_type == "singvals" | length(vectors[[i]]) == 1) {
      breaks <- get_breaks(p)
    } else {
      breaks <- c(0, 0)
    }
    
    xlims <- get_lims(p, "x")
    ylims <- get_lims(p, "y")
    
    if (diff(xlims) > diff(best_xlims)) {
      best_xlims <- xlims
      best_breaks <- breaks
    }

    if (diff(ylims) > diff(best_ylims)) {
      best_ylims <- ylims
    }
  }
  
  list(
    breaks = best_breaks,
    xlims = best_xlims,
    ylims = best_ylims
  )
}

# Plot dimension reduced representation of expression matrix (PCA or UMAP)
# PCA or UMAP has to be present in the SingleCellExperiment as "PCA_[assayname]"
# or "UMAP_[assayname]"
plot_data <- function(sce, assayname, type = "PCA", alpha = 0.7, size = 1, legend = TRUE) {  # extended to also include UMAP plots lol
  if (type == "PCA")  {
    df <- data.frame(reducedDim(sce, paste0("PCA_", assayname)))
    colname <- "PC"
  } else if (type == "UMAP") {
    df <- data.frame(reducedDim(sce, paste0("UMAP_", assayname)))
    names(df)[1:2] <- paste0("UMAP", 1:2)
    colname <- "UMAP"
  }
  
  df$Batch <- sce$Batch
  df$Group <- sce$Group
  
  if (is.null(df$Batch)) Batch <- NULL
  if (is.null(df$Group)) Group <- NULL
  
  p <- scramble(df, rows = TRUE) %>%
    ggplot(aes_string(x = paste0(colname, 1), y = paste0(colname, 2), color = "Batch", shape = "Group")) +
    geom_point(size = size, alpha = alpha) +
    scale_color_brewer(palette = "Set1") + 
    theme_bw() +
    theme(panel.grid = element_blank())
  
  if (legend == FALSE) {
    return(p + theme(legend.position = "none"))
  } else {
    p + theme(
      legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
      legend.box.background = element_rect(colour = "black")
    )
  }
}

# (Not used) plot all t- or AD-test p-values for batch effect association
plot_pvals <- function(SVD, batch, adjust = FALSE) {
  pvals <- identify_batch_vectors(SVD$v, batch, adjust = adjust, idx = FALSE)
  scale_factor <- max(c(dim(SVD$u), dim(SVD$v)))
  singvals <- SVD$d / sqrt(scale_factor - 1)
  p_pvals <- data.frame(pvalue = pvals, singvals) %>%
    ggplot(aes(x = singvals, y = pvals)) +
    geom_point(size = 0.2) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    ylab("p-value") +
    theme_bw() +
    theme(
      panel.grid.major.y = element_line(), panel.grid.major.x = element_blank(), 
      panel.grid.minor.x = element_blank(), axis.title.x = element_blank(), 
      axis.text.x = element_blank()
    )
  
  p_pvals
}

# Compute KL divergence 
# P = global batch size proportions
# Q = local batch size proprtions
KL <- function(P, Q) sum(log((P/Q)^P))
KL2 <- function(P, Q) rowSums(log((P/Q)^P))  # faster version

# (Not used in favor of FNN package) Compute KNNs among vector elements
knn_1d <- function(x, k) {
  n <- length(x)
  dist <- matrix(0, n, n)
  for (i in 1:n) {
    dist[, i] <- x - x[i]
  }
  dist <- abs(dist)
  
  knn <- apply(dist, 1, order)[1:k, ]
  return(knn)
}

# Compute local KL divergence based on kNN matrix
local_KL <- function(knns, labels) {
  labels <- as.numeric(as.factor(labels))
  if (nrow(knns) > ncol(knns)) knns <- t(knns)  # to get knns as columns
  n <- length(labels)
  k <- nrow(knns)
  # Get local labels
  locals <- matrix(labels[knns], nrow = k)
  
  p <- tabulate(labels) / n
  m <- length(p)
  P <- matrix(rep(p, each = n), ncol = m)
  if (m < 10) {  # speed consideration
    Q <- sapply(1:m, function(x) colMeans(locals == x))
  } else {
    Q <- t(apply(locals, 2, function(x) tabulate(x, nbins = m)) / k)
  }
  KLs <- KL2(Q, P)
  KLs
}

# Compute kNNs and compute local KL metrics
compute_localKL <- function(sce, assayname, k) {
  X <- reducedDim(sce, paste0("PCA_", assayname))
  knns <- FNN::get.knn(X, k = k)$nn.index
  KLs <- local_KL(knns, sce$Batch)
  
  KLs
}

# Compute kBET metric
compute_kbet <- function(sce, assayname, k) {
  X <- reducedDim(sce, paste0("PCA_", assayname))
  knns <- get.knn(X, k = k)$nn.index
  out <- suppressWarnings(
    kBET(X, sce$Batch, k0 = k, knn = knns, do.pca = FALSE, plot = FALSE)
  )
  
  out$stats$kBET.observed
}

# Compute all metrics given k nearest neighbors
compute_metrics <- function(sce, assayname, k) {
  kbets <- data.frame(metric = "kBET", value = compute_kbet(sce, assayname, k))
  lisis <- data.frame(metric = "LISI", value = compute_lisi(
    reducedDim(sce, paste0("PCA_", assayname)), 
    meta_data = data.frame(Batch = sce$Batch),
    label_colnames = "Batch", perplexity = k/3
  )$Batch)
  KLs <- data.frame(metric = "KL", value = local_KL(sce, assayname, k))
  
  out <- rbind(kbets, lisis, KLs)
}

# Filter Highly Variable Genes 
filter_HVGs <- function(sce, percent_topHVGs = 0.10) {
  gene_var <- scran::modelGeneVar(sce, block = sce$Batch)
  topHVGs <- order(gene_var$p.value)[1:(round(nrow(gene_var) * percent_topHVGs))]
  
  sce[topHVGs, ]
}

# Perform seurat batch correction
# As Seurat works with a different class to represent the data, it first has
# to be converted to SingleCellExperiment.
# Computes PCA afterwards.
seurat_wrapper <- function(sce, d = 30) {
  sce_seur <- sce
  altExps(sce_seur) <- NULL  # remove MT genes in case separated into alternative experiment
  seur <- as.Seurat(sce_seur, counts = "counts", data = "logcounts")
  message("Performing Seurat batch correction ...")
  batch_list <- SplitObject(seur, split.by = "Batch")
  # Feature selection already performed with scran
  anchors <- FindIntegrationAnchors(
    batch_list, anchor.features = rownames(seur), dims = 1:d
    #verbose = FALSE
  )
  # Compute corrections
  combined <- IntegrateData(anchors)#, verbose = FALSE)
  # Put in sce object
  seurat_corr <- combined@assays$integrated@data[, colnames(sce)]  # reorder by sce_hvg colnames
  assay(sce, "seurat", withDimnames = FALSE) <- seurat_corr
  reducedDim(sce, "PCA_seurat") <- calculatePCA(seurat_corr, ncomponents = d, ntop = nrow(sce))  # calculate d principal components
  
  sce
}

# (Deprecated) Perform fastMNN batch correction
fastmnn_wrapper <- function(sce, k = 20, d = 50) {
  mnn <- fastMNN(
    sce, 
    batch = sce$Batch,
    d = d,
    k = k
  )
  reducedDim(sce, "PCA_mnn") <- reducedDim(mnn, "corrected")
  colnames(reducedDim(sce, "PCA_mnn")) <- paste0("PC", 1:d)
  assay(sce, "mnn") <- assay(mnn, "reconstructed")
  
  sce
}

# Perform fastMNN batch correction
# After batch correction, a (near-)full-rank expression matrix is reconstructed.
fastmnn_wrapper2 <- function(sce, k = 20, d = 50) {
  
  # Manually perform cosine normalization
  cosnorm <- apply(logcounts(sce), 2, function(x) x / sqrt(sum(x^2)))
  message(paste0("Performing fastMNN batch correction with k = ", k, ", d = ", d, " ..."))
  
  # Manually perform PCA
  PCA <- prcomp(t(cosnorm))  # Avoiding multiBatchNorm
  
  # Perform fastMNN using handmade PCs
  mnn <- fastMNN(
    t(PCA$x[, 1:d, drop = FALSE]),
    batch = sce$Batch,
    cos.norm = FALSE,  # already normalized
    d = NA,  # use the given PCA
    k = k
  )
  
  # Reconstruct expression matrix
  assay(sce, "mnn") <- BiocSingular::LowRankMatrix(
    rotation = PCA$rotation,
    components = cbind(reducedDim(mnn, "corrected"), PCA$x[, -(1:d)])
  )
  
  #reducedDim(sce, "PCA_mnn") <- calculatePCA(sce, exprs_values = "mnn", ntop = nrow(sce))
  
  # Perform PCA again on reconstructed matrix
  reducedDim(sce, "PCA_mnn") <- prcomp(t(assay(sce, "mnn")))$x
  
  sce
}

# Perform Harmony batch correction 
harmony_wrapper <- function(sce, d = 50) {
  
  # For when d = 1
  nclust <- NULL
  if (d == 1) nclust <- 1
  
  # Manually perform PCA
  PCA <- prcomp(t(assay(sce, "logcounts")))
  message("Performing Harmony batch correction with d = ", d, " ...")
  
  corrected <- HarmonyMatrix(
    data_mat = PCA$x[, 1:d, drop=FALSE], 
    meta_data = sce$Batch,
    do_pca = FALSE,
    nclust = nclust
    #verbose = FALSE
  )
  # reconstruct using loadings matrix
  assay(sce, "harmony") <- BiocSingular::LowRankMatrix(
    rotation = PCA$rotation,
    components = cbind(corrected, PCA$x[, -(1:d)])
  )
  
  # Perform PCA again on reconstructed matrix
  reducedDim(sce, "PCA_harmony") <- prcomp(t(assay(sce, "harmony")))$x
  
  sce
}


