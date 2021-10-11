# This file contains the code for generating the histogram in Figure 2.1.

source("helpers.R")  # some handmade functions
library(tidyverse)  # for pipe operator and ggplot

# Specify dimensions
N <- 1000
M <- 2000

# Generate noise matrix
set.seed(2)
X <- matrix(rnorm(N*M), nrow = M, ncol = N)

# Create and add rank-1 perturbation
u <- rep(1, M)
v <- rep(c(0, 0.1), each = N/2)
P <- u %*% t(v)
X_tilde <- X + P

SVD <- svd(X_tilde)
singvals <- SVD$d / sqrt(max(N, M)-1)
bins <- 100
fill <- "#777777"

# Get MP bounds
c <- min(N/M, M/N)
a <- (1 - sqrt(c))^2 # for the min boundary
b <- (1 + sqrt(c))^2  # for the max boundary

p <- 
  data.frame(singvals) %>%
  ggplot(aes(x = singvals)) +
  geom_histogram(
    aes(y = ..density..), 
    fill = fill, bins = bins
  ) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_function(
    fun = dmp, 
    args = list(N = N, P = M),
    n = 1001,
    xlim = sqrt(c(a, b)),
    color = 'black',
    size = 0.8
  ) +
  xlim(c(0, 3.2)) +
  ylim(c(0, 1.2)) +
  xlab(bquote(sigma / sqrt(P - 1)))

ggsave(
  "./plots_report/MP_sketch.pdf", 
  p, 
  width = 120, height = 80, 
  units = "mm", 
  device = "pdf"
)
