# PIRF

Title: Phylogeny-Informed Random Forest

Version: 1.0

Date: 2025-09-18

Author: Hyunwook Koh

Maintainer: Hyunwook Koh <hyunwook.koh@stonybrook.edu>

Description: This R package implements Phylogeny-Informed Random Forest (PIRF), a method designed to improve predictive accuracy in human microbiome studies. PIRF adopts a localized approach: rather than treating all features as competing globally to be selected or weighted, it identifies informative features within each phylogenetic cluster - a localized group of evolutionarily and functionally related microbial features. This strategy enriches functional representations while reducing tree correlation. Specifically, PIRF partitions the microbial feature space into multiple phylogenetic clusters using phylogenetic tree information, computes feature importance scores within each cluster, and converts them into cluster-specific probabilities. Finally, these cluster-specific probabilities are integrated across all clusters to derive community-level selection probabilities which are used for feature selection and weighting in the Random Forest (Breiman, 2001) algorithm.

NeedsCompilation: no

Depends: R(>= 4.4.1), cluster, phyloseq, ranger

License: GPL v3.0

NeedsCompilation: no

## Reference

* Koh H. Phylogeny-Informed Random Forests for Human Microbiome Studies. (_In Review_)


## Troubleshooting Tips

If you have any problems for using this R package, please report in Issues (https://github.com/hk1785/PIRF/issues) or email Hyunwook Koh (hyunwook.koh@stonybrook.edu).

## Prerequites

cluster, ranger
```
install.packages(c("cluster", "ranger"))
```
phyloseq
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")
```

## Installation

```
library(devtools)
install_github("hk1785/PRIF", force=T)
```

---------------------------------------------------------------------------------------------------------------------------------------

## :mag: pirf

### Description
This function streamlines the entire implementation process for Phylogeny-Informed Random Forest (PIRF). Its core routines are built upon the ranger package (Wright and Ziegler, 2017), which is implemented in C++ and supports multi-core parallel computation.

### Syntax
```
pirf(X, y, phy.tree, num.trees = 1000, num.threads = 1, prop.ran.sel.features = c(1/10, "sqrt", "log"), oob.err = TRUE, ...)
```

### Arguments
* _X_ - A data frame containing feature abundances, where rows represent subjects (i = 1, ..., n) and columns represent features (j = 1, ... p).
* _y_ - A numeric vector containing the output values: class labels (e.g., 0: healthy, 1: diseased) for classification, or continuous response values for regression.
* _phy.tree_ - A rooted phylogenetic tree. The tip labels of this tree should match the column names of X.
* _num.trees_ - The number of bagged decision trees (Default: 1000).
* _num.threads_ - The number of threads for multi-core parallel computation (Default: 1).
* _prop.ran.sel.features_ - A vector containing the proportions of randomly selected features (Default: c(1/10, "sqrt", "log")). Here, "sqrt" refers to the square root of the number of features, and "log" refers to the base-2 logarithm of the number of features.
* _oob.err_ - A logical value indicating whether to compute out-of-bag prediction errors (TRUE) or not (FALSE). (Default: oob.err = TRUE).
* _..._ - Additional arguments passed from the ranger package (Wright and Ziegler, 2017).

### Values
A list containing multiple components, each corresponding to an element in prop.ran.sel.features. If prop.ran.sel.features contains only one element, the returned list will have a single component.

Each component is again a list with the following components: 
* _fit_ - The final fitted model
* _rev.features_ - Features eliminated (i.e., features with zero selection probabilities)
* _mtry_ - The number of randomly selected features used in the final model
* _sel.prob_ - Community-level selection probabilities
* _clust_ - Cluster labels for phylogenetic clusters
* _oob.err_ - Out-of-bag prediction errors (returned when oob.err = TRUE)

### Example (Classification)
Import requisite R packages
```
library(cluster)
library(phyloseq)
library(ranger)
```
Example Data: Gut microbiome data in the phyloseq format for the classification tasks of type 1 diabetes (Zhang et al., 2018).
```
data(t1d)

t1d

X <- as.data.frame(otu_table(t1d))
y <- as.numeric(unlist(sample_data(t1d)))
phy.tree <- phy_tree(t1d)
```
Perform PIRF.
```
out.cla <- pirf(X, y, phy.tree, num.trees = 1000, num.threads = 4, prop.ran.sel.features = c(1/10, "sqrt", "log"))
```
Check out-of-bag (OOB) prediction errors.
```
out.cla[["0.1"]]$oob.err
out.cla[["sqrt"]]$oob.err
out.cla[["log"]]$oob.err
```
Check community-level selection probabilities, cluster labels, and fitted model.
```
out.cla[["0.1"]]$sel.prob
out.cla[["0.1"]]$clust
out.cla[["0.1"]]$fit
```
Compute predicted responses.
```
predict(out.cla[["0.1"]]$fit, data = X)$predictions



#######################
# Ordination Plotting #
#######################

##### First kernel (Jaccard) 

# Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"

cols <- c("green2", "blue2", "purple2", "red2")

# Singular value decomposition

svd.out <- svd(out$Ks[[1]]) # First kernel (Jaccard)
pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
eivals <- svd.out$d[svd.out$d > 0]

# 2D plot

plot(pcs[,1], pcs[,2], 
     xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     main = "Jaccard", 
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)

# 3D plot

plot3d( 
  x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
  type = "s", 
  zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
  xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
  ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
  cex.axis = 0.5)

rglwidget()

### Second kernel (Bray-Curtis) 

# Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"

cols <- c("green2", "blue2", "purple2", "red2")

# Singular value decomposition

svd.out <- svd(out$Ks[[2]]) # Second kernel (Bray-Curtis)
pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
eivals <- svd.out$d[svd.out$d > 0]

# 2D plot

plot(pcs[,1], pcs[,2], 
     xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     main = "Bray-Curtis", 
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)

# 3D plot

plot3d( 
  x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
  type = "s", 
  zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
  xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
  ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
  cex.axis = 0.5)

rglwidget()

### Third kernel (UUniFrac) 

# Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"

cols <- c("green2", "blue2", "purple2", "red2")

# Singular value decomposition

svd.out <- svd(out$Ks[[3]]) # Third kernel (UUniFrac)
pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
eivals <- svd.out$d[svd.out$d > 0]

# 2D plot

plot(pcs[,1], pcs[,2], 
     xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     main = "UUniFrac", 
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)

# 3D plot

plot3d( 
  x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
  type = "s", 
  zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
  xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
  ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
  cex.axis = 0.5)

rglwidget()

### Fourth kernel (GUniFrac 0.25) 

# Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"

cols <- c("green2", "blue2", "purple2", "red2")

# Singular value decomposition

svd.out <- svd(out$Ks[[4]]) # Fourth kernel (GUniFrac 0.25)
pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
eivals <- svd.out$d[svd.out$d > 0]

# 2D plot

plot(pcs[,1], pcs[,2], 
     xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     main = "GUniFrac (0.25)", 
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)

# 3D plot

plot3d( 
  x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
  type = "s", 
  zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
  xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
  ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
  cex.axis = 0.5)

rglwidget()

### Fifth kernel (GUniFrac 0.5) 

# Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"

cols <- c("green2", "blue2", "purple2", "red2")

# Singular value decomposition

svd.out <- svd(out$Ks[[5]]) # Fifth kernel (GUniFrac 0.5)
pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
eivals <- svd.out$d[svd.out$d > 0]

# 2D plot

plot(pcs[,1], pcs[,2], 
     xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     main = "GUniFrac (0.5)", 
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)

# 3D plot

plot3d( 
  x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
  type = "s", 
  zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
  xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
  ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
  cex.axis = 0.5)

rglwidget()

### Sixth kernel (GUniFrac 0.75) 

# Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"

cols <- c("green2", "blue2", "purple2", "red2")

# Singular value decomposition

svd.out <- svd(out$Ks[[6]]) # Sixth kernel (GUniFrac 0.75)
pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
eivals <- svd.out$d[svd.out$d > 0]

# 2D plot

plot(pcs[,1], pcs[,2], 
     xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     main = "GUniFrac (0.75)", 
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)

# 3D plot

plot3d( 
  x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
  type = "s", 
  zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
  xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
  ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
  cex.axis = 0.5)

rglwidget()

### Seventh kernel (WUniFrac) 

# Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"

cols <- c("green2", "blue2", "purple2", "red2")

# Singular value decomposition

svd.out <- svd(out$Ks[[7]]) # Seventh kernel (WUniFrac)
pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
eivals <- svd.out$d[svd.out$d > 0]

# 2D plot

plot(pcs[,1], pcs[,2], 
     xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
     main = "WUniFrac", 
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)

# 3D plot

plot3d( 
  x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
  type = "s", 
  zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
  xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
  ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
  cex.axis = 0.5)

rglwidget()
```
