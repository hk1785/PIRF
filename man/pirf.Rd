\name{pirf}
\alias{pirf}
\title{
Phylogeny-Informed Random Forest (PIRF)
}
\description{
The pirf function streamlines the entire implementation process for Phylogeny-Informed Random Forest (PIRF). Its core routines are built upon the ranger package (Wright and Ziegler, 2017), which is implemented in C++ and supports multi-core parallel computation.
}
\usage{
pirf(X, y, phy.tree, num.trees = 1000, num.threads = 1, prop.ran.sel.features = c(1/10, "sqrt", "log"), oob.err = TRUE, ...)
}
\arguments{
  \item{X}{
A data frame containing feature abundances, where rows represent subjects (i = 1, ..., n) and columns represent features (j = 1, ... p).
}
  \item{y}{
A numeric vector containing the output values: class labels (e.g., 0: healthy, 1: diseased) for classification, or continuous response values for regression.
}
  \item{phy.tree}{
A rooted phylogenetic tree. The tip labels of this tree should match the column names of X.
}
  \item{num.trees}{
The number of bagged decision trees (Default: 1000).
}
  \item{num.threads}{
The number of threads for multi-core parallel computation (Default: 1).
}
  \item{prop.ran.sel.features}{
A vector containing the proportions of randomly selected features (Default: c(1/10, "sqrt", "log")). Here, "sqrt" refers to the square root of the number of features, and "log" refers to the base-2 logarithm of the number of features.
}
  \item{oob.err}{
A logical value indicating whether to compute out-of-bag prediction errors (TRUE) or not (FALSE). (Default: oob.err = TRUE).
}
  \item{\dots}{
Additional arguments passed from the ranger package (Wright and Ziegler, 2017).
}
}

\value{
A list containing multiple components, each corresponding to an element in prop.ran.sel.features. If prop.ran.sel.features contains only one element, the returned list will have a single component.

Each component is again a list with the following components: 
fit (the final fitted model); 
rev.features (features eliminated; that is, features with zero selection probabilities);
mtry (the number of randomly selected features used in the final model);
sel.prob (community-level selection probabilities);
clust (cluster labels for phylogenetic clusters);
oob.err (out-of-bag prediction errors (returned when oob.err = TRUE)).
}
}
\references{
Koh H. Phylogeny-Informed Random Forests for Human Microbiome Studies (In Review)
}
\author{
Hyunwook Koh
}

\examples{
## Classification

data(t1d)

t1d

X <- as.data.frame(otu_table(t1d))
y <- as.numeric(unlist(sample_data(t1d)))
phy.tree <- phy_tree(t1d)

out.cla <- pirf(X, y, phy.tree, num.trees = 1000, 
num.threads = 4, prop.ran.sel.features = c(1/10, "sqrt", "log"))
  
out.cla[["0.1"]]$oob.err
out.cla[["sqrt"]]$oob.err
out.cla[["log"]]$oob.err

out.cla[["0.1"]]$sel.prob
out.cla[["0.1"]]$clust

out.cla[["0.1"]]$fit

predict(out.cla[["0.1"]]$fit, data = X)$predictions

## Regression

data(age.oral)

age.oral

X <- as.data.frame(otu_table(age.oral))
y <- as.numeric(unlist(sample_data(age.oral)))
phy.tree <- phy_tree(age.oral)

out.reg <- pirf(X, y, phy.tree, num.trees = 1000, 
num.threads = 4, prop.ran.sel.features = c(1/10, "sqrt", "log"))

sqrt(out.reg[["0.1"]]$oob.err)
sqrt(out.reg[["sqrt"]]$oob.err)
sqrt(out.reg[["log"]]$oob.err)

out.reg[["sqrt"]]$sel.prob
out.reg[["sqrt"]]$clust

out.reg[["sqrt"]]$fit

predict(out.reg[["sqrt"]]$fit, data = X)$predictions
}


