pirf <-
function(X, y, phy.tree, num.trees = 1000, num.threads = 1, prop.ran.sel.features = c(1/10, "sqrt", "log"), oob.err = TRUE, ...) {
  
  out <- list()
  
  for (i in 1:length(prop.ran.sel.features)) {
    prsf <- prop.ran.sel.features[i]
    if (!is.na(as.numeric(prsf))) {
      if (prsf <= 0 | prsf > 1) {
        stop("the proportion of randomly selected features must be a number between 0 and 1")
      } else {
        nrsf <- ceiling(ncol(X)*as.numeric(prsf))
      }
    }
    if (prsf == "sqrt") {
      nrsf <- ceiling(sqrt(ncol(X)))
    }
    if (prsf == "log") {
      nrsf <- ceiling(log2(ncol(X)))
    }
    
    norm.tree <- phy.tree
    br.len <- norm.tree$edge.length
    norm.tree$edge.length <- br.len/sum(br.len)
    cop <- ape::cophenetic.phylo(norm.tree)
    
    pam.fit.list <- list()
    sil.widths <- numeric()
    for (k in 2:10) {
      pam.fit.list[[k-1]] <- pam.fit <- pam(as.dist(cop), k = k, diss = TRUE)
      sil.widths[k-1] <- pam.fit$silinfo$avg.width
    }
    
    ind.min <- which.min(sil.widths)
    clusts <- unique(pam.fit.list[[ind.min]]$clustering)
    X.clust.list <- list()
    fis.clust.list <- list()
    non.neg.fis.clust.list <- list()
    cum.fis.clust.list <- list()
    rev.features.list <- list()
    for (j in 1:length(clusts)) {
      ind.clust <- which(pam.fit.list[[ind.min]]$clustering == j)
      if (length(ind.clust) == 1) {
        X.clust <- as.matrix(X[,ind.clust])
        colnames(X.clust) <- names(ind.clust)
        X.clust.list[[j]] <- X.clust 
        fis.clust <- 1
        names(fis.clust) <- names(ind.clust)
        fis.clust.list[[j]] <- fis.clust
        non.neg.fis.clust.list[[j]] <- fis.clust
        cum.fis.clust.list[[j]] <- fis.clust
        rev.features.list[[j]] <- rev.features <- ""
      } else {
        X.clust.list[[j]] <- X.clust <- X[,ind.clust]
        dat.clust <- as.data.frame(cbind(y, X.clust))
        imp.fit.clust.list <- list()
        for (k in 1:length(prop.ran.sel.features)) {
          prsf.clust <- prop.ran.sel.features[k] 
          if (!is.na(as.numeric(prsf.clust))) {
            nrsf.clust <- ceiling(ncol(X.clust)*as.numeric(prsf.clust))
          }
          if (prsf.clust == "sqrt") {
            nrsf.clust <- ceiling(sqrt(ncol(X.clust)))
          }
          if (prsf.clust == "log") {
            nrsf.clust <- ceiling(log2(ncol(X.clust)))
          }
          imp.fit.clust.list[[k]] <- ranger(y ~ ., data = dat.clust, num.trees = num.trees, mtry = nrsf.clust, importance = "permutation", num.threads = num.threads, ...)$variable.importance
        }
        if (length(prop.ran.sel.features) == 1) {
          fis.clust.list[[j]] <- fis.clust <- unlist(imp.fit.clust.list)
        }
        if (length(prop.ran.sel.features) > 1) {
          fis.clust.list[[j]] <- fis.clust <- rowMeans(as.data.frame(imp.fit.clust.list))
        }
        ind.rev <- which(fis.clust < 0)
        if (length(ind.rev) > 0) {
          fis.clust[ind.rev] <- 0
          non.neg.fis.clust.list[[j]] <- non.neg.fis.clust <- fis.clust
          cum.fis.clust.list[[j]] <- cum.fis.clust <- (non.neg.fis.clust/sum(non.neg.fis.clust))*length(fis.clust)
          rev.features.list[[j]] <- rev.features <- names(ind.rev)
        } else {
          non.neg.fis.clust.list[[j]] <- non.neg.fis.clust <- fis.clust
          cum.fis.clust.list[[j]] <- cum.fis.clust <- (non.neg.fis.clust/sum(non.neg.fis.clust))*length(fis.clust)
          rev.features.list[[j]] <- rev.features <- ""
        }
      }
    }
    
    X.clusts <- X.clust.list[[1]]
    for (l in 2:length(clusts)) {
      X.clusts <- cbind(X.clusts, X.clust.list[[l]])
    }
    non.neg.fis.clust <- unlist(non.neg.fis.clust.list)
    cum.fis.clust <- unlist(cum.fis.clust.list)/sum(unlist(cum.fis.clust.list))
    
    dat.clusts <- as.data.frame(cbind(y, X.clusts))
    
    if (oob.err) {
      fit <- ranger(y ~ ., data = dat.clusts, num.trees = num.trees, mtry = ceiling((ncol(X) - sum(unlist(rev.features.list) != ""))*(nrsf/ncol(X))), 
                    split.select.weights = cum.fis.clust, oob.error = TRUE, num.threads = num.threads, ...) 
      
      oob.err <- fit$prediction.error
      
      fit <- ranger(y ~ ., data = dat.clusts, num.trees = num.trees, mtry = ceiling((ncol(X) - sum(unlist(rev.features.list) != ""))*(nrsf/ncol(X))), 
                    split.select.weights = cum.fis.clust, num.threads = num.threads, ...) 
      
      out[[i]] <- list(fit = fit, rev.features = unlist(rev.features.list), mtry = nrsf, sel.prob = cum.fis.clust, clust = pam.fit.list[[ind.min]]$clustering, oob.err = oob.err)
    } else {
      fit <- ranger(y ~ ., data = dat.clusts, num.trees = num.trees, mtry = ceiling((ncol(X) - sum(unlist(rev.features.list) != ""))*(nrsf/ncol(X))), 
                    split.select.weights = cum.fis.clust, num.threads = num.threads, ...) 
      
      out[[i]] <- list(fit = fit, rev.features = unlist(rev.features.list), mtry = nrsf, sel.prob = cum.fis.clust, clust = pam.fit.list[[ind.min]]$clustering)
    }
    
  }
  
  names(out) <- prop.ran.sel.features
  
  return(out)
  
}
