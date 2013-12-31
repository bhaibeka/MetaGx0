########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################


`subtypeCorrelation` <- 
function (eset, geneid, plot=TRUE, method=c("pearson", "spearman"), weighted=TRUE, condensed=TRUE, resdir="cache", nthread=1) {
  ## assess (weighted) correlation between gene expression with respect to subtypes
  #
  # Arga:
  #   eset: an expressionSet object
  #   geneid: vector of Entrez Gene IDs. If missing, all genes will be considered.
  #   plot: should the correlation heatmap be plotted?
  #   method: method for correlation
  #   resdir
  #   nthread:
  #
  # Returns
  #   list containing p-values for comparisons
  #   Kruskal-Wallist to test whether the expression of the genes(s) of interest is dependent on the molecular subtypes
  #   pairwise wilcoxon rank sum test p-values, is the expression of the gene(s) of interest higher in the subtype in rows compared to the subtype in column?
  
  ######################
  
  wcor <- function (d, w, na.rm=TRUE) {
    s <- sum(w, na.rm=na.rm)
    m1 <- sum(d[ , 1L] * w, na.rm=na.rm) / s
    m2 <- sum(d[ , 2L] * w, na.rm=na.rm) / s
    res <- (sum(d[ , 1L] * d[ , 2L] * w, na.rm=na.rm) / s - m1 * m2) / sqrt((sum(d[ , 1L]^2 * w, na.rm=na.rm) / s - m1^2) * (sum(d[ , 2L]^2 * w, na.rm=na.rm) / s - m2^2))
    return (res)
  }
  
  ######################
  
  if (class(eset) != "ExpressionSet") {
    stop("Handling list of expressionSet objects is not implemented yet")
  }
  
  if (missing(geneid)) {
    gened <- Biobase::fData(eset)[ , "ENTREZID"]
  }
  
  if (!file.exists(file.path(resdir))) { dir.create(file.path(resdir), showWarnings=FALSE, recursive=TRUE) }
  
  ## for a single expressionSet object
  
  ## extract subtypes
  sbts <- Biobase::pData(eset)[ , "subtype"]
  names(sbts) <- rownames(Biobase::pData(eset))
  sbtu <- levels(sbts)
  if (sum(table(sbts) > 3) < 2) {
    warning("Not enough tumors in each subtype")
    return(NULL)
  }
  sbts.proba <- Biobase::pData(eset)[ , sprintf("subtyproba.%s", sbtu), drop=FALSE]
  if (!weighted) {
    sbts.proba <- t(apply(sbts.proba, 1, function (x) {
      xx <- array(0, dim=length(x))
      xx[which.max(x)] <- 1
      return (xx)
    }))
  }
  sbts.proba <- cbind(1, sbts.proba)
  colnames(sbts.proba) <- c("ALL", sbtu)
  
  ## extract genes
  gid <- paste("geneid", intersect(geneid, Biobase::fData(eset)[ , "ENTREZID"]), sep=".")
  gsymb <- Biobase::fData(eset)[gid, "SYMBOL"]
  names(gsymb) <- gid
  if (length(gid) == 0) {
    stop("Genes not in the expressionSet object")
  }
  if (length(gid) < length(geneid)) {
    warning(sprintf("%i/%i genes were present in the expressionSet object", length(gid), length(geneid)))
  }
  expr <- Biobase::exprs(eset)[gid, , drop=FALSE]
  if (method == "spearman") {
    expr <- t(apply(expr, 1, rank))
  }
  
  ## compute subtype-specific pairwise correlation across query genes
  pairs <- t(combn(1:length(gid), 2, simplify=TRUE))
  splitix <- parallel::splitIndices(nx=nrow(pairs), ncl=nthread)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, ...) {    
    res <- apply(pairs[x, , drop=FALSE], 1, function (x, ...) {
      res <- apply(sbts.proba, 2, function (w, x, expr) {
        return (wcor(d=t(expr[x, , drop=FALSE]), w=w))
      }, x=x, expr=expr)
      return (res)
    })
    return (res)
  }, expr=expr, sbts.proba=sbts.proba)
  res <- t(do.call(cbind, mcres))
  
  rr <- unlist(apply(res, 2, function (x, y, gid) {
    rr <- matrix(NA, nrow=length(gid), ncol=length(gid), dimnames=list(gid, gid))
    rr[y] <- x
    rr[y[ , 2:1]] <- x
    diag(rr) <- 1
    
    return(list(rr))
  }, y=pairs, gid), recursive=FALSE)
  
  dd <- lapply(rr, data.frame)
  if (condensed) {
    WriteXLS::WriteXLS(x="dd", ExcelFileName=file.path(resdir, sprintf("subtype_correlation_%s.xls", method)), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
  } else {
    mapply(function(x, y, method, resdir) {
       WriteXLS::WriteXLS("x", ExcelFileName=file.path(resdir, sprintf("subtype_association_%s_%s.xls", method, y)), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
     }, x=dd, y=names(dd), method=method, resdir=resdir)
  }
  return(rr)
}


