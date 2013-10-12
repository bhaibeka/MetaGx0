########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

         
`duplicateFinder` <- 
function (eset, var.genes=1000, dupl.cor=0.95, nthread=1) {
  ## find duplicates based on correlation of gene expresison profiles
  #
  # Arga:
  #   eset: an expressionSet object
  #   var.genes: number of most variant genes used to define the expression profiles
  #
  # Returns
  #   list of duplicates sample names
  require(Biobase)
  if(nthread > 1) {
    require(parallel)
  }
  ## select the most variant genes
  ## at least in 80% of the datasets
  iix <- apply(exprs(eset), 1, function (x, y) {
    return ((sum(is.na(x)) / length(x)) < (1-y))
  }, y=0.8)
  varg <- Biobase::featureNames(eset)[iix][order(apply(exprs(eset)[iix, , drop=FALSE], 1, var, na.rm=TRUE), decreasing=TRUE)[1:var.genes]]
  
  splitix <- parallel::splitIndices(nx=length(sampleNames(eset)), ncl=nthread)
  mcres <- parallel::mclapply(splitix, function(splitix, data) {
      cores <- cor(x=data[ , splitix, drop=FALSE], y=data, use="complete.obs")
    }, data=exprs(eset)[varg, , drop=FALSE])
  cor.samples <- do.call(rbind, mcres)
  diag(cor.samples) <- NA
  ## create list of duplictaes for each sample
  duplix <- apply(cor.samples, 1, function (x, y) {
    res <- names(x)[!is.na(x) & x > y]
    return (res)
  }, y=dupl.cor)
  duplix <- duplix[sapply(duplix, length) > 0]
  return(duplix)
}

