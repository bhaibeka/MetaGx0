########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

`platformMerging` <- 
function (esets) {
  #Open every GSE, specify in the configuration CSV, from InsilicoDB and create a list of Eset
  #structure call gselist
  #
  # Args:
  #     esets: list of ExpressionSet objects
  #   
  # Returns:     
  #     The eset with all platforms merged
  
  if (class(esets) == "ExpressionSet") {
    ## a single eset, no platform to merge
    return(esets)
  }
  if (is.list(esets)) {
    esets.check <- sapply(esets, function(x) { return(class(x) == "ExpressionSet") })
    if (any(!esets.check)) { stop("Some esets in the list are not ExpressionSet") }
    if (length(esets) == 1) {
      eset <- esets[[1]]
    } else {
      for (j in 2:length(esets)) {
        ## merge expression data
        Biobase::exprs(eset) <- rbind(Biobase::exprs(eset), Biobase::exprs(esets[[j]]))
        ## merge probe annotations
        Biobase::featureData(eset)@data <- rbind(Biobase::featureData(eset)@data, Biobase::featureData(esets[[j]])@data)
      }
      ## remove duplicated Entrez Gene IDs
      duplix <- duplicated(as.character(Biobase::featureData(eset)@data[ , "ENTREZID"]))
      Biobase::exprs(eset) <- Biobase::exprs(eset)[!duplix, , drop=FALSE]
    }
    return(eset)
  }
}

