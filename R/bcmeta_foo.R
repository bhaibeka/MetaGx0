########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

checkNames <- function (eset) {
  ## check feature and sample names in esets
  #
  # Args:
  #   esets: expressionSet object
  #
  # Returns:
  #   esets: updated expressionSet object with atching feature and sample names
  
  require(Biobase)
  ## check feature names
  check.feature <- intersect(rownames(exprs(eset)), rownames(fData(eset)))
  if (length(check.feature) == 0) {
    warning("Names of features do not match between expressions and annotations")
    return (NULL)
  } else {
    if (length(check.feature) != nrow(exprs(eset)) || length(check.feature) != nrow(fData(eset))) {
      warning("Some features are missing between expressions and annotations")
    }
  }
  ## check sample names
  check.sample <- intersect(colnames(exprs(eset)), rownames(pData(eset)))
  if (length(check.sample) == 0) {
    warning("Names of samples do not match between expressions and phenotypes")
    return (NULL)
  } else {
    if (length(check.sample) != ncol(exprs(eset)) || length(check.sample) != nrow(pData(eset))) {
      warning("Some samples are missing between expressions and phenotypes")
    }
  }
  exprs(eset) <- exprs(eset)[check.feature, check.sample, drop=FALSE]
  fData(eset) <- fData(eset)[check.feature, , drop=FALSE]
  pData(eset) <- pData(eset)[check.sample, , drop=FALSE]
  return (eset)
}


probeGeneMapping <- function (eset, platform=c("MISC", "GPL8300", "GPL96", "GPL97", "GPL570", "GPL1352"), method=c("variance", "jetset")){
  ## probe-gene mapping: which package to use for which platform
  #
  # Args:
  #    eset: ExpressionSet object
  #    platform: identifier of the microarray platform (usually a GPL id from GEO)
  #    method: either the most variant probes or jetset for Affymetrix platform
  #   
  # Returns:     
  #     updated ExpressionSet object with single probe per Entrez gene id
  
  require(org.Hs.eg.db)
  
  platform <- match.arg(platform)
  method <- match.arg(method)
  
  platf.map <- rbind(c("MISC", "variance", ""),
    c("GPL8300", "jetset", "hgu95av2"),
    c("GPL96", "jetset", "hgu133a"),
    c("GPL97", "jetset", "hgu133b"),
    c("GPL570", "jetset", "hgu133plus2"),
    c("GPL1352", "jetset", "u133x3p"))
  dimnames(platf.map) <- list(platf.map[ , 1], c("platform", "method", "parameters"))
  if (!is.element(method, platf.map[platf.map[ , "platform"], "method"])) {
    stop(sprintf("Method %s cannot be applied on platform %s\nUse the following method(s) instead: %s", method, platform, paste(x=platf.map[platf.map[ , "platform"] == platform, "method"], collapse=", ")))
  }
  params <- platf.map[which(platform == platf.map[ , "platform"]), "parameters"]
  
  ## keep only ENTREZID and SYMBOL in feature annotation
  Biobase::fData(eset) <- Biobase::fData(eset)[ , c("ENTREZID", "SYMBOL"), drop=FALSE]
  switch(method,
    "jetset" = {
      require(jetset.bhk)
      js <- jetset.bhk::jscores(chip=params, probeset=rownames(Biobase::exprs(eset)))
      js <- js[rownames(Biobase::exprs(eset)), , drop=FALSE]
      ## identify the best probeset for each Entrez Gene ID
      geneid1 <- as.character(js[ ,"EntrezID"])
      names(geneid1) <- rownames(js)
      geneid2 <- sort(unique(geneid1))
      names(geneid2) <- paste("geneid", geneid2, sep=".")
      gix1 <- !is.na(geneid1)
      gix2 <- !is.na(geneid2)
      geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
      ## probes corresponding to common gene ids
      gg <- names(geneid1)[is.element(geneid1, geneid.common)]
      gid <- geneid1[is.element(geneid1, geneid.common)]
      ## duplicated gene ids
      gid.dupl <- unique(gid[duplicated(gid)])
      gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
      ## unique gene ids
      gid.uniq <- gid[!is.element(gid, gid.dupl)]
      gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
      ## which are the best probe for each gene
      js <- data.frame(js, "best"=FALSE)
      js[gg.uniq, "best"] <- TRUE
      ## data for duplicated gene ids
      if(length(gid.dupl) > 0) {	
      	## use jetset oevrall score to select the best probesets
      	myscore <- js[gg.dupl,"overall"]
      	myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
      	myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
      	myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
      	js[myscore[ ,"probe"], "best"] <- TRUE
      }
      ## update the esets
      probes <- rownames(Biobase::exprs(eset))[js[ , "best"]]
      names(probes) <- paste("geneid", js[js[ , "best"], "EntrezID"], sep=".")
      gid <- js[js[ , "best"], "EntrezID"]
      gsymb <- js[js[ , "best"], "symbol"]
      Biobase::exprs(eset) <- Biobase::exprs(eset)[probes, , drop=FALSE]
      rownames(Biobase::exprs(eset)) <- names(probes)
      Biobase::featureData(eset)@data <- Biobase::featureData(eset)@data[probes, , drop=FALSE]
      rownames(Biobase::featureData(eset)@data) <- names(probes)
      Biobase::featureData(eset)@data[ , "ENTREZID"] <- gid
      Biobase::featureData(eset)@data[ , "SYMBOL"] <- gsymb
    },
    "variance" = {
      ## other platform, select the most variant probe per Entrez Gene ID
      gid <- as.character(Biobase::featureData(eset)@data[ , "ENTREZID"])
      names(gid) <- rownames(Biobase::exprs(eset))
      ugid <- sort(unique(gid))
      rr <- genefu::geneid.map(geneid1=gid, data1=t(Biobase::exprs(eset)), geneid2=ugid)
      probes <- colnames(rr$data1)
      names(probes) <- paste("geneid", rr$geneid1, sep=".")
      Biobase::exprs(eset) <- Biobase::exprs(eset)[probes, , drop=FALSE]
      rownames(Biobase::exprs(eset)) <- names(probes)
      Biobase::featureData(eset)@data <- Biobase::featureData(eset)@data[probes, , drop=FALSE]
      rownames(Biobase::featureData(eset)@data) <- names(probes)
      ## get the gene symbols from entrez gene id using org.Hs.eg.db
      gs <- toTable(org.Hs.egSYMBOL)
      gs <- gs[!duplicated(gs[ , "gene_id"]), , drop=FALSE]
      rownames(gs) <- gs[ , "gene_id"]
      gs <- gs[as.character(Biobase::featureData(eset)@data[ , "ENTREZID"]), "symbol"]
      Biobase::featureData(eset)@data[ , "SYMBOL"] <- as.character(gs)
    },
    {
      stop(sprintf("Unknow method for probe-gene mapping for platform %s", platform))
    }
  )
  return(eset)
}


platformMerging <- function (esets){
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


subtypeClassification <- function (eset, model=c("scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003")) {
  # Classify GSE subtype and give out the probability of being that subtype
  #
  # Args:
  #   eset: expression set to classify
  #
  # Returns:
  #   eset: expression set with updated subtype information (subtype claling and subtype.probabilities)
  
  require(genefu)
  require(Biobase)
  model <- match.arg(model)
  
  ## convert SCM to SSP nomenclature
  sbt.conv <- rbind(c("ER+/HER2- Low Prolif", "LumA"),
    c("ER+/HER2- High Prolif", "LumB"),
    c("HER2+", "Her2"),
    c("ER-/HER2-", "Basal"))
  colnames(sbt.conv) <- c("SCM.nomenclature", "SSP.nomenclature")
    
  sbtn2 <- c("LumA", "LumB", "Her2", "Basal", "Normal")
  
  datage <- t(Biobase::exprs(eset))   
  annotge <- cbind("probe"=rownames(Biobase::featureData(eset)@data), "EntrezGene.ID"=as.character(Biobase::featureData(eset)@data[ , "ENTREZID"]))
  rownames(annotge) <- as.character(annotge[ , "probe"])
  
  switch(model,
    "scmgene" = {
      sbts <- genefu::subtype.cluster.predict(sbt.model=scmgene.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype2", "subtype.proba2")]
      names(sbts) <- c("subtype", "subtype.proba")
      ## use SSP nomenclature
      ## update subtype calling
      ss <- factor(x=sbts$subtype)
      levels(ss)[match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], levels(ss))] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
      sbts$subtype <- as.character(ss)
      ## update subtype probabilities
      iix <- match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], colnames(sbts$subtype.proba))
      colnames(sbts$subtype.proba)[iix] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
    },
    "scmod1" = {
      sbts <- genefu::subtype.cluster.predict(sbt.model=scmod1.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype2", "subtype.proba2")]
      names(sbts) <- c("subtype", "subtype.proba")
      ## use SSP nomenclature
      ## update subtype calling
      ss <- factor(x=sbts$subtype)
      levels(ss)[match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], levels(ss))] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
      sbts$subtype <- as.character(ss)
      ## update subtype probabilities
      iix <- match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], colnames(sbts$subtype.proba))
      colnames(sbts$subtype.proba)[iix] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
    },
    "scmod2" = {
      sbts <- genefu::subtype.cluster.predict(sbt.model=scmod2.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype2", "subtype.proba2")]
      names(sbts) <- c("subtype", "subtype.proba")
      ## use SSP nomenclature
      ## update subtype calling
      ss <- factor(x=sbts$subtype)
      levels(ss)[match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], levels(ss))] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
      sbts$subtype <- as.character(ss)
      ## update subtype probabilities
      iix <- match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], colnames(sbts$subtype.proba))
      colnames(sbts$subtype.proba)[iix] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
    },
    "pam50" = {
      sbts <- genefu::intrinsic.cluster.predict(sbt.model=pam50.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype", "subtype.proba")]
    },
    "ssp2006" = {
      sbts <- genefu::intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype", "subtype.proba")]
    },
    "ssp2003" = {
      sbts <- genefu::intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype", "subtype.proba")]
    },
    {
      stop("Unknown subtype classification model")
    })
  
  ## merge clinical information and subtype classification
  colnames(sbts$subtype.proba) <- paste("subtyproba", colnames(sbts$subtype.proba), sep=".")
  Biobase::phenoData(eset)@data <- cbind(Biobase::phenoData(eset)@data, "subtype"=sbts$subtype, sbts$subtype.proba)
  
  return(eset)
}


datasetMerging <- function (esets, method=c("union", "intersect"), standardization=c("quantile", "robust.scaling", "scaling", "none"), nthread=1) {
  ## function merging all individual esets and merging them into a big eset
  #
  # Args:
  #   gselist: The list containing all GSE file that need to be merged.
  #   STL: List containing all the subtype information for each GSE dataset.        
  #   duplication.checker: A marker, either TRUE or FALSE if you you want to verify
  #                        wheter or not you have duplicate samples into your master 
  #                        gene expression matrix.
  #   survdata: For t.fs and e.fs
  #   time.cens: maximum follow up (years)
  #   Method: either "unique" or "intersect" is use to for selecting geneid
  #
  # Returns:     
  #     The merging eset
  
  require(Biobase)
  require(genefu)
  if(nthread > 1) {
    require(parallel)
  }
  
  method <- match.arg(method)
  standardization <- match.arg(standardization)
  
  ## all unique Entrez gene ids
  ## gene ids
  ugid <- lapply(esets, function(x) { return(Biobase::featureData(x)@data) })
  ugid <- do.call(rbind, ugid)
  ugid <- ugid[!is.na(ugid[ , "ENTREZID"]) & !duplicated(as.character(ugid[ , "ENTREZID"])), , drop=FALSE]
  rownames(ugid) <- gsub(sprintf("(%s).", paste(names(esets), collapse="|")), "", rownames(ugid))
  switch (method,
    "union" = {
      feature.merged <- ugid
    },
    "intersect" = {
      feature.merged <- lapply(esets, function(x) { return(as.character(Biobase::featureData(x)@data[ , "ENTREZID"])) })
      feature.merged <- table(unlist(feature.merged))
      feature.merged <- names(feature.merged)[feature.merged == length(esets)]
      feature.merged <- ugid[match(feature.merged, as.character(ugid[ , "ENTREZID"])), , drop=FALSE]
    },
    {
      stop("Unknown method")
    }
  )
  ## expression data
  exprs.merged <- lapply(esets, function (x, y) {
    ee <- Biobase::exprs(x)
    eem <- matrix(NA, nrow=length(y), ncol=ncol(ee), dimnames=list(y, colnames(ee)))
    eem[rownames(ee), colnames(ee)] <- ee
    return (eem)
  }, y=rownames(feature.merged))
  exprs.merged <- do.call(cbind, exprs.merged)
  ## clinical info
  ucid <- lapply(esets, function(x) { return(colnames(phenoData(x)@data)) })
  ucid <- table(unlist(ucid))
  ucid <- names(ucid)[ucid == length(esets)]
  clinicinfo.merged <- lapply(esets, function (x , y) {
    ee <- Biobase::pData(x)[ , y, drop=FALSE]
  }, y=ucid)
  clinicinfo.merged <- do.call(rbind, clinicinfo.merged)
  rownames(clinicinfo.merged) <- gsub(sprintf("(%s).", paste(names(esets), collapse="|")), "", rownames(clinicinfo.merged))
  ## create a merged expressionSet object
  eset.merged <- ExpressionSet(assayData=exprs.merged, phenoData=AnnotatedDataFrame(data=clinicinfo.merged), featureData=AnnotatedDataFrame(data=feature.merged))
  experimentData(eset.merged)@preprocessing <- list("normalization"="mixed", package="unspecified", version="0")
  annotation(eset.merged) <- "mixed"
  
  ## standardization
  switch(standardization,
    "none" = {
      ## do nothing
    },
    "quantile" = {
      require(limma)
      require(genefu)
      ## robust scaling followed by quantile normalization
      ee <- exprs(eset.merged)
      # ee <- apply(ee, 2, genefu::rescale)
      splitix <- parallel::splitIndices(nx=ncol(ee), ncl=nthread)
      mcres <- parallel::mclapply(splitix, function(x, data) {
        res <- apply(data[ , x, drop=FALSE], 2, function (dx) {
          return ((genefu::rescale(dx, q=0.05, na.rm=TRUE) - 0.5) * 2)
        })
        return (res)
      }, data=ee, mc.cores=nthread)
      ee <- do.call(cbind, mcres)
      ## quantile normalization
      ee <- limma::normalizeBetweenArrays(object=ee, method="quantile")
      exprs(eset.merged) <- ee
    },
    "robust.scling" = {
      ## robust scaling
      require(genefu)
      ## robust scaling
      ee <- exprs(eset.merged)
      # ee <- apply(ee, 2, genefu::rescale)
      splitix <- parallel::splitIndices(nx=ncol(ee), ncl=nthread)
      mcres <- parallel::mclapply(splitix, function(x, data) {
        res <- apply(data[ , x, drop=FALSE], 2, function (dx) {
          return ((genefu::rescale(dx, q=0.05, na.rm=TRUE) - 0.5) * 2)
        })
        return (res)
      }, data=ee, mc.cores=nthread)
      ee <- do.call(cbind, mcres)
      exprs(eset.merged) <- ee
    },
    "scaling" = {
      ## traditional scaling
      # robust scaling
      ee <- exprs(eset.merged)
      # ee <- apply(ee, 2, genefu::rescale)
      splitix <- parallel::splitIndices(nx=ncol(ee), ncl=nthread)
      mcres <- parallel::mclapply(splitix, function(x, data) {
        return(apply(data[ , x, drop=FALSE], 2, scale))
      }, data=ee, mc.cores=nthread)
      ee <- do.call(cbind, mcres)
      exprs(eset.merged) <- ee
    },
    {
      stop("Unknown data standardization method")
    }
  )
    
  return (eset.merged)
}
         
duplicateFinder <- function (eset, var.genes=1000, dupl.cor=0.95, nthread=1) {
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
  
  
  
## end
