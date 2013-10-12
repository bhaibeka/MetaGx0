########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

`checkNames` <- 
function (eset) {
  ## check feature and sample names in esets
  #
  # Args:
  #   esets: expressionSet object
  #
  # Returns:
  #   esets: updated expressionSet object with atching feature and sample names
  
  require(Biobase)
  ## check feature names
  check.feature <- intersect(rownames(Biobase::exprs(eset)), rownames(Biobase::fData(eset)))
  if (length(check.feature) == 0) {
    warning("Names of features do not match between expressions and annotations")
    return (NULL)
  } else {
    if (length(check.feature) != nrow(Biobase::exprs(eset)) || length(check.feature) != nrow(Biobase::fData(eset))) {
      warning("Some features are missing between expressions and annotations")
    }
  }
  ## check sample names
  check.sample <- intersect(colnames(Biobase::exprs(eset)), rownames(Biobase::pData(eset)))
  if (length(check.sample) == 0) {
    warning("Names of samples do not match between expressions and phenotypes")
    return (NULL)
  } else {
    if (length(check.sample) != ncol(Biobase::exprs(eset)) || length(check.sample) != nrow(Biobase::pData(eset))) {
      warning("Some samples are missing between expressions and phenotypes")
    }
  }
  Biobase::exprs(eset) <- Biobase::exprs(eset)[check.feature, check.sample, drop=FALSE]
  Biobase::fData(eset) <- Biobase::fData(eset)[check.feature, , drop=FALSE]
  Biobase::pData(eset) <- Biobase::pData(eset)[check.sample, , drop=FALSE]
  Biobase::pData(eset)[ , "samplename"] <- rownames(Biobase::pData(eset))
  return (eset)
}


`probeGeneMapping` <- 
function (eset, platform=c("MISC", "GPL8300", "GPL96", "GPL97", "GPL570", "GPL1352"), method=c("variance", "jetset")){
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


`subtypeClassification` <- 
function (eset, model=c("scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003")) {
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


`datasetMerging` <- 
function (esets, method=c("union", "intersect"), standardization=c("quantile", "robust.scaling", "scaling", "none"), nthread=1) {
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

`subtypeAssociation` <- 
function (eset, geneid, boxp=TRUE, subtype.col, resdir, nthread=1) {
  ## assess association between gene expression and subtypes
  #
  # Arga:
  #   eset: an expressionSet object
  #   gene: vector of Entrez Gene IDs. If missing, all genes will be considered.
  #
  # Returns
  #   list containing p-values and effect sizes
  
  if (class(eset) != "ExpressionSet") {
    stop("Handling list of expressionSet objects is not implemented yet")
  }
  
  if (missing(geneid)) {
    gened <- as.character(Biobase::fData(eset)[ , "ENTREZID"])
  }
  
  ## for a single expressionSet object
  
  ## extract subtypes
  sbts <- Biobase::pData(eset)[ , "subtype"]
  if (sum(table(sbts) > 3) < 2) {
    warning("Not enough tumors in each subtype")
    return(NULL)
  }
  sbtu <- levels(sbts)
  if (missing(subtype.col)) {
    subtype.col <- rainbow(length(sbtu), alpha=0.6)
  } else {
    if (length(subtype.col) < length(sbtu)) {
      stop(sprintf("Not enough color for %i subtypes", length(sbtu)))
    }  
  }
  
  ## extract genes
  gid <- intersect(geneid, as.character(Biobase::fData(eset)[ , "ENTREZID"]))
  if (length(gid) == 0) {
    stop("Genes not in the expressionSet object")
  }
  if (length(gid) < length(geneid)) {
    warning(sprintf("%i/%i genes were present in the expressionSet object", length(gid), length(geneid)))
  }
  
  splitix <- parallel::splitIndices(nx=length(gid), ncl=nthread)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, ...) {    
    pp <- lapply(gid[x], function (gid, eset, sbts, boxp, resdir) {
      ## gene symbol
      gsymb <- Biobase::fData(eset)[match(gid, Biobase::fData(eset)[ , "ENTREZID"]), "SYMBOL"]
      xx <- Biobase::exprs(eset)[paste("geneid", gid, sep="."), ]
      ## kruskal-wallis test
      kt <- kruskal.test(x=xx, g=sbts)$p.value
      ## pairwise wilcoxon test
      wt <- matrix(NA, nrow=length(sbtu), ncol=length(sbtu), dimnames=list(sbtu, sbtu))
      wt1 <- pairwise.wilcox.test(x=xx, g=sbts, p.adjust.method="none", paired=FALSE, alternative="greater")$p.value
      wt2 <- pairwise.wilcox.test(x=xx, g=sbts, p.adjust.method="none", paired=FALSE, alternative="less")$p.value
      nix <- !is.na(wt1)
      wt[rownames(wt1), colnames(wt1)][nix] <- wt1[nix]
      nix <- !is.na(t(wt2))
      wt[colnames(wt2), rownames(wt2)][nix] <- t(wt2)[nix]
      diag(wt) <- 1
      if (boxp) {
        pdf(file.path(resdir, sprintf("subtype_association_boxplot_%s.pdf", gsymb)))
        par(las=2, mar=c(5, 4, 4, 2) + 0.1, xaxt="n")
        # dd <- c(list(" "=NA), list("  "=NA), dd2)
  			mp <- boxplot(xx ~ sbts, las=3, outline=FALSE, ylim=c(-2,2), main=sprintf("%s", gsymb), col=subtype.col)
        axis(1, at=1:length(mp$names), tick=TRUE, labels=T)
        text(x=1:length(mp$names), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=mp$names, srt=45, xpd=NA, font=2, col=c("black"))
        text(x=1:length(mp$names), y=par("usr")[3], pos=3, labels=sprintf("n=%i", table(sbts)[mp$names]), col=c("black"))
        # title(sub=sprintf("Kruskal-Wallis p-value = %.1E", kt))
        # legend("topleft", legend=c(paste(mp$names, sprintf("(%i)", table(sbts)[mp$names]), sep=" "), "", sprintf("K-W p = %.1E", kt)), bty="n", cex=0.8)
        legend("topleft", legend=sprintf("Kruskal-Wallis p-value = %.1E", kt), bty="n")
        dev.off()
      }
      return(list("kruskal.pvalue"=kt, "wilcoxon.pvalue"=wt))
    }, eset=eset, sbts=sbts, boxp=boxp, resdir=resdir)
  }, gid=gid, eset=eset, sbts=sbts, boxp=boxp, resdir=resdir)
  pp <- do.call(c, mcres)
  names(pp) <- Biobase::fData(eset)[match(gid, Biobase::fData(eset)[ , "ENTREZID"]), "SYMBOL"]
  
  ## write spreadsheets
  dd <- sapply(pp, function(x) { return(x[[1]])})
  dd <- data.frame("Kruskal.Wallis.pvalue"=dd, "Kruskal.Wallis.fdr"=p.adjust(dd, method="fdr"), Biobase::fData(eset)[match(gid, Biobase::fData(eset)[ , "ENTREZID"]), ])
  write.csv(dd, file=file.path(resdir, "subtype_association_kruskal.csv"))
  mapply(function(x, y, resdir) {
    write.csv(x, file=file.path(resdir, sprintf("subtype_association_wilcoxon_%s.csv", y)))
  }, x=lapply(pp, function(x) { return(x[[2]])}), y=names(pp), resdir=resdir)
  
  return(pp)
}

`runPipeline` <- 
function (sbt.model=c("scmgene", "scmod2", "scmod1", "pam50", "ssp2006", "ssp2003"), resdir="cache", probegene.method, remove.duplicates=TRUE, topvar.genes=1000, duplicates.cor=0.975, datasets, nthread=1, verbose=TRUE) {  

  badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

  ## directory where all the analysis results will be stored
  if(!file.exists(resdir)) { dir.create(resdir, showWarnings=FALSE, recursive=TRUE) }
  if(!file.exists(file.path(resdir, "processed"))) { dir.create(file.path(resdir, "processed"), showWarnings=FALSE, recursive=TRUE) }

  ## number of cpu cores available for the analysis pipeline
  ## set to 'NULL' if all the available cores should be used
  availcore <- parallel::detectCores()
  if (nthread > availcore) { nthread <- availcore }
  options("mc.cores"=nthread)

  ## probe gene mapping method for each platform
  if (missing(probegene.method)) {
    probegene.method <- rbind(
      c("MISC", "variance"),
      c("GPL8300", "jetset"),
      c("GPL96", "jetset"),
      c("GPL97", "jetset"),
      c("GPL570", "jetset"),
      c("GPL1352", "jetset"))
    colnames(probegene.method) <- c("platform", "method")
  }

  ## which subtype classification model
  ## scmgene, scmod1, scmod2, pam50, ssp2006, ssp2003
  sbt.model <- match.arg(sbt.model)

  ## number of top variant genes to use when comparing gene expression profiles

  ## maximum correlation between gene expression profiles of different samples

  ## read info about datasets
  if (missing(datasets)) {
    datasets <- read.csv(system.file(file.path("extdata", "datasets.csv"), package="bcmeta"), stringsAsFactors=FALSE)
  }
  datasets <- datasets[datasets[ , "Include"], , drop=FALSE]

  # clinical info
  clin.info <- c("samplename", "id", "series", "dataset", "age", "node", "er", "e.rfs", "e.os", "e.dmfs", "grade", "size", "pgr", "her2", "t.rfs", "t.os", "t.dmfs", "treatment", "tissue")

  ## log in InSilicoDB
  InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")

  ## clinical information
  cinfo <- clin.info

  ## download all datasets
  eset.all <- NULL
  for (i in 1:nrow(datasets)) {
    ddn <- as.character(datasets[i, "Dataset.ID"])
    if (verbose) {
      message(sprintf("Get dataset %s", ddn))
    }
    dataset.fn <- file.path(resdir, "processed", sprintf("%s_processed.RData", ddn))
    if (!file.exists(dataset.fn)) {
      ## get dataset
      # inSilicoDb2::getCurationInfo(dataset=as.character(datasets[i, "Dataset.ID"]))
      platf <- inSilicoDb2::getPlatforms(dataset=ddn)
      esets <- inSilicoDb2::getDatasets(dataset=ddn, norm=as.character(datasets[i, "Normalization"]), curation=datasets[i, "Curation.ID"], features="PROBE")
      if (is.null(unlist(esets))) {
        stop("Rerun the script when data are ready to download from InSilicoDB")
      }
      platf <- unlist(platf[!sapply(esets, is.null)])
      esets <- esets[!sapply(esets, is.null)]
      ## check sample and feature names
      if (verbose) {
        message("\tCheck feature and sample names")
      }
      esets <- lapply(esets, checkNames)

      ## probe gene mapping
      if (verbose) {
        message("\tProbe-gene mapping")
      }
      platf2 <- platf
      platf2[!is.element(platf, probegene.method[ , "platform"])] <- "MISC"
      esets <- mapply(probeGeneMapping, eset=esets, platform=platf2, method=probegene.method[match(platf2, probegene.method[ , "platform"]), "method"])
  
      ## merge datasets if they are on different platforms
      eset <- platformMerging(esets)
      ## eset is a gene-centric expressionSet object
  
      ## format clinical info
      colnames(Biobase::phenoData(eset)@data) <- gsub(" ", ".", colnames(Biobase::phenoData(eset)@data))
      cinfo <- intersect(cinfo, colnames(Biobase::phenoData(eset)@data))
      ## transform factors into characters
      Biobase::phenoData(eset)@data <- data.frame(apply(Biobase::phenoData(eset)@data, 2, function(x) {
        if (is.factor(x)) { x <- as.character(x) }
        return(x)
        }), stringsAsFactors=FALSE)
      Biobase::phenoData(eset)@data[Biobase::phenoData(eset)@data == "NA"] <- NA
      save(list=c("eset"), compress=TRUE, file=dataset.fn)
      if (verbose) {
        message("")
      }
    } else {
      load(file=dataset.fn)
    }
  
    ## special action may be required for some datasets
    switch (datasets[i, "Dataset.ID"],
      "GSE2034" = {
        if ("GSE5327" %in% names(eset.all)) {
          ## merge datasets GSE2034 and GSE5327
          if (verbose) {
            message("\tMerge GSE2034 and GSE5327")
          }
          ## update gene expression data
          cg <- union(rownames(Biobase::exprs(eset.all[["GSE5327"]])), rownames(Biobase::exprs(eset)))
          cs <- c(colnames(Biobase::exprs(eset.all[["GSE5327"]])), colnames(Biobase::exprs(eset)))
          ee <- matrix(NA, nrow=length(cg), ncol=length(cs), dimnames=list(cg, cs))
          ee[rownames(Biobase::exprs(eset.all[["GSE5327"]])), colnames(Biobase::exprs(eset.all[["GSE5327"]]))] <- Biobase::exprs(eset.all[["GSE5327"]])
          ee[rownames(Biobase::exprs(eset)), colnames(Biobase::exprs(eset))] <- Biobase::exprs(eset)
          Biobase::exprs(eset.all[["GSE5327"]]) <- ee
          ## update feature data
          cf <- c(colnames(Biobase::featureData(eset.all[["GSE5327"]])@data), colnames(Biobase::featureData(eset)@data))
          ee <- data.frame(matrix(NA, nrow=length(cg), ncol=length(cf), dimnames=list(cg, cf)))
          ee[rownames(Biobase::featureData(eset.all[["GSE5327"]])@data), colnames(Biobase::featureData(eset.all[["GSE5327"]])@data)] <- Biobase::featureData(eset.all[["GSE5327"]])@data
          ee[rownames(Biobase::featureData(eset)@data), colnames(Biobase::featureData(eset)@data)] <- Biobase::featureData(eset)@data
          Biobase::featureData(eset.all[["GSE5327"]]) <- ee
          ## update pheno data
          cp <- intersect(colnames(Biobase::phenoData(eset.all[["GSE5327"]])@data), colnames(Biobase::phenoData(eset)@data))
          ee <- data.frame(matrix(NA, nrow=length(cs), ncol=length(cp), dimnames=list(cs, cp)))
          ee[rownames(Biobase::phenoData(eset.all[["GSE5327"]])@data), cp] <- Biobase::phenoData(eset.all[["GSE5327"]])@data[ , cp]
          ee[rownames(Biobase::phenoData(eset)@data), cp] <- Biobase::phenoData(eset)@data[ , cp]
          Biobase::phenoData(eset.all[["GSE5327"]]) <- ee
          ## rename the merged dataset
          names(eset.all)[names(eset.all) == "GSE2034"] <- "GSE2034.GSE5327"
        } else {
          eset.all <- c(eset.all, eset)
          names(eset.all)[length(eset.all)] <- datasets[i, "Dataset.ID"]
        }
      },
      "GSE5327" = {
        if ("GSE2034" %in% names(eset.all)) {
          ## merge datasets GSE2034 and GSE5327
          if (verbose) {
            message("\tMerge GSE2034 and GSE5327")
          }
          ## update gene expression data
          cg <- union(rownames(Biobase::exprs(eset.all[["GSE2034"]])), rownames(Biobase::exprs(eset)))
          cs <- c(colnames(Biobase::exprs(eset.all[["GSE2034"]])), colnames(Biobase::exprs(eset)))
          ee <- matrix(NA, nrow=length(cg), ncol=length(cs), dimnames=list(cg, cs))
          ee[rownames(Biobase::exprs(eset)), colnames(Biobase::exprs(eset))] <- Biobase::exprs(eset)
          ee[rownames(Biobase::exprs(eset.all[["GSE2034"]])), colnames(Biobase::exprs(eset.all[["GSE2034"]]))] <- Biobase::exprs(eset.all[["GSE2034"]])
          Biobase::exprs(eset.all[["GSE2034"]]) <- ee
          ## update feature data
          cf <- intersect(colnames(Biobase::featureData(eset.all[["GSE2034"]])@data), colnames(Biobase::featureData(eset)@data))
          ee <- data.frame(matrix(NA, nrow=length(cg), ncol=length(cf), dimnames=list(cg, cf)))
          ee[rownames(Biobase::featureData(eset)@data), colnames(Biobase::featureData(eset)@data)] <- Biobase::featureData(eset)@data
          ee[rownames(Biobase::featureData(eset.all[["GSE2034"]])@data), colnames(Biobase::featureData(eset.all[["GSE2034"]])@data)] <- Biobase::featureData(eset.all[["GSE2034"]])@data
          Biobase::featureData(eset.all[["GSE2034"]])@data <- ee
          ## update pheno data
          cp <- intersect(colnames(Biobase::phenoData(eset.all[["GSE2034"]])@data), colnames(Biobase::phenoData(eset)@data))
          ee <- data.frame(matrix(NA, nrow=length(cs), ncol=length(cp), dimnames=list(cs, cp)))
          ee[rownames(Biobase::phenoData(eset.all[["GSE2034"]])@data), cp] <- Biobase::phenoData(eset.all[["GSE2034"]])@data[ , cp]
          ee[rownames(Biobase::phenoData(eset)@data), cp] <- Biobase::phenoData(eset)@data[ , cp]
          Biobase::phenoData(eset.all[["GSE2034"]])@data <- ee
          ## rename the merged dataset
          names(eset.all)[names(eset.all) == "GSE2034"] <- "GSE2034.GSE5327"
        } else {
          eset.all <- c(eset.all, eset)
          names(eset.all)[length(eset.all)] <- datasets[i, "Dataset.ID"]
        }
      },
      "GSE2109" = {
        ## select only breast tumors
        iix <- Biobase::phenoData(eset)@data[ , "Anatomical.site"] == "breast"
        Biobase::exprs(eset) <- Biobase::exprs(eset)[ , iix, drop=FALSE]
        Biobase::phenoData(eset)@data <- Biobase::phenoData(eset)@data[iix, , drop=FALSE]
        eset.all <- c(eset.all, eset)
        names(eset.all)[length(eset.all)] <- datasets[i, "Dataset.ID"]
      },
      {
        ## no action is necessary
        eset.all <- c(eset.all, eset)
        names(eset.all)[length(eset.all)] <- datasets[i, "Dataset.ID"]
      })
    ## eset.all contains a list of gene-centric expression sets
  }

  ## align clinical information
  if (verbose) {
    message("Update clinical information")
  }
  for (j in 1:length(eset.all)) {
    ## same order for all datasets
    Biobase::phenoData(eset.all[[j]])@data <- Biobase::phenoData(eset.all[[j]])@data[ , cinfo, drop=FALSE]
    ## ensure age, survival data are numeric
    tt <- intersect(colnames(Biobase::phenoData(eset.all[[j]])@data), c("age", "size", "er", "her2", "pgr", "grade", c(paste(c("t", "e"), rep(c("rfs", "dmfs", "os"), each=2), sep="."))))
    Biobase::phenoData(eset.all[[j]])@data[ , tt] <- data.frame(apply(Biobase::phenoData(eset.all[[j]])@data[ , tt, drop=FALSE], 2, as.numeric), stringsAsFactors=FALSE)
    ## create dfs data: use rfs when available, dmfs otherwise
    surv.time <- Biobase::phenoData(eset.all[[j]])@data[ , "t.rfs"]
    surv.time[is.na(Biobase::phenoData(eset.all[[j]])@data[ , "t.rfs"])] <- Biobase::phenoData(eset.all[[j]])@data[is.na(Biobase::phenoData(eset.all[[j]])@data[ , "t.rfs"]), "t.dmfs"]
    surv.event <- Biobase::phenoData(eset.all[[j]])@data[ , "e.rfs"]
    surv.event[is.na(Biobase::phenoData(eset.all[[j]])@data[ , "e.rfs"])] <- Biobase::phenoData(eset.all[[j]])@data[is.na(Biobase::phenoData(eset.all[[j]])@data[ , "e.rfs"]), "e.dmfs"]
    Biobase::phenoData(eset.all[[j]])@data <- cbind(Biobase::phenoData(eset.all[[j]])@data, "t.dfs"=surv.time, "e.dfs"=surv.event)
  }

  ## annotate with subtypes
  if (verbose) {
    message("Subtype classification")
  }
  for (j in 1:length(eset.all)) {
    eset.all[[j]] <- subtypeClassification(eset=eset.all[[j]], model=sbt.model)
  }

  ## merge expressionSet objects
  if (verbose) {
    message("Merge datasets")
  }
  eset.merged <- datasetMerging(esets=eset.all, nthread=nthread)

  ## identify potential duplicated samples
  duplicates <- duplicateFinder(eset=eset.merged, var.genes=topvar.genes, dupl.cor=duplicates.cor)
  ## annotate the separate esets and the merged eset
  tt <- sapply(duplicates, paste, collapse="///")
  ## merged eset
  Biobase::pData(eset.merged) <- cbind(Biobase::pData(eset.merged), "duplicates"=NA)
  Biobase::pData(eset.merged)[names(tt), "duplicates"] <- tt
  ## individual esets
  eset.all <- lapply(eset.all, function (x, y) {
    Biobase::pData(x) <- cbind(Biobase::pData(x), "duplicates"=NA)
    nn <- intersect(rownames(pData(x)), y)
    Biobase::pData(x)[nn, "duplicates"] <- y[nn]
    return(x)
  }, y=tt)

  ## remove duplicates
  if (remove.duplicates) {
    ## duplicates are removed by order of datasets
    ## select sample names to remove
    rmix <- duplicates
    ii <- 1
    while (length(rmix) > ii) {
      rmix <- rmix[!is.element(names(rmix), rmix[[ii]])]
      ii <- ii + 1
    }
    rmix <- unique(unlist(rmix))
    ## merged eset
    keepix <- setdiff(sampleNames(eset.merged), rmix)
    Biobase::exprs(eset.merged) <- Biobase::exprs(eset.merged)[ , keepix, drop=FALSE]
    Biobase::pData(eset.merged) <- Biobase::pData(eset.merged)[keepix, , drop=FALSE]
    ## individual esets
    eset.all <- lapply(eset.all, function (x, y) {
      keepix <- setdiff(sampleNames(x), y)
      Biobase::exprs(x) <- Biobase::exprs(x)[ , keepix, drop=FALSE]
      Biobase::pData(x) <- Biobase::pData(x)[keepix, , drop=FALSE]
      return(x)
    }, y=rmix)
  }

  ## log out from InSilicoDB
  InSilicoLogout()
  
  return(list("merged"=eset.merged, "each"=eset.all))
}

## end
