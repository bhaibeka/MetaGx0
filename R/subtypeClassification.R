########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################


`subtypeClassification` <- 
function (eset, model=c("scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003")) {
  # Classify GSE subtype and give out the probability of being that subtype
  #
  # Args:
  #   eset: expression set to classify
  #
  # Returns:
  #   eset: expression set with updated subtype information (subtype claling and subtype.probabilities)
  
  model <- match.arg(model)
  
  ## convert SCM to SSP nomenclature
  sbt.conv <- rbind(c("ER+/HER2- Low Prolif", "LumA"),
    c("ER+/HER2- High Prolif", "LumB"),
    c("HER2+", "Her2"),
    c("ER-/HER2-", "Basal"))
  colnames(sbt.conv) <- c("SCM.nomenclature", "SSP.nomenclature")
    
  sbtn2 <- c("LumA", "LumB", "Her2", "Basal", "Normal")
  
  datage <- t(Biobase::exprs(eset))   
  annotge <- cbind("probe"=rownames(Biobase::fData(eset)), "EntrezGene.ID"=stripWhiteSpace(as.character(Biobase::fData(eset)[ , "ENTREZID"])))
  rownames(annotge) <- stripWhiteSpace(as.character(annotge[ , "probe"]))
  
  switch(model,
    "scmgene" = {
      sbts <- genefu::subtype.cluster.predict(sbt.model=genefu::scmgene.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype2", "subtype.proba2")]
      names(sbts) <- c("subtype", "subtype.proba")
      ## use SSP nomenclature
      ## update subtype calling
      ss <- factor(x=sbts$subtype)
      levels(ss)[match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], levels(ss))] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
      sbts$subtype <- stripWhiteSpace(as.character(ss))
      ## update subtype probabilities
      iix <- match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], colnames(sbts$subtype.proba))
      colnames(sbts$subtype.proba)[iix] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
    },
    "scmod1" = {
      sbts <- genefu::subtype.cluster.predict(sbt.model=genefu::scmod1.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype2", "subtype.proba2")]
      names(sbts) <- c("subtype", "subtype.proba")
      ## use SSP nomenclature
      ## update subtype calling
      ss <- factor(x=sbts$subtype)
      levels(ss)[match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], levels(ss))] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
      sbts$subtype <- stripWhiteSpace(as.character(ss))
      ## update subtype probabilities
      iix <- match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], colnames(sbts$subtype.proba))
      colnames(sbts$subtype.proba)[iix] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
    },
    "scmod2" = {
      sbts <- genefu::subtype.cluster.predict(sbt.model=genefu::scmod2.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype2", "subtype.proba2")]
      names(sbts) <- c("subtype", "subtype.proba")
      ## use SSP nomenclature
      ## update subtype calling
      ss <- factor(x=sbts$subtype)
      levels(ss)[match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], levels(ss))] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
      sbts$subtype <- stripWhiteSpace(as.character(ss))
      ## update subtype probabilities
      iix <- match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], colnames(sbts$subtype.proba))
      colnames(sbts$subtype.proba)[iix] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
    },
    "pam50" = {
      sbts <- genefu::intrinsic.cluster.predict(sbt.model=genefu::pam50.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype", "subtype.proba")]
    },
    "ssp2006" = {
      sbts <- genefu::intrinsic.cluster.predict(sbt.model=genefu::ssp2006.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype", "subtype.proba")]
    },
    "ssp2003" = {
      sbts <- genefu::intrinsic.cluster.predict(sbt.model=genefu::ssp2003.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype", "subtype.proba")]
    },
    {
      stop("Unknown subtype classification model")
    })
  
  ## merge clinical information and subtype classification
  colnames(sbts$subtype.proba) <- paste("subtyproba", colnames(sbts$subtype.proba), sep=".")
  Biobase::pData(eset) <- cbind(Biobase::pData(eset), "subtype"=sbts$subtype, sbts$subtype.proba)
  
  return(eset)
}

