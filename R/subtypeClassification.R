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
  sbt.conv <- rbind(c("ER-/HER2-", "Basal"),
    c("HER2+", "Her2"),
    c("ER+/HER2- High Prolif", "LumB"),
    c("ER+/HER2- Low Prolif", "LumA")
  )
  colnames(sbt.conv) <- c("SCM.nomenclature", "SSP.nomenclature")
    
  sbtn.ssp <- c("Basal", "Her2", "LumB", "LumA", "Normal")
  sbtn2.ssp <- c("Basal", "Her2", "Lums", "LumB", "LumA", "Normal")
  
  datage <- t(Biobase::exprs(eset))   
  annotge <- cbind("probe"=rownames(Biobase::fData(eset)), "EntrezGene.ID"=stripWhiteSpace(as.character(Biobase::fData(eset)[ , "ENTREZID"])))
  rownames(annotge) <- stripWhiteSpace(as.character(annotge[ , "probe"]))
  
  ## SCM family
  if (model %in% c("scmgene", "scmod1", "scmod2")) {
    switch(model,
      "scmgene" = {
        sbts <- genefu::subtype.cluster.predict(sbt.model=genefu::scmgene.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype2", "subtype.proba2")]
      },
      "scmod1" = {
        sbts <- genefu::subtype.cluster.predict(sbt.model=genefu::scmod1.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype2", "subtype.proba2")]
      },
      "scmod2" = {
        sbts <- genefu::subtype.cluster.predict(sbt.model=genefu::scmod2.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype2", "subtype.proba2")]
      }
    )
    names(sbts) <- c("subtype", "subtype.proba")
    ## use SSP nomenclature
    ## update subtype calling
    ss <- factor(x=sbts$subtype)
    levels(ss)[match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], levels(ss))] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
    sbts$subtype <- factor(as.character(ss), levels=sbt.conv[ , "SSP.nomenclature"])
    ## update subtype probabilities
    iix <- match(sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SCM.nomenclature"], colnames(sbts$subtype.proba))
    colnames(sbts$subtype.proba)[iix] <- sbt.conv[!is.na(sbt.conv[ , "SCM.nomenclature"]), "SSP.nomenclature"]
    ## compute crisp classification
    sbts$subtype.crisp <- t(apply(sbts$subtype.proba, 1, function (x) {
      xx <- array(0, dim=length(x), dimnames=list(names(x)))
      xx[which.max(x)] <- 1
      return (xx)
    }))
    ## merge LumA and LumB
    lums.proba <- apply(sbts$subtype.proba[ , c("LumB", "LumA"), drop=FALSE], 1, sum, na.rm=TRUE)
    sbts$subtype.proba <- cbind(sbts$subtype.proba, "Lums"=lums.proba)
    lums.crisp <- as.numeric(is.element(sbts$subtype, c("LumA", "LumB")))
    sbts$subtype.crisp <- cbind(sbts$subtype.crisp, "Lums"=lums.crisp)
    ## reorder columns
    ss <- sbtn2.ssp[is.element(sbtn2.ssp, colnames(sbts$subtype.proba))]
    sbts$subtype.proba <- sbts$subtype.proba[ , ss, drop=FALSE]
    sbts$subtype.crisp <- sbts$subtype.crisp[ , ss, drop=FALSE]
    ## set the proper names
    names(sbts$subtype) <- rownames(sbts$subtype.proba) <- rownames(sbts$subtype.crisp)<- rownames(datage)
  }
  
  ## SSP family
  if (model %in% c("ssp2003", "ssp2006", "pam50")) {
    switch(model,
      "pam50" = {
        sbts <- genefu::intrinsic.cluster.predict(sbt.model=genefu::pam50.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype", "subtype.proba")]
      },
      "ssp2006" = {
        sbts <- genefu::intrinsic.cluster.predict(sbt.model=genefu::ssp2006.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype", "subtype.proba")]
      },
      "ssp2003" = {
        sbts <- genefu::intrinsic.cluster.predict(sbt.model=genefu::ssp2003.robust, data=datage, annot=annotge, do.mapping=TRUE)[c("subtype", "subtype.proba")]
      }
    )
    sbts$subtype <- factor(as.character(sbts$subtype), levels=sbtn.ssp)
    ## compute crisp classification
    sbts$subtype.crisp <- t(apply(sbts$subtype.proba, 1, function (x) {
      xx <- array(0, dim=length(x), dimnames=list(names(x)))
      xx[which.max(x)] <- 1
      return (xx)
    }))
    ## merge LumA and LumB
    lums.proba <- apply(sbts$subtype.proba[ , c("LumB", "LumA"), drop=FALSE], 1, sum, na.rm=TRUE)
    sbts$subtype.proba <- cbind(sbts$subtype.proba, "Lums"=lums.proba)
    lums.crisp <- as.numeric(is.element(sbts$subtype, c("LumA", "LumB")))
    sbts$subtype.crisp <- cbind(sbts$subtype.crisp, "Lums"=lums.crisp)
    ## reorder columns
    ss <- sbtn2.ssp[is.element(sbtn2.ssp, colnames(sbts$subtype.proba))]
    sbts$subtype.proba <- sbts$subtype.proba[ , ss, drop=FALSE]
    sbts$subtype.crisp <- sbts$subtype.crisp[ , ss, drop=FALSE]
    ## set the proper names
    names(sbts$subtype) <- rownames(sbts$subtype.proba) <- rownames(sbts$subtype.crisp)<- rownames(datage)
  }
  
  ## merge clinical information and subtype classification
  eset <- setSubtype(eset=eset, subtype.class=sbts$subtype, subtype.crisp=sbts$subtype.crisp, subtype.fuzzy=sbts$subtype.proba)
  return (eset)
}

