########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

## Note run this code from the STSP directory
## source(file.path("R", "STSP_isdb_data"))

## remove all existing objects from the workspace
rm(list=ls(all=TRUE))

require(inSilicoDb2)
require(genefu)

## this version of jetset allow for simple probe-gene mapping for affymetrix chips hgu95av2 (GPL8300), hgu133a (GPL96), hgu133b (GPL97), hgu133plus2 (GPL570), u133x3p (GPL1352)

source(file.path("R", "STSP_foo.R"))

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## directory where all the analysis results will be stored
saveres <- "cache"
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(file.path(saveres, "processed"))) { dir.create(file.path(saveres, "processed"), showWarnings=FALSE, recursive=TRUE) }
  


## number of cpu cores available for the analysis pipeline
## set to 'NULL' if all the available cores should be used
nbcore <- 16
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { ncore <- availcore }
options("mc.cores"=nbcore)

## probe gene mapping method for each platform
probegene.method <- rbind(
  c("MISC", "variance"),
  c("GPL8300", "jetset"),
  c("GPL96", "jetset"),
  c("GPL97", "jetset"),
  c("GPL570", "jetset"),
  c("GPL1352", "jetset"))
colnames(probegene.method) <- c("platform", "method")

## which subtype classification model
## scmgene, scmod1, scmod2, pam50, ssp2006, ssp2003
sbt.model <- "scmgene"

## remove duplicates
rm.dupl <- TRUE
## number of top variant genes to use when comparing gene expression profiles
topvarg <- 1000
## maximum correlation between gene expression profiles of different samples
duplcor <- 0.975

## read info about datasets
datasets <- read.csv(file.path("inst", "extdata", "datasets.csv"), stringsAsFactors=FALSE)
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
  message(sprintf("Get dataset %s", ddn))
  dataset.fn <- file.path(saveres, "processed", sprintf("%s_processed.RData", ddn))
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
    message("\tCheck feature and sample names")
    esets <- lapply(esets, checkNames)

    ## probe gene mapping
    message("\tProbe-gene mapping")
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
    message("")
  } else {
    load(file=dataset.fn)
  }
  
  ## special action may be required for some datasets
  switch (datasets[i, "Dataset.ID"],
    "GSE2034" = {
      if ("GSE5327" %in% names(eset.all)) {
        ## merge datasets GSE2034 and GSE5327
        message("\tMerge GSE2034 and GSE5327")
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
        message("\tMerge GSE2034 and GSE5327")
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
message("Update clinical information")
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
message("Subtype classification")
for (j in 1:length(eset.all)) {
  eset.all[[j]] <- subtypeClassification(eset=eset.all[[j]], model=sbt.model)
}

## merge expressionSet objects
message("Merge datasets")
eset.merged <- datasetMerging(esets=eset.all, nthread=nbcore)

## identify potential duplicated samples
duplicates <- duplicateFinder(eset=eset.merged, var.genes=topvarg, dupl.cor=duplcor)
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
if (rm.dupl) {
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
  }, y=rmix)
}

## log out from InSilicoDB
InSilicoLogout()
