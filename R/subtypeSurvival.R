########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################


`subtypeSurvival` <- 
function (eset, geneid, plot=TRUE, weighted=TRUE, time.cens, condensed=TRUE, resdir="cache", nthread=1) {

  ######################
  
  concIndex <- function (stime, sevent, strat, tau, inpudata, cv.fold, itr=1000, seed=12345) {
    
    if (missing(tau))) { tau <- max(stime, na.rm=TRUE) }
    if (missing(strat)) { strat <- array(1, dim=length(stime), dimnames=list(names(stime))) }
    if (length(stime) != length(sevent) || length(stime) != length(strata) || length(stime) != nrow(inpudata)) { stop("Wrong dimensions dimensions for stime, sevent, strat and inpudata") }
    if (!missing(cv.fold) && (cv.fold <= 1 || cv.fold > nrow(inpudata))) { cv.fold <- nrow(inpudata) }
    if (!is.matrix(inpudata)) { inpudata <- cbind("score"=inpudata) }
    ## use Uno C index
    if (missing(cv.fold)) {
      cindex <- survC1::Inf.Cval(mydata=cbind("stime"=stime, "sevent"-sevent, inpudata), tau=tau, itr=itr, seed=seed)
  }
  
  ######################
  
  if (class(eset) != "ExpressionSet") {
    stop("Handling list of expressionSet objects is not implemented yet")
  }
  
  if (missing(geneid)) {
    geneid <- Biobase::fData(eset)[ , "ENTREZID"]
  }
  
  if (!file.exists(file.path(resdir))) { dir.create(file.path(resdir), showWarnings=FALSE, recursive=TRUE) }
  
  ## for a single expressionSet object
  
  survC1::inf.Cval
}


