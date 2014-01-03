########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################


`subtypeSurvival` <- 
function (eset, geneid, plot=TRUE, weighted=TRUE, time.cens, condensed=TRUE, resdir="cache", nthread=1) {

  ######################
  
  concIndex <- function (x, stime, sevent, strat, weights, tau, alternative=c("two.sided", "less", "greater")) {
    if (missing(tau)) {
      tau <- max(stime, na.rm=TRUE)
    } else {
      ss <- survcomp::censor.time(surv.time=stime, surv.event=sevent, time.cens=tau)  
      stime <- ss[[1]]
      sevent <- ss[[2]]
    }
    if (missing(strat)) { strat <- array(1, dim=length(stime), dimnames=list(names(stime))) }
    if (length(stime) != length(sevent) || length(stime) != length(strat) || length(stime) != length(x)) { stop("stime, sevent, strat and x must have the same length") }
    rr <-  mRMRe::correlate(X=x, Y=Surv(stime, sevent), method="cindex", strata=strat, weights=weights)
    rr <- rr[c("estimate", "se", "lower", "upper", "p")]
    rr <- c(list("Dxy"= 2 * (rr[["estimate"]] - 0.5)), rr)
    names(rr) <- c("Dxy", "cindex", "se", "lower", "upper", "p.value")
    return (rr)
  }
  
  dIndex <- function (x, stime, sevent, strat, weights, tau, alternative=c("two.sided", "less", "greater")) {
    if (missing(tau)) {
      tau <- max(stime, na.rm=TRUE)
    } else {
      ss <- survcomp::censor.time(surv.time=stime, surv.event=sevent, time.cens=tau)  
      stime <- ss[[1]]
      sevent <- ss[[2]]
    }
    if (missing(strat)) { strat <- array(1, dim=length(stime), dimnames=list(names(stime))) }
    if (length(stime) != length(sevent) || length(stime) != length(strat) || length(stime) != length(x)) { stop("stime, sevent, strat and x must have the same length") }
    rr <-  survcomp::D.index(x=x, surv.time=stime, surv.event=sevent, strat=strat, weights=weights, method.test="logrank", na.rm=TRUE)
    rr <- rr[c("d.index", "coef", "se", "lower", "upper", "p.value")]
    names(rr) <- c("Dindex", "coef", "se", "lower", "upper", "p.value")
    return (rr)
  }
  
  ######################
  
  if (class(eset) != "ExpressionSet") {
    stop("Handling list of expressionSet objects is not implemented yet")
  }
  
  if (missing(geneid)) {
    geneid <- Biobase::fData(eset)[ , "ENTREZID"]
  }
  
  if (!file.exists(file.path(resdir))) { dir.create(file.path(resdir), showWarnings=FALSE, recursive=TRUE) }
  
  if (missing(time.cens)) {
    time.cens <- max(stime, na.rm=TRUE)
  }
  
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
  if (length(gid) == 0) {
    stop("Genes not in the expressionSet object")
  }
  if (length(gid) < length(geneid)) {
    warning(sprintf("%i/%i genes were present in the expressionSet object", length(gid), length(geneid)))
  }
  glabel <- Biobase::fData(eset)[gid, "SYMBOL"]
  glabel[is.na(glabel)] <- paste("ENTREZID", Biobase::fData(eset)[gid, "ENTREZID"][is.na(glabel)], sep=".")
  names(glabel) <- gid
  expr <- Biobase::exprs(eset)[gid, , drop=FALSE]
  stime <- Biobase::pData(eset)[ , "t.dfs"] / 365
  time.cens <- time.cens / 365
  sevent <- Biobase::pData(eset)[ , "e.dfs"]
  strat <- as.factor(Biobase::pData(eset)[ , "dataset"])
  
  ## concordance index
  splitix <- parallel::splitIndices(nx=length(gid), ncl=nthread)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, gid, expr, stime, sevent, strat, sbts.proba) {    
    ci <- lapply(gid[x], function (x, expr, stime, sevent, strat, sbts.proba) {
      res <- t(apply(sbts.proba, 2, function (w, xx, stime, sevent, strat) {
        return (unlist(concIndex(x=xx, stime=stime, sevent=sevent, strat=strat, weights=w, alternative="two.sided")))
      }, xx=expr[x, ], stime=stime, sevent=sevent, strat=strat))
      return (res)
    }, expr=expr, stime=stime, sevent=sevent, strat=strat, sbts.proba=sbts.proba)
  }, gid=gid, expr=expr, stime=stime, sevent=sevent, strat=strat, sbts.proba=sbts.proba)
  rr <- unlist(mcres, recursive=FALSE)
  names(rr) <- glabel
  ## save results
  dd <- lapply(rr, data.frame)
  if (condensed) {
    WriteXLS::WriteXLS(x="dd", ExcelFileName=file.path(resdir, sprintf("subtype_cindex.xls")), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
  } else {
    mapply(function(x, y, method, resdir) {
       WriteXLS::WriteXLS("x", ExcelFileName=file.path(resdir, sprintf("subtype_cindex_%s.xls", y)), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
     }, x=dd, y=names(dd), method=method, resdir=resdir)
  }
  cindices <- rr
  
  ## D index (hazard ratio)
  splitix <- parallel::splitIndices(nx=length(gid), ncl=nthread)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, gid, expr, stime, sevent, strat, sbts.proba) {    
    ci <- lapply(gid[x], function (x, expr, stime, sevent, strat, sbts.proba) {
      res <- t(apply(sbts.proba, 2, function (w, xx, stime, sevent, strat) {
         return (unlist(dIndex(x=xx, stime=stime, sevent=sevent, strat=strat, weights=w, alternative="two.sided")))
      }, xx=expr[x, ], stime=stime, sevent=sevent, strat=strat))
      return (res)
    }, expr=expr, stime=stime, sevent=sevent, strat=strat, sbts.proba=sbts.proba)
  }, gid=gid, expr=expr, stime=stime, sevent=sevent, strat=strat, sbts.proba=sbts.proba)
  rr <- unlist(mcres, recursive=FALSE)
  names(rr) <- glabel
  ## save results
  dd <- lapply(rr, data.frame)
  if (condensed) {
    WriteXLS::WriteXLS(x="dd", ExcelFileName=file.path(resdir, sprintf("subtype_gene_dindex.xls")), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
    ## per subtype
    dd2 <- lapply(sbtu, function (x, y) {
      res <- t(sapply(y, function (y, x) {
        return (y[x, ])
      }, x=x))
      return (data.frame(res))
    }, y=rr)
    names(dd2) <- sbtu
     WriteXLS::WriteXLS(x="dd2", ExcelFileName=file.path(resdir, sprintf("subtype_dindex.xls")), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
  } else {
    mapply(function(x, y, method, resdir) {
       WriteXLS::WriteXLS("x", ExcelFileName=file.path(resdir, sprintf("subtype_dindex_%s.xls", y)), AdjWidth=FALSE, BoldHeaderRow=FALSE, row.names=TRUE, col.names=TRUE, FreezeRow=1, FreezeCol=1)
     }, x=dd, y=names(dd), method=method, resdir=resdir)
  }
  dindices <- rr
  
  ## kaplan-meier survival curves
  if (plot) {
    figsize <- 6
    nc <- 3
    nr <- ceiling(ncol(sbts.proba) / nc)
    if (condensed) { pdf(file.path(resdir, "subtype_surv_curves.pdf"), height=nr * figsize, width=nc * figsize) }
    lapply(gid, function (x, expr, stime, sevent, strat, condensed, sbts.proba, glabel, subtype.col, resdir, nc, nr) {
      if (!condensed) { pdf(file.path(resdir, sprintf("subtype_surv_curves_%s.pdf", x$symbol)), height=nr * figsize, width=nc * figsize) }
      par(mfrow=c(nr, nc))
      cc <- quantile(expr[x, ], probs=c(0, 0.33, 0.66, 1), na.rm=TRUE)
      xx <- factor(cut(x=expr[x, ], breaks=cc, labels=FALSE))
      for (i in 1:ncol(sbts.proba)) {
        w <- sbts.proba[ , i]
        dd <- data.frame("stime"=stime, "sevent"=sevent, "risk"=xx, strangsAsFactors=FALSE)
        survcomp::km.coxph.plot(formula.s=Surv(stime, sevent) ~ risk, data.s=dd, weight.s=w, sub.s="all", x.label="time (years)", y.label="probability of disease-free survival", main.title=sprintf("%s\n%s", glabel[x], colnames(sbts.proba)[i]), leg.text=paste(c("Low", "Intermediate", "High"), "     ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen", "darkred"), .lty=c(1,1,1), show.n.risk=TRUE, n.risk.step=2, n.risk.cex=0.85, bty="n", leg.bty="n", verbose=FALSE)
      }     
      if (!condensed) { dev.off() }
    }, expr=expr, stime=stime, sevent=sevent, strat=strat, condensed=condensed, sbts.proba=sbts.proba, glabel=glabel, subtype.col=subtype.col, resdir=resdir, nc=nc, nr=nr)
    if (condensed) { dev.off() }
  }

}


