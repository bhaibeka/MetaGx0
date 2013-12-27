########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################


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
    gened <- stripWhiteSpace(as.character(Biobase::fData(eset)[ , "ENTREZID"]))
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
  gid <- intersect(geneid, stripWhiteSpace(as.character(Biobase::fData(eset)[ , "ENTREZID"])))
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
  dd <- data.frame("Kruskal.Wallis.pvalue"=dd, "Kruskal.Wallis.fdr"=p.adjust(dd, method="fdr"), Biobase::fData(eset)[match(gid, Biobase::fData(eset)[ , "ENTREZID"]), ], stringsAsFactors=FALSE)
  write.csv(dd, file=file.path(resdir, "subtype_association_kruskal.csv"))
  mapply(function(x, y, resdir) {
    write.csv(x, file=file.path(resdir, sprintf("subtype_association_wilcoxon_%s.csv", y)))
  }, x=lapply(pp, function(x) { return(x[[2]])}), y=names(pp), resdir=resdir)
  
  return(pp)
}


