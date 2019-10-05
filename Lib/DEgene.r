


getDifferentialGenes2=function (cm,type = "counts", clusterType = NULL, groups = NULL, 
          name = "customClustering", z.threshold = 3, upregulated.only = FALSE, 
          verbose = FALSE) 
{
  "Determine differentially expressed genes, comparing each group against all others using Wilcoxon rank sum test\n\n       - type data type (currently only default 'counts' is supported)\n\n       - clusterType optional cluster type to use as a group-defining factor\n\n       - groups explicit cell group specification - a named cell factor (use NA in the factor to exclude cells from the comparison)\n\n       - name name slot to store the results in\n\n       - z.threshold minimal absolute Z score (adjusted) to report\n\n       - upregulated.only whether to report only genes that are expressed significantly higher in each group\n\n       - verbose verbose flag\n\n       return a list, with each element of the list corresponding to a cell group in the provided/used factor (i.e. factor levels); Each element of a list is a data frame listing the differentially epxressed genes (row names), with the following columns: Z - adjusted Z score, with positive values indicating higher expression in a given group compare to the rest; M - log2 fold change; highest- a boolean flag indicating whether the expression of a given gene in a given vcell group was on average higher than in every other cell group; fe - fraction of cells in a given group having non-zero expression level of a given gene;\n\n       Examples:\n\n\n         result <- r$getDifferentialGenes(groups=r$clusters$PCA$multilevel);\n\n         str(r$diffgenes)"
  if (is.null(groups)) {
    if (is.null(clusterType)) {
      cols <- clusters[[type]][[1]]
    }
    else {
      cols <- clusters[[type]][[clusterType]]
      if (is.null(cols)) {
        stop("clustering ", clusterType, " for type ", 
             type, " doesn't exist")
      }
    }
  }
  else {
    cols <- groups
  }

  if (!all(rownames(cm) %in% names(cols))) {
    warning("cluster vector doesn't specify groups for all of the cells, dropping missing cells from comparison")
  }
  valid.cells <- rownames(cm) %in% names(cols)[!is.na(cols)]
  if (!all(valid.cells)) {
    cm <- cm[valid.cells, ]
  }
  cols <- as.factor(cols[match(rownames(cm), names(cols))])
  cols <- as.factor(cols)
  if (verbose) {
    cat("running differential expression with ", length(levels(cols)), 
        " clusters ... ")
  }
  lower.lpv.limit <- -100
  xr <- pagoda2:::sparse_matrix_column_ranks(cm)
  grs <- pagoda2:::colSumByFac(xr, as.integer(cols))[-1, , drop = F]
  xr@x <- numeric(length(xr@x)) + 1
  gnzz <- pagoda2:::colSumByFac(xr, as.integer(cols))[-1, , drop = F]
  group.size <- as.numeric(tapply(cols, cols, length))[1:nrow(gnzz)]
  group.size[is.na(group.size)] <- 0
  gnz <- (group.size - gnzz)
  zero.ranks <- (nrow(xr) - diff(xr@p) + 1)/2
  ustat <- t((t(gnz) * zero.ranks)) + grs - group.size * (group.size + 
                                                            1)/2
  n1n2 <- group.size * (nrow(cm) - group.size)
  usigma <- sqrt(n1n2 * (nrow(cm) + 1)/12)
  usigma <- sqrt((nrow(cm) + 1 - (gnz^3 - gnz)/(nrow(cm) * 
                                                  (nrow(cm) - 1))) * n1n2/12)
  x <- t((ustat - n1n2/2)/usigma)
  if (verbose) {
    cat("adjusting p-values ... ")
  }
  x <- matrix(qnorm(pagoda2:::bh.adjust(pnorm(as.numeric(abs(x)), lower.tail = FALSE, 
                                    log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE), 
              ncol = ncol(x)) * sign(x)
  rownames(x) <- colnames(cm)
  colnames(x) <- levels(cols)[1:ncol(x)]
  if (verbose) {
    cat("done.\n")
  }
  log.gene.av <- log2(Matrix::colMeans(cm))
  group.gene.av <- pagoda2:::colSumByFac(cm, as.integer(cols))[-1, , 
                                                     drop = F]/(group.size + 1)
  log2.fold.change <- log2(t(group.gene.av)) - log.gene.av
  f.expressing <- t(gnzz/group.size)
  max.group <- max.col(log2.fold.change)
  if (upregulated.only) {
    ds <- lapply(1:ncol(x), function(i) {
      z <- x[, i]
      vi <- which(z >= z.threshold)
      r <- data.frame(Z = z[vi], M = log2.fold.change[vi, 
                                                      i], highest = max.group[vi] == i, fe = f.expressing[vi, 
                                                                                                          i])
      rownames(r) <- rownames(x)[vi]
      r <- r[order(r$Z, decreasing = T), ]
      r
    })
  }
  else {
    ds <- lapply(1:ncol(x), function(i) {
      z <- x[, i]
      vi <- which(abs(z) >= z.threshold)
      r <- data.frame(Z = z[vi], M = log2.fold.change[vi, 
                                                      i], highest = max.group[vi] == i, fe = f.expressing[vi, 
                                                                                                          i])
      rownames(r) <- rownames(x)[vi]
      r <- r[order(r$Z, decreasing = T), ]
      r
    })
  }
  names(ds) <- colnames(x)
 
  return(ds)
}


matCSC <- as(t(tmp), "dgCMatrix")
f=getDifferentialGenes2(matCSC,groups=conosCluster)
