



getDifferentialGenes2=function (scon2,clustering = NULL, groups = NULL, z.threshold = 3, 
          upregulated.only = F, verbose = T, plot = FALSE, n.genes.to.show = 10, 
          inner.clustering = FALSE, append.specifisity.metrics = TRUE, 
          append.auc = FALSE, n.cores = NULL) 
{
  samples=scon2$samples
  if (!is.null(clustering)) {
    groups <- clusters[[clustering]]$groups
  }
  if (is.null(groups)) {
    if (length(clusters) == 0) 
      stop("Either 'groups' be provided or clustering must be estimated (see 'findCommunities')")
    groups <- clusters[[1]]$groups
  }
  groups %>%  as.factor()
  if (class(samples[[1]]) != "Pagoda2") 
    stop("Only Pagoda objects are supported for marker genes")
  n.jobs <- if (is.null(n.cores)) 
    .self$n.cores
  else n.cores
  de.genes <- conos:::getDifferentialGenesP2(samples, groups = groups, 
                                             z.threshold = z.threshold, upregulated.only = upregulated.only, 
                                             verbose = verbose, n.cores = n.jobs)
  de.genes <- de.genes[levels(groups)]
  if (plot) {
    plotDEGenes(de.genes, samples, groups = groups, n.genes.to.show = n.genes.to.show, 
                inner.clustering = inner.clustering)
  }
  if (append.specifisity.metrics) {
    lapply.func <- if (verbose) 
      function(...) pbapply::pblapply(..., cl = n.jobs)
    else function(...) papply(..., n.cores = n.jobs)
    if (verbose) 
      cat("Estimating specifisity metrics\n")
    cm.merged <- lapply(samples, conos:::getRawCountMatrix, transposed = T) %>% 
      conos:::mergeCountMatrices(transposed = T)
    cm.merged=cm.merged[names(groups),]
    
    if (length(intersect(rownames(cm.merged), names(groups))) != 
        nrow(cm.merged)) 
      stop("`groups` must contain values for all cells in the samples")
    if (nrow(cm.merged) < length(groups)) {
      groups %<>% .[rownames(cm.merged)]
    }
    de.genes %<>% names() %>% setNames(., .) %>% lapply.func(function(n) appendSpecifisityMetricsToDE2(de.genes[[n]], 
                                                                                                       groups, n, p2.counts = cm.merged, append.auc = append.auc))
  }
  if (verbose) 
    cat("All done!\n")
  return(de.genes)
}




appendSpecifisityMetricsToDE2=function (de.df, clusters, cluster.id, p2.counts, low.expression.threshold = 0, 
                                        append.auc = FALSE) 
{
  cluster.mask <- setNames(clusters == cluster.id, names(clusters))
  gs=intersect(de.df$Gene,colnames(p2.counts))
  
  de.df=de.df[match(gs,de.df$Gene),]
  
  counts.bin <- (p2.counts[names(cluster.mask), gs, 
                           drop = F] > low.expression.threshold)
  counts.bin.sums <- Matrix::colSums(counts.bin)
  counts.bin.clust.sums <- Matrix::colSums(counts.bin & cluster.mask)
  if (append.auc) {
    if (requireNamespace("pROC", quietly = TRUE)) {
      de.df$AUC <- apply(counts.bin, 2, function(col) pROC::auc(as.integer(cluster.mask), 
                                                                as.integer(col)))
    }
    else {
      warning("You have to install pROC package to use append.auc")
    }
  }
  de.df$Specifisity <- (length(cluster.mask) - counts.bin.sums)/(length(cluster.mask) - 
                                                                   counts.bin.clust.sums)
  de.df$Precision <- counts.bin.clust.sums/counts.bin.sums
  de.df$ExpressionFraction <- Matrix::colMeans(counts.bin[cluster.mask, 
                                                          , drop = F])
  return(de.df)
}

