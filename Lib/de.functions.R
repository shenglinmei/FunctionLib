
#' @param correction.method 'varianceweighted' or 'mean' specifies way to merge the fold changes from different cell types
getCorrectionVector_raw_list <- function(conObj, groups=NULL, sampleGroups=NULL, cooksCutoff=FALSE, independentFiltering = FALSE, n.cores=1, cluster.sep.chr = '+', return.details=FALSE,de.init=NULL,exclude.celltypes=c(),correction.method='varianceweighted',reflevel=NULL) {
    ## conObj <- con; sampleGroups <- sampleGroupsTvsW; cooksCutoff<-F; independentFiltering <- F; n.cores <- 1; cluster.sep.chr <- '+'
    ## groups <- as.factor(jcl3.coarse)
    ## Check arguments
    if ( is.null(groups) ) stop('groups must be specified');
    if ( is.null(sampleGroups) ) stop('sampleGroups must be specified')
    if ( class(sampleGroups) != 'list' ) stop('sampleGroups must be a list');
    if ( length(sampleGroups) != 2 ) stop('sampleGroups must be of length 2');
    if ( ! all(unlist(lapply(sampleGroups, function(x) class(x) == 'character'))) )
        stop('sampleGroups must be a list of character vectors');
    if ( ! all(unlist(lapply(sampleGroups, function(x) length(x) > 0))) )
        stop('sampleGroups entries must be on length greater or equal to 1')
    
    ## todo: check samplegrousp are named
    if(is.null(names(sampleGroups))) stop('sampleGroups must be named')
    if(class(groups) != 'factor') stop('groups must be a factor')
    if(any(grepl(cluster.sep.chr, names(conObj$samples),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any sample name')
    if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE))) 
        stop('cluster.sep.chr must not be part of any cluster name')
    if(is.null(reflevel)) stop('reference level is not defined')
    ## Main function
    if(is.null(de.init)) {
        de.init <- getPerCellTypeDE_raw_list(conObj, groups=groups, sampleGroups=sampleGroups,
                                    cooksCutoff=cooksCutoff, independentFiltering=independentFiltering,
                                    n.cores=n.cores, cluster.sep.chr=cluster.sep.chr,return.details=FALSE,
                                    reflevel = reflevel);
    }
    allfcs <- lapply(de.init, function(x) {
        if(!is.error(x)) {
            fc <- x$log2FoldChange
            names(fc) <- rownames(x)
            fc
        } else {
            NULL
        }
    })
    allfcs <- allfcs[!unlist(lapply(allfcs, is.null))]
    genes <- Reduce(intersect, lapply(allfcs, names))
    ## Matrix of fold changes
    fc.mat <- do.call(rbind, lapply(allfcs, function(x) {x[genes]}))
    fc.mat <- fc.mat[!rownames(fc.mat) %in% exclude.celltypes,]
    if (correction.method == 'mean') {
        correction <- apply(fc.mat, 2, mean, na.rm=TRUE)
    } else if (correction.method == 'varianceweighted') {
        mu <- apply(fc.mat, 2, mean, na.rm=TRUE)
        var <- apply(fc.mat, 2, function(x) {var(x,na.rm=TRUE)})
        weight <- 1 - pchisq(q=var,df=nrow(fc.mat)-1)
        correction <- mu * weight
    } else {
        error(paste0('unknown correction method: ', correction.method))
    }
    ## return
    if (!return.details) {
        correction
    } else {
        list(correction.vector=correction,de.init=de.init)
    }        
}








getPerCellTypeDECorrected_raw_list <- function(conObj, groups=NULL, sampleGroups=NULL, cooksCutoff = FALSE,
                             independentFiltering = FALSE, n.cores=1,cluster.sep.chr = '+',
                             correction=NULL, return.details=TRUE, reflevel=NULL) {
    ## Check arguments
    if ( is.null(correction) ) stop("Correction can't by null'")
    if ( is.null(groups) ) stop('groups must be specified');
    if ( is.null(sampleGroups) ) stop('sampleGroups must be specified')
    if ( class(sampleGroups) != 'list' ) stop('sampleGroups must be a list');
    if ( length(sampleGroups) != 2 ) stop('sampleGroups must be of length 2');
    if ( ! all(unlist(lapply(sampleGroups, function(x) class(x) == 'character'))) )
        stop('sampleGroups must be a list of character vectors');
    if ( ! all(unlist(lapply(sampleGroups, function(x) length(x) > 0))) )
        stop('sampleGroups entries must be on length greater or equal to 1')
    if ( ! all(unlist(lapply(sampleGroups, function(x) {all(x %in% names(conObj))}))) )
        stop('sampleGroups entries must be names of samples in the conos object')
    if (is.null(reflevel)) stop('reference level is not defined')
    ## todo: check samplegrousp are named
    if(is.null(names(sampleGroups))) stop('sampleGroups must be named')
    if(class(groups) != 'factor') stop('groups must be a factor')
    if(any(grepl(cluster.sep.chr, names(conObj),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any sample name')
    if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE))) 
        stop('cluster.sep.chr must not be part of any cluster name')
    ## Generate a summary dataset collapsing the cells of the same type in each sample
    ## and merging everything in one matrix
    samples.used <- unlist(sampleGroups)
    ## Generate an aggregated matrix
    raw.mats <- lapply(conObj[samples.used], function(p2) {
        p2
    })
    common.genes <- Reduce(intersect,lapply(raw.mats, colnames))
    raw.mats <- lapply(raw.mats, function(x) {x[,common.genes]})
    aggr2 <- lapply(raw.mats, function(x) {
        g1 <- groups[intersect(names(groups), rownames(x))]
        aggr <- Matrix.utils::aggregate.Matrix(x, g1)
        aggr
    })
    aggr2 <- lapply(names(aggr2), function(n) {
        x <- aggr2[[n]]
        rownames(x) <- paste0(n,cluster.sep.chr,rownames(aggr2[[n]]))
        x
    })
    aggr2 <- t(do.call(rbind, aggr2))
    rm(raw.mats); gc()
    ## For every cell type get differential expression results
    de.res <- mclapply(namedLevels(groups), function(l) {
        try({
            ## Get count matrix
            cm <- aggr2[,strpart(colnames(aggr2),cluster.sep.chr,2,fixed=TRUE) == l]
            ## Generate metadata
            meta <- data.frame(
                sample.id= colnames(cm),
                group= as.factor(unlist(lapply(colnames(cm), function(y) {
                    y <- strpart(y,cluster.sep.chr,1,fixed=TRUE)
                    names(sampleGroups)[unlist(lapply(sampleGroups,function(x) any(x %in% y)))]
                })))
            )
            meta$group <- relevel(meta$group, ref=reflevel)
            if (length(unique(as.character(meta$group))) < 2)
                stop('The cluster is not present in both conditions')
            library(DESeq2)
            dds1 <- DESeqDataSetFromMatrix(cm,meta,design=~group)
            dds1 <- estimateSizeFactors(dds1)
            sf <- sizeFactors(dds1)
            if(!(all(rownames(cm) %in% names(correction)) & all(names(correction) %in% rownames(cm))))
                stop('incompatible matrices')
            nf.tmp <- matrix(rep(sf, nrow(cm)),nrow=nrow(cm),byrow=TRUE)
            rownames(nf.tmp) <- rownames(cm);
            colnames(nf.tmp) <- colnames(cm)
            gene.scale.factors <- 2^(correction[rownames(nf.tmp)])
            baselevel <- levels(colData(dds1)$group)[1]
            x <- do.call(cbind, lapply(colData(dds1)$group, function(x) {
                if (x == baselevel) {
                    rep(1, length(gene.scale.factors))
                } else {
                    gene.scale.factors
                }
            }))
            rownames(x) <- rownames(nf.tmp);
            colnames(x) <- colnames(nf.tmp)
            nf.tmp <- nf.tmp * x
            x2 <- plyr::aaply(nf.tmp, 1, function(x) {x / exp(mean(log(x)))})
            normalizationFactors(dds1) <- x2
            dds1 <- DESeq(dds1)
            res1 <- results(dds1, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering)
            if (return.details) {
                list(res=res1,cm=cm,sampleGroups=sampleGroups)
            } else {
                res1
            }
        })
    }, mc.cores=n.cores)
    de.res
}






#' Do differential expression for each cell type in a conos object between the specified subsets of apps
#' @param conObj conos object
#' @param groups factor specifying cell types
#' @param sampleGroups a list of two character vector specifying the app groups to compare
#' @param cookscutoff cookscugoff for DESeq2
#' @param independentFiltering independentFiltering for DESeq2
#' @param n.cores number of cores
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param return.details return detals
getPerCellTypeDE_raw_list <- function(conObj, groups=NULL, sampleGroups=NULL, cooksCutoff = FALSE, reflevel = NULL,
                             independentFiltering = FALSE, n.cores=1,cluster.sep.chr = '+',return.details=TRUE) {
    ## Check arguments
    if ( is.null(groups) ) stop('groups must be specified');
    if ( is.null(sampleGroups) ) stop('sampleGroups must be specified')
    if ( class(sampleGroups) != 'list' ) stop('sampleGroups must be a list');
    if ( length(sampleGroups) != 2 ) stop('sampleGroups must be of length 2');
    if ( ! all(unlist(lapply(sampleGroups, function(x) class(x) == 'character'))) )
        stop('sampleGroups must be a list of character vectors');
    if ( ! all(unlist(lapply(sampleGroups, function(x) length(x) > 0))) )
        stop('sampleGroups entries must be on length greater or equal to 1')
    if ( ! all(unlist(lapply(sampleGroups, function(x) {all(x %in% names(conObj))}))) )
        stop('sampleGroups entries must be names of samples in the conos object')
    if ( is.null(reflevel) ) stop('reference level is not defined')
    ## todo: check samplegrousp are named
    if(is.null(names(sampleGroups))) stop('sampleGroups must be named')
    if(class(groups) != 'factor') stop('groups must be a factor')
    if(any(grepl(cluster.sep.chr, names(conObj),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any sample name')
    if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE))) 
        stop('cluster.sep.chr must not be part of any cluster name')
    ## Generate a summary dataset collapsing the cells of the same type in each sample
    ## and merging everything in one matrix
    samples.used <- unlist(sampleGroups)
    ## Generate an aggregated matrix
    raw.mats <- lapply(conObj[samples.used], function(p2) {
        p2
    })
    common.genes <- Reduce(intersect,lapply(raw.mats, colnames))
    raw.mats <- lapply(raw.mats, function(x) {x[,common.genes]})
    aggr2 <- lapply(raw.mats, function(x) {
        g1 <- groups[intersect(names(groups), rownames(x))]
        aggr <- Matrix.utils::aggregate.Matrix(x, g1)
        aggr
    })
    aggr2 <- lapply(names(aggr2), function(n) {
        x <- aggr2[[n]]
        rownames(x) <- paste0(n,cluster.sep.chr,rownames(aggr2[[n]]))
        x
    })
    aggr2 <- t(do.call(rbind, aggr2))
    rm(raw.mats); gc()
    ## For every cell type get differential expression results
    de.res <- mclapply(namedLevels(groups), function(l) {
        try({
            ## Get count matrix
            cm <- aggr2[,strpart(colnames(aggr2),cluster.sep.chr,2,fixed=TRUE) == l]
            ## Generate metadata
            meta <- data.frame(
                sample.id= colnames(cm),
                group= as.factor(unlist(lapply(colnames(cm), function(y) {
                    y <- strpart(y,cluster.sep.chr,1,fixed=TRUE)
                    names(sampleGroups)[unlist(lapply(sampleGroups,function(x) any(x %in% y)))]
                })))
            )
            meta$group <- relevel(meta$group, ref=reflevel)
            if (length(unique(as.character(meta$group))) < 2)
                stop('The cluster is not present in both conditions')
            dds1 <- DESeqDataSetFromMatrix(cm,meta,design=~group)
            dds1 <- DESeq(dds1)
            res1 <- results(dds1, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering)
            ##
            if(return.details) {
                list(res=res1, cm=cm, sampleGroups=sampleGroups)
            } else {
                res1
            }
        })
    }, mc.cores=n.cores)
    de.res
}








#' Do differential expression for each cell type in a conos object between the specified subsets of apps
#' @param conObj conos object
#' @param groups factor specifying cell types
#' @param sampleGroups a list of two character vector specifying the app groups to compare
#' @param cookscutoff cookscugoff for DESeq2
#' @param independentFiltering independentFiltering for DESeq2
#' @param n.cores number of cores
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param return.details return detals
getPerCellTypeDE <- function(conObj, groups=NULL, sampleGroups=NULL, cooksCutoff = FALSE, reflevel = NULL,
                             independentFiltering = FALSE, n.cores=1,cluster.sep.chr = '+',return.details=TRUE) {
    ## Check arguments
    if ( class(conObj) != 'Conos') stop('conObj must be a conos object')
    if ( is.null(groups) ) stop('groups must be specified');
    if ( is.null(sampleGroups) ) stop('sampleGroups must be specified')
    if ( class(sampleGroups) != 'list' ) stop('sampleGroups must be a list');
    if ( length(sampleGroups) != 2 ) stop('sampleGroups must be of length 2');
    if ( ! all(unlist(lapply(sampleGroups, function(x) class(x) == 'character'))) )
        stop('sampleGroups must be a list of character vectors');
    if ( ! all(unlist(lapply(sampleGroups, function(x) length(x) > 0))) )
        stop('sampleGroups entries must be on length greater or equal to 1')
    if ( ! all(unlist(lapply(sampleGroups, function(x) {all(x %in% names(conObj$samples))}))) )
        stop('sampleGroups entries must be names of samples in the conos object')
    if ( is.null(reflevel) ) stop('reference level is not defined')
    ## todo: check samplegrousp are named
    if(is.null(names(sampleGroups))) stop('sampleGroups must be named')
    if(class(groups) != 'factor') stop('groups must be a factor')
    if(any(grepl(cluster.sep.chr, names(conObj$samples),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any sample name')
    if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE))) 
        stop('cluster.sep.chr must not be part of any cluster name')
    ## Generate a summary dataset collapsing the cells of the same type in each sample
    ## and merging everything in one matrix
    samples.used <- unlist(sampleGroups)
    ## Generate an aggregated matrix
    raw.mats <- lapply(conObj$samples[samples.used], function(p2) {
        p2$misc$rawCounts
    })
    common.genes <- Reduce(intersect,lapply(raw.mats, colnames))
    raw.mats <- lapply(raw.mats, function(x) {x[,common.genes]})
    aggr2 <- lapply(raw.mats, function(x) {
        g1 <- groups[intersect(names(groups), rownames(x))]
        aggr <- Matrix.utils::aggregate.Matrix(x, g1)
        aggr
    })
    aggr2 <- lapply(names(aggr2), function(n) {
        x <- aggr2[[n]]
        rownames(x) <- paste0(n,cluster.sep.chr,rownames(aggr2[[n]]))
        x
    })
    aggr2 <- t(do.call(rbind, aggr2))
    rm(raw.mats); gc()
    ## For every cell type get differential expression results
    de.res <- mclapply(namedLevels(groups), function(l) {
        try({
            ## Get count matrix
            cm <- aggr2[,strpart(colnames(aggr2),cluster.sep.chr,2,fixed=TRUE) == l]
            ## Generate metadata
            meta <- data.frame(
                sample.id= colnames(cm),
                group= as.factor(unlist(lapply(colnames(cm), function(y) {
                    y <- strpart(y,cluster.sep.chr,1,fixed=TRUE)
                    names(sampleGroups)[unlist(lapply(sampleGroups,function(x) any(x %in% y)))]
                })))
            )
            meta$group <- relevel(meta$group, ref=reflevel)
            if (length(unique(as.character(meta$group))) < 2)
                stop('The cluster is not present in both conditions')
            dds1 <- DESeqDataSetFromMatrix(cm,meta,design=~group)
            dds1 <- DESeq(dds1)
            res1 <- results(dds1, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering)
            ##
            if(return.details) {
                list(res=res1, cm=cm, sampleGroups=sampleGroups)
            } else {
                res1
            }
        })
    }, mc.cores=n.cores)
    de.res
}

#' Obtain a correction vector for removing the constant effect between the same clusters of two different apps
#' @param conObj conos object
#' @param groups a vector specifying clusters
#' @param sampleGroups the groups of the samples
#' @param cooksCutoff cooksCutoff distance for DESeq2
#' @param independentFiltering independentFiltering for DESeq2
#' @param n.cores number of cores
#' @param cluster.sep.chr separator for cluster and sample name
#' @param return.details logical, if TRUE return internal sturcuters
#' @param de.init if specified reuses existing differential expression results
#' @param exclude.celltypes names of cell types to exclude from the generation of the vecotr
#' @param correction.method 'varianceweighted' or 'mean' specifies way to merge the fold changes from different cell types
getCorrectionVector <- function(conObj, groups=NULL, sampleGroups=NULL, cooksCutoff=FALSE, independentFiltering = FALSE, n.cores=1, cluster.sep.chr = '+', return.details=FALSE,de.init=NULL,exclude.celltypes=c(),correction.method='varianceweighted',reflevel=NULL) {
    ## conObj <- con; sampleGroups <- sampleGroupsTvsW; cooksCutoff<-F; independentFiltering <- F; n.cores <- 1; cluster.sep.chr <- '+'
    ## groups <- as.factor(jcl3.coarse)
    ## Check arguments
    if ( class(conObj) != 'Conos') stop('conObj must be a conos object')
    if ( is.null(groups) ) stop('groups must be specified');
    if ( is.null(sampleGroups) ) stop('sampleGroups must be specified')
    if ( class(sampleGroups) != 'list' ) stop('sampleGroups must be a list');
    if ( length(sampleGroups) != 2 ) stop('sampleGroups must be of length 2');
    if ( ! all(unlist(lapply(sampleGroups, function(x) class(x) == 'character'))) )
        stop('sampleGroups must be a list of character vectors');
    if ( ! all(unlist(lapply(sampleGroups, function(x) length(x) > 0))) )
        stop('sampleGroups entries must be on length greater or equal to 1')
    if ( ! all(unlist(lapply(sampleGroups, function(x) {all(x %in% names(conObj$samples))}))) )
        stop('sampleGroups entries must be names of samples in the conos object')
    ## todo: check samplegrousp are named
    if(is.null(names(sampleGroups))) stop('sampleGroups must be named')
    if(class(groups) != 'factor') stop('groups must be a factor')
    if(any(grepl(cluster.sep.chr, names(conObj$samples),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any sample name')
    if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE))) 
        stop('cluster.sep.chr must not be part of any cluster name')
    if(is.null(reflevel)) stop('reference level is not defined')
    ## Main function
    if(is.null(de.init)) {
        de.init <- getPerCellTypeDE(conObj, groups=groups, sampleGroups=sampleGroups,
                                    cooksCutoff=cooksCutoff, independentFiltering=independentFiltering,
                                    n.cores=n.cores, cluster.sep.chr=cluster.sep.chr,return.details=FALSE,
                                    reflevel = reflevel);
    }
    allfcs <- lapply(de.init, function(x) {
        if(!is.error(x)) {
            fc <- x$log2FoldChange
            names(fc) <- rownames(x)
            fc
        } else {
            NULL
        }
    })
    allfcs <- allfcs[!unlist(lapply(allfcs, is.null))]
    genes <- Reduce(intersect, lapply(allfcs, names))
    ## Matrix of fold changes
    fc.mat <- do.call(rbind, lapply(allfcs, function(x) {x[genes]}))
    fc.mat <- fc.mat[!rownames(fc.mat) %in% exclude.celltypes,]
    if (correction.method == 'mean') {
        correction <- apply(fc.mat, 2, mean, na.rm=TRUE)
    } else if (correction.method == 'varianceweighted') {
        mu <- apply(fc.mat, 2, mean, na.rm=TRUE)
        var <- apply(fc.mat, 2, function(x) {var(x,na.rm=TRUE)})
        weight <- 1 - pchisq(q=var,df=nrow(fc.mat)-1)
        correction <- mu * weight
    } else {
        error(paste0('unknown correction method: ', correction.method))
    }
    ## return
    if (!return.details) {
        correction
    } else {
        list(correction.vector=correction,de.init=de.init)
    }        
}

#' Do differential expression for each cell type in a conos object between the specified subsets of apps
#' applying the specified correction vector
#' @param conObj conos object
#' @param groups factor specifying cell types
#' @param sampleGroups a list of two character vector specifying the app groups to compare
#' @param cookscutoff cookscugoff for DESeq2
#' @param independentFiltering independentFiltering for DESeq2
#' @param n.cores number of cores
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param correction correction vector obtained from getCorrectionVector
getPerCellTypeDECorrected <- function(conObj, groups=NULL, sampleGroups=NULL, cooksCutoff = FALSE,
                             independentFiltering = FALSE, n.cores=1,cluster.sep.chr = '+',
                             correction=NULL, return.details=TRUE, reflevel=NULL) {
    ## Check arguments
    if ( is.null(correction) ) stop("Correction can't by null'")
    if ( class(conObj) != 'Conos') stop('conObj must be a conos object')
    if ( is.null(groups) ) stop('groups must be specified');
    if ( is.null(sampleGroups) ) stop('sampleGroups must be specified')
    if ( class(sampleGroups) != 'list' ) stop('sampleGroups must be a list');
    if ( length(sampleGroups) != 2 ) stop('sampleGroups must be of length 2');
    if ( ! all(unlist(lapply(sampleGroups, function(x) class(x) == 'character'))) )
        stop('sampleGroups must be a list of character vectors');
    if ( ! all(unlist(lapply(sampleGroups, function(x) length(x) > 0))) )
        stop('sampleGroups entries must be on length greater or equal to 1')
    if ( ! all(unlist(lapply(sampleGroups, function(x) {all(x %in% names(conObj$samples))}))) )
        stop('sampleGroups entries must be names of samples in the conos object')
    if (is.null(reflevel)) stop('reference level is not defined')
    ## todo: check samplegrousp are named
    if(is.null(names(sampleGroups))) stop('sampleGroups must be named')
    if(class(groups) != 'factor') stop('groups must be a factor')
    if(any(grepl(cluster.sep.chr, names(conObj$samples),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any sample name')
    if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE))) 
        stop('cluster.sep.chr must not be part of any cluster name')
    ## Generate a summary dataset collapsing the cells of the same type in each sample
    ## and merging everything in one matrix
    samples.used <- unlist(sampleGroups)
    ## Generate an aggregated matrix
    raw.mats <- lapply(conObj$samples[samples.used], function(p2) {
        p2$misc$rawCounts
    })
    common.genes <- Reduce(intersect,lapply(raw.mats, colnames))
    raw.mats <- lapply(raw.mats, function(x) {x[,common.genes]})
    aggr2 <- lapply(raw.mats, function(x) {
        g1 <- groups[intersect(names(groups), rownames(x))]
        aggr <- Matrix.utils::aggregate.Matrix(x, g1)
        aggr
    })
    aggr2 <- lapply(names(aggr2), function(n) {
        x <- aggr2[[n]]
        rownames(x) <- paste0(n,cluster.sep.chr,rownames(aggr2[[n]]))
        x
    })
    aggr2 <- t(do.call(rbind, aggr2))
    rm(raw.mats); gc()
    ## For every cell type get differential expression results
    de.res <- mclapply(namedLevels(groups), function(l) {
        try({
            ## Get count matrix
            cm <- aggr2[,strpart(colnames(aggr2),cluster.sep.chr,2,fixed=TRUE) == l]
            ## Generate metadata
            meta <- data.frame(
                sample.id= colnames(cm),
                group= as.factor(unlist(lapply(colnames(cm), function(y) {
                    y <- strpart(y,cluster.sep.chr,1,fixed=TRUE)
                    names(sampleGroups)[unlist(lapply(sampleGroups,function(x) any(x %in% y)))]
                })))
            )
            meta$group <- relevel(meta$group, ref=reflevel)
            if (length(unique(as.character(meta$group))) < 2)
                stop('The cluster is not present in both conditions')
            library(DESeq2)
            dds1 <- DESeqDataSetFromMatrix(cm,meta,design=~group)
            dds1 <- estimateSizeFactors(dds1)
            sf <- sizeFactors(dds1)
            if(!(all(rownames(cm) %in% names(correction)) & all(names(correction) %in% rownames(cm))))
                stop('incompatible matrices')
            nf.tmp <- matrix(rep(sf, nrow(cm)),nrow=nrow(cm),byrow=TRUE)
            rownames(nf.tmp) <- rownames(cm);
            colnames(nf.tmp) <- colnames(cm)
            gene.scale.factors <- 2^(correction[rownames(nf.tmp)])
            baselevel <- levels(colData(dds1)$group)[1]
            x <- do.call(cbind, lapply(colData(dds1)$group, function(x) {
                if (x == baselevel) {
                    rep(1, length(gene.scale.factors))
                } else {
                    gene.scale.factors
                }
            }))
            rownames(x) <- rownames(nf.tmp);
            colnames(x) <- colnames(nf.tmp)
            nf.tmp <- nf.tmp * x
            x2 <- plyr::aaply(nf.tmp, 1, function(x) {x / exp(mean(log(x)))})
            normalizationFactors(dds1) <- x2
            dds1 <- DESeq(dds1)
            res1 <- results(dds1, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering)
            if (return.details) {
                list(res=res1,cm=cm,sampleGroups=sampleGroups)
            } else {
                res1
            }
        })
    }, mc.cores=n.cores)
    de.res
}

saveDEasCSV <- function(de.results=NULL,saveprefix=NULL,gene.metadata=NULL) {
    if(is.null(de.results)) stop('de.results has not been specified')
    if(is.null(saveprefix)) stop('saveprefix has not bee specified')
    ## find errors
    n.error <- sum(unlist(lapply(de.results,is.error)))
    if(n.error > 0)
        cat("Warning: ", n.error, " of ", length(de.results), ' results have returned an error; ignoring...\n')
    de.results <- de.results[!unlist(lapply(de.results,is.error))]
    ## 
    x <- lapply(namedNames(de.results), function(ncc) {
        res.celltype <- de.results[[ncc]]
        res.table <- as.data.frame(res.celltype$res)
        ## append gene names
        res.table$gene <- rownames(res.table)
        ## append singificance
        res.table$significant <- res.table$padj < 0.05
        res.table$log2FoldChange[is.na(res.table$log2FoldChange)] <- 0 
        ## Append Z scores and rowid
        res.table$Z <- qnorm(1 - (res.table$pval/2))
        res.table$Z[is.na(res.table$Z)] <- 0
        res.table$Za <- qnorm(1 - (res.table$padj/2))
        res.table$Za[is.na(res.table$Za)] <- 0
        res.table$Z <- res.table$Z  * sign(res.table$log2FoldChange)
        res.table$Za <- res.table$Za  * sign(res.table$log2FoldChange)
        ## match order to metadata table
        mo <- match(as.character(gene.metadata$geneid),as.character(res.table$gene))
        ## drop gene id column
        keep.cols <- colnames(gene.metadata)[colnames(gene.metadata) != 'geneid']
        names(keep.cols) <- keep.cols
        res.table <- cbind(res.table, gene.metadata[mo,keep.cols,drop=FALSE])
        file <- paste0(saveprefix,make.names(ncc),'.csv')
        write.table(x=res.table,file=file)
        res.table
    })
    invisible(x)
}

saveDEasJSON <- function(de.results = NULL, saveprefix = NULL, gene.metadata = NULL) {
    ## ### DEVEL
    ## de.results <- all.percl.TvsW
    ## saveprefix <- 'json/'
    ## rm(de.results, saveprefix)
    ## ##
    ## Check input
    if(is.null(de.results)) stop('de.results have not been specified')
    if(is.null(saveprefix)) stop('saveprefix has not been specified')
    ## Find de instances that didn't work (usually because cell type is absent from one or more sample types)
    n.error <- sum(unlist(lapply(de.results, is.error)))
    if(n.error > 0)
        cat("Warning: ", n.error,' of ', length(de.results) ,' results have returned an error; ignoring...\n')
    ## get the de results that worked
    de.results <- de.results[!unlist(lapply(de.results, is.error))]
    ## Generate structure and save JSON
    lapply(namedNames(de.results), function(ncc) {
        res.celltype <- de.results[[ncc]]
        ## Get results table as df
        res.table <- as.data.frame(res.celltype$res)
        ## append gene names
        res.table$gene <- rownames(res.table)
        ## append singificance
        res.table$significant <- res.table$padj < 0.05
        res.table$log2FoldChange[is.na(res.table$log2FoldChange)] <- 0 
        ## Append Z scores and rowid
        res.table$Z <- qnorm(1 - (res.table$pval/2))
        res.table$Z[is.na(res.table$Z)] <- 0
        res.table$Za <- qnorm(1 - (res.table$padj/2))
        res.table$Za[is.na(res.table$Za)] <- 0
        res.table$Z <- res.table$Z  * sign(res.table$log2FoldChange)
        res.table$Za <- res.table$Za  * sign(res.table$log2FoldChange)
        res.table$rowid <- 1:nrow(res.table)
        ## match order to metadata table
        mo <- match(as.character(gene.metadata$geneid),as.character(res.table$gene))
        ## drop gene id column
        keep.cols <- colnames(gene.metadata)[colnames(gene.metadata) != 'geneid']
        names(keep.cols) <- keep.cols
        res.table <- cbind(res.table, gene.metadata[mo,keep.cols,drop=FALSE])
        ## get names of all the genes
        all.genes <- rownames(res.table)
        ## Get the count matrix
        cm <-res.celltype$cm
        ## remove the cell type suffix
        ## TODO make '+' a parameter
        colnames(cm) <- strpart(colnames(cm),'+',1,fixed=TRUE)
        ## ilev entry (submatrices of cps)
        ilev <- lapply(res.celltype$sampleGroups, function(sg) {
            ## In certain cases columns may be missing,skip
            sg <- sg[sg %in% colnames(cm)]
            ## keep only cols of interest
            cm.tmp <- cm[,sg]
            ## convert to matrix
            cm.tmp <- as.matrix(cm.tmp)
            rownames(cm.tmp) <- rownames(cm)
            ## calculate cpm
            cpm <- sweep(cm.tmp, 2, apply(cm.tmp,2, sum), FUN='/')
            cpm <- log10(cpm * 1e6 + 1)
            ## 
            snames1 <- colnames(cpm)
            ## Put genes in order
            cpm <- cpm[all.genes,]
            colnames(cpm) <- NULL;
            rownames(cpm) <- NULL;
            ## return
            list(snames=snames1, val=as.matrix(cpm))
        })
        ## snames entry (samplenames)
        snames <- names(res.celltype$sampleGroups)
        ## convert to json
        tojson <- list(
            res = res.table,
            genes = all.genes,
            ilev = ilev,
            snames = snames
        )
        y <- jsonlite::toJSON(tojson)
        ## File to save to 
        file <- paste0(saveprefix,make.names(ncc),'.json')
        ## create the json file
        write(y,file)
        NULL
    })
    invisible(NULL)
}









# input fc or normalized Z score 




library('liger')
library(qusage)



runGSEA=function(fc){
  msigdb.h <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c2.cp.kegg.v6.2.symbols.gmt')
  msigdb.c2 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c2.all.v6.2.symbols.gmt')
  msigdb.c5 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c5.all.v6.2.symbols.gmt')
  msigdb.c6 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c6.all.v6.2.symbols.gmt')
  msigdb.c7 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c7.all.v6.2.symbols.gmt')
  sign.sets <- list(kegg=msigdb.h,curatedGeneSet=msigdb.c2,GO=msigdb.c5,onocogene=msigdb.c6,immunologic=msigdb.c7)
  fc <- fc[!is.na(fc)]
  
  fc=fc[order(fc,decreasing=T)]
  r=lapply(namedNames(sign.sets), function(n) {
    cat(paste0('\tProcessing set ',n,'\n'))
    signset <- sign.sets[[n]]
    gsea.res <- iterative.bulk.gsea(fc, set.list=signset, mc.cores=15)
    gsea.res[order(gsea.res$q.val),]
  })
  return(r)
}




runDAVID2<-function(genelist,appname){
  genelist=genelist[!is.na(genelist)]
  ENTREZID=unlist(mget(genelist, org.Hs.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID=ENTREZID[!is.na(ENTREZID)]
  
  tmpF=paste(appname,'.ENTREZID',sep='')
  fout=paste(appname,'.DAVID.txt',sep='')
  write.table(ENTREZID,tmpF,sep='\t',row.names=F,col.names=F,quote=F)
  cmd=paste("python /d0/home/meisl/Workplace/Archive/DAVID/DAVID_annotation_chart.py --genelist=", tmpF ,' --email="ma.tongji@gmail.com" --out=',fout,sep='')  
  system(cmd)
  
  cmd2=paste('Rscript /d0/home/meisl/Workplace/Archive/DAVID/DAVID_to_figure_vison.r ',fout,' --termNumber=18  -o ',appname,sep='')
  # system(cmd2)
  
  print(cmd2)
  ENTREZID2=names(ENTREZID)
  names(ENTREZID2)=ENTREZID
  
  tmp=read.csv(fout,sep='\t',header=T)
  tmp$symb=apply(tmp,1,function(x) paste(ENTREZID2[strsplit(as.character(unlist(x['Genes'][1])),', ')[[1]]],collapse=','))
  
  write.table(tmp,fout,sep='\t',row.names=F,col.names=T,quote=F)
  return(tmp)
}



DAVID_topG=function(fc,appname){
  fc <- fc[!is.na(fc)]
  fc=fc[order(fc,decreasing=T)]
  glis=list()
  for ( num in c(200,300,400,500)){
    genelist=names(fc[1:num])
    tmpn=paste(appname,'_',num,sep='')
    glis[[tmpn]]=runDAVID2(genelist,tmpn)
  }
  return(glis)
}



#categrey in c('KEGG_PATHWAY','GOTERM_BP_DIRECT')

DAVID_merge=function(appname,categrey){
  
  fo=paste(appname,'.DAVID.rds',sep='')
  rlis=readRDS(fo)
  
  plist=list()
  splist=list()
  
  
  for (i in names(rlis)){
    res=rlis[[i]]
    res=res[res[,1]==categrey,]
    pval=res$Pvalue
    names(pval)=as.character(res[,2])
    
    plist[[i]]=pval
    num=15
    if(length(pval)<15){
      num=length(pval)
    }
    splist[[i]]=names(pval)[1:num]
    
  }

  
  Pmatrix = listTomatrix(splist)
  Pmatrix =Pmatrix[rowSums(Pmatrix)>2,]
  
  for ( i in names(splist)){
    inter=intersect(rownames(Pmatrix),names(plist[[i]]))
    Pmatrix[inter,i]=plist[[i]][inter]
  }
  
  tmp=Pmatrix
  
  library(pheatmap)
  library(gplots)
  
  tmp[is.na(tmp)]=0
  tmp=(-log(tmp))
  tmp[is.na(tmp)]=0
  tmp[is.nan(tmp)]=0
  tmp[!is.finite(tmp)]=0
  
  
  rgb.palette <- colorRampPalette(c("white","red"), space = "rgb" )
  mi=min(tmp)
  ma=max(tmp)
  
  pdf=paste(appname,'.',categrey,'.pdf',sep='')
  a=pheatmap(tmp,filename=pdf,color=rgb.palette(100),scale='none',border_color='NA',cluster_rows=T,cluster_cols=T,
             show_colnames=T,fontsize_col=6,fontsize_row=8,height=0.6*nrow(tmp),
             breaks = c(mi,seq(2,ma,length.out = 99)))
  
  return(rowMeans(tmp))
}



GO_Data=function(fc,appname){
  
  allR=runGSEA(fc)
  fo=paste(appname,'.gsea.new.rds',sep='')
  saveRDS(allR,fo)
  
  rlis=DAVID_topG(fc,appname)

  fo=paste(appname,'.DAVID.rds',sep='')
  saveRDS(rlis,fo)
  
}





listTomatrix_Value=function(res){
  # list to a matrix,  row is union genes
  
  unionPos=NULL
  for ( i in names(res)){
    unionPos=union(unionPos,names(res[[i]]))
  }
  #  mutation matix
  lc=length(res)
  lr=length(unionPos)
  stat=matrix(rep(0,lc*lr),lr,lc)
  for( i in seq(lc)){
    #\tindex=match(unionPos,tmp$pos)
    #\tstat[!is.na(index),i]=1
    index=match(names(res[[i]]),unionPos)
    stat[index,i]=res[[i]]
    #\tstat[index,i]=1
  }
  colnames(stat)=names(res)
  rownames(stat)=unionPos
  return(stat)
}





prepare_rawList=function(bigM2,group4,atype,allano){
  cell_ano=group4[group4==atype]
  cname=names(cell_ano)
  anoSample=allano[['sample']][cname]
 
  index=grepl('BMET11-',names(anoSample))
  anoSample=anoSample[!index]
  index=grepl('BMET10-',names(anoSample))
  anoSample=anoSample[!index]
  
  #index=grepl('BMET1-',names(anoSample))
  #anoSample=anoSample[!index] 
  gname=names(anoSample)

  t.bigM2=bigM2[,gname]
  
  t.bigM2=t.bigM2[rowSums(t.bigM2)>15,]
  tab=table(allano[['sample']][gname])
  nnames=names(tab[tab>15])
  
  raw.mats=list()
  
  n.anoSample=allano[['sample']][gname]
  p2.lis=list()
  for (n in nnames){
    cd=t.bigM2[,n.anoSample==n]
    n1=n
    raw.mats[[n1]]=cd
  }
	nn=unlist(lapply(raw.mats, function(x) ncol(x)))
	print(atype)
	print(nn)
	return(raw.mats)
  
}





RunDE_cellType=function(raw.mats2,jcl3.coarse,appname){

  source('/home/meisl/bin/FunctionLib/Lib/pagodaLib.r')
  
  allano <- readRDS('/d0-mendel/home/meisl/Workplace/BMME/a.data/jcl3.coarse.all.rds')
  allclean <- readRDS('/d0-mendel/home/meisl/Workplace/BMME/a.data/jcl3.cleanup.all.rds')
  
  tmp=readRDS('/d0-mendel/home/meisl/Workplace/BMME/a.data/Selected_Joint_embdding/anCell_12_28.rds')
  allano[['cell']]=tmp
  
  
  
  cytokine=readRDS('/home/meisl/bin/data/cytokine.rds')
  surface=readRDS('/home/meisl/bin/data/membrane.rds')
  ano=readRDS('/home/meisl/bin/data/gene.annotation.rds')
  
  
  ## run differential expressed genes 
  membrane <- readRDS('/d0-mendel/home/barkasn/work/extracellularProteins/output/hs.membrane.hgnc.rds')
  
  source('/home/meisl/bin/FunctionLib/Lib/de.functions.R')
  membrane=cytokine
  
  
  
  ## membrane information
  getMeta <- function(res) {
    all.genes <- unique(unlist(lapply(res, function(x) {
      if(!is.error(x)){
        rownames(as.data.frame(x$res))
      } else {
        NULL
      }
    })))
    gene.metadata <- data.frame(
      geneid=all.genes,
      membrane=all.genes %in% membrane
    )
    gene.metadata
  }
  
  
  nn=apply(data.frame(names(raw.mats2)),1,function(x) strsplit(x,'-')[[1]][2]) %>% table()
  fname=names(nn[nn>1])
  
  
  ## 
  allres=list()
  Vec=readRDS('/home/meisl/Workplace/BMME/b.newAno/diffGene/rds/correctedVec.rds')
  
  for ( i in seq(length(fname))){
    for (j in seq( length(fname))){
      if (fname[i]!=fname[j]){
        print('##')
        print(fname[i])
        print(fname[j])
        
        sampleGroups <- list(
          treat=names(raw.mats2)[grepl(fname[i],names(raw.mats2))],
          control=names(raw.mats2)[grepl(fname[j],names(raw.mats2))]
        )
        
        n1=paste(fname[i],'_vs_',fname[j],sep='')      
        
        
        if (grepl('Tumor',n1)){
          print('adjust')
          fc.TvsW=Vec[[n1]]
          gs=intersect(names(fc.TvsW),colnames(raw.mats2[[1]]))
          raw.mats3=lapply(raw.mats2,function(x) x[,gs])
          
          all.percl.TvsW.nocorr <- getPerCellTypeDECorrected_raw_list(conObj=raw.mats3,
                                                                      groups=jcl3.coarse,
                                                                      sampleGroups=sampleGroups,
                                                                      n.cores=32,correction=fc.TvsW[gs],reflevel='control')
          
        }else{
          ## Do differential expression without correction
          all.percl.TvsW.nocorr <- getPerCellTypeDE_raw_list(conObj=raw.mats2,
                                                             groups=jcl3.coarse,
                                                             sampleGroups=sampleGroups,
                                                             reflevel='control',
                                                             n.cores=22)
        }
        
        allres[[n1]]=all.percl.TvsW.nocorr
        #      fo=paste(fname[i],'_vs_',fname[j],'.rds',sep='')
        #      saveRDS(all.percl.TvsW.nocorr,fo)
      }
    }
    
    
  }
  
  fo=paste(appname,'.DESeq.with.BMET1.rds',sep='')
  saveRDS(allres,fo)
  
}


GETnomrlizedCount=function(raw.mats,appname){

  count1=do.call(cbind,lapply(raw.mats,function(x) rowSums(x)))
  
  tab1=apply(data.frame(colnames(count1)),1,function(x) strsplit(x,'-')[[1]][2])
  table(tab1)
  
  cm=count1
  meta <- data.frame(
    sample.id= colnames(cm),
    group= tab1
  )
  dds1 <- DESeqDataSetFromMatrix(cm,meta,design=~group)
  dds1 <- DESeq(dds1)
  
  dat_DEseq=counts(dds1, normalized=TRUE)
  
  cm.tmp=count1
  cm.tmp <- as.matrix(cm.tmp)
  rownames(cm.tmp) <- rownames(cm)
  ## calculate cpm
  cpm <- sweep(cm.tmp, 2, apply(cm.tmp,2, sum), FUN='/')
  dat_cpm <- log10(cpm * 1e6 + 1)

  r=list('dat_DEseq'=dat_DEseq,'dat_cpm'=dat_cpm,'raw'=raw.mats)
  
  fo=paste(appname,'.normlized.count.rds',sep='')
  saveRDS(r,fo)
  
  return(r)
}




