simulate.data <- function(G=5, N=100, M=1000, initmean=0, initvar=10, upreg=10, upregvar=10, ng=100, seed=0, plot=TRUE) {
    set.seed(seed)
    mat <- matrix(rnorm(N*M*G, initmean, initvar), M, N*G)
    rownames(mat) <- paste0('gene', 1:M)
    colnames(mat) <- paste0('cell', 1:(N*G))
    group <- factor(sapply(1:G, function(x) rep(paste0('group', x), N)))
    names(group) <- colnames(mat)

    diff <- lapply(1:G, function(x) {
        diff <- rownames(mat)[(((x-1)*ng)+1):(((x-1)*ng)+ng)]
        mat[diff, group==paste0('group', x)] <<- mat[diff, group==paste0('group', x)] + rnorm(ng, upreg, upregvar)
        return(diff)
    })
    names(diff) <- paste0('group', 1:G)

    diff2 <- lapply(2:(G-1), function(x) {
        y <- x+G
        diff <- rownames(mat)[(((y-1)*ng)+1):(((y-1)*ng)+ng)]
        mat[diff, group %in% paste0("group", 1:x)] <<- mat[diff, group %in% paste0("group", 1:x)] + rnorm(ng, upreg, upregvar)
        return(diff)
    })

    mat[mat<0] <- 0
    mat <- t(t(mat)/colSums(mat))
    mat <- log10(mat*1e6+1)

    if(plot) {
        heatmap(mat, Rowv=NA, Colv=NA, col=colorRampPalette(c('blue', 'white', 'red'))(100), scale="none", ColSideColors=rainbow(G)[group], labCol=FALSE, labRow=FALSE)
    }

    return(mat)
}

fast.pca <- function(m,nPcs=2,tol=1e-10,scale=F,center=F,transpose=F) {
    require(irlba)
    ## note transpose is meant to speed up calculations when neither scaling nor centering is required
    if(transpose) {
        if(center) { m <- m-Matrix::rowMeans(m)}; if(scale) { m <- m/sqrt(Matrix::rowSums(m*m)); }
        a <- irlba(tcrossprod(m)/(ncol(m)-1), nu=0, nv=nPcs,tol=tol);
        a$l <- t(t(a$v) %*% m)
    } else {
        if(scale||center) { m <- scale(m,scale=scale,center=center) }
        a <- irlba(crossprod(m)/(nrow(m)-1), nu=0, nv=nPcs,tol=tol);
        a$l <- m %*% a$v
    }
    a
}

get.knn.membership <- function(mat, k, method=igraph::cluster_walktrap, PCA=FALSE, nPC=100, verbose=TRUE) {
    ## approximate nearest neighbor
    if(PCA) {
        ## use PCA to reduce dimensionality
        if(verbose) {
            print(paste0("identifying ", nPC, " principal components ..."))
        }

        pcs <- fast.pca(t(mat), nPC)
        m <- t(pcs$l)
        colnames(m) <- colnames(mat)
        rownames(m) <- paste0('PC', seq_len(nPC))
        mat <- m
    }

    require(RANN)
    if(verbose) {
        print("finding approximate nearest neighbors ...")
    }
    ##require(nabor)
    knn <- RANN::nn2(t(mat), k=k)[[1]]
    ##knn <- nabor::knn(t(mat), k=k)

    ## convert to adjacency matrix
    adj <- matrix(0, ncol(mat), ncol(mat))
    rownames(adj) <- colnames(adj) <- colnames(mat)
    invisible(lapply(seq_len(ncol(mat)), function(i) {
        adj[i,colnames(mat)[knn[i,]]] <<- 1
    }))

    ## convert to graph for clustering
    if(verbose) {
        print("calculating clustering by walktrap ...")
    }
    require(igraph)
    g <- igraph::graph.adjacency(adj, mode="undirected")
    g <- simplify(g)
    km <- method(g)
    if(verbose) {
        print(paste0('graph modularity: ', modularity(km)))
    }

    ## community membership
    com <- km$membership
    names(com) <- km$names
    if(verbose) {
        print("identifying cluster membership ...")
        print(table(com))
    }
    return(com)
}

model.lda <- function(mat, com, nfeatures=min(1000, nrow(mat)), verbose=TRUE, retest=TRUE) {
    ## filter
    if(nfeatures < nrow(mat)) {
        vi <- sample(seq_len(nrow(mat)), nfeatures)
        mat <- mat[vi,]
    }
    ## make data frame
    df <- data.frame(celltype=com, as.matrix(t(mat)))

    if(verbose) {
        print("calculating LDA ...")
    }
    require(MASS)
    model <- lda(celltype ~ ., data=df)

    if(retest) {
        ## predict our data based on model
        model.output <- predict(model, df)
        if(verbose) {
            print("LDA prediction accuracy ...")
            print(table(model.output$class==com))
        }
    }

    return(model)
}

tsne.lda <- function(mat, model, com, perplexity=30, ncores=10, verbose=TRUE, plot=TRUE, ...) {
    if(verbose) {
        print(paste0("Running Rtsne multicore with perplexity ", perplexity, " on ", ncores, " cores"))
    }

    ## compute LDs
    #reduction <- t(mat[rownames(model$scaling),]) %*% model$scaling
    reduction <- predict(model, data.frame(t(as.matrix(mat))))$x

    ## tSNE
    require(Rtsne.multicore)
    d <- dist(reduction)
    emb <- Rtsne.multicore::Rtsne.multicore(d, is_distance=TRUE, perplexity=perplexity, num_threads=ncores, verbose=verbose)$Y
    rownames(emb) <- labels(d)

    if(plot) {
        par(mfrow=c(1,1))
        plotEmbedding(emb, groups=com, ...)
    }

    return(emb)
}

plotEmbedding <- function(emb, clusterType=NULL, groups=NULL, colors=NULL, do.par=T, cex=0.6, alpha=0.4, gradientPalette=NULL, zlim=NULL, s=1, v=0.8, min.group.size=1, show.legend=FALSE, mark.clusters=FALSE, mark.cluster.cex=2, shuffle.colors=F, legend.x='topright', gradient.range.quantile=0.95, quiet=F, unclassified.cell.color='gray70', group.level.colors=NULL, ...) {
    factor.mapping=FALSE;
    if(is.null(colors) && is.null(groups)) {
        ## look up the clustering based on a specified type
        if(is.null(clusterType)) {
            ## take the first one
            groups <- clusters[[type]][[1]]
        } else {
            groups <- clusters[[type]][[clusterType]]
            if(is.null(groups)) { stop("clustering ",clusterType," for type ", type," doesn't exist")}
        }

        groups <- as.factor(groups[rownames(emb)]);
        if(min.group.size>1) { groups[groups %in% levels(groups)[unlist(tapply(groups,groups,length))<min.group.size]] <- NA; groups <- as.factor(groups); }
        factor.colors <- fac2col(groups,s=s,v=v,shuffle=shuffle.colors,min.group.size=min.group.size,level.colors=group.level.colors,return.details=T)
        cols <- factor.colors$colors[rownames(emb)]
        factor.mapping=TRUE;
    } else {
        if(!is.null(colors)) {
            ## use clusters information
            if(!all(rownames(emb) %in% names(colors))) { warning("provided cluster vector doesn't list colors for all of the cells; unmatched cells will be shown in gray. ")}
            if(all(areColors(colors))) {
                if(!quiet) cat("using supplied colors as is\n")
                cols <- colors[match(rownames(emb),names(colors))]; cols[is.na(cols)] <- unclassified.cell.color;
            } else {
                if(is.numeric(colors)) { # treat as a gradient
                    if(!quiet) cat("treating colors as a gradient")
                    if(is.null(gradientPalette)) { # set up default gradients
                        if(all(sign(colors)>=0)) {
                            gradientPalette <- colorRampPalette(c('gray80','red'), space = "Lab")(1024)
                        } else {
                            gradientPalette <- colorRampPalette(c("blue", "grey70", "red"), space = "Lab")(1024)
                        }
                    }
                    if(is.null(zlim)) { # set up value limits
                        if(all(sign(colors)>=0)) {
                            zlim <- as.numeric(quantile(colors,p=c(1-gradient.range.quantile,gradient.range.quantile)))
                            if(diff(zlim)==0) {
                                zlim <- as.numeric(range(colors))
                            }
                        } else {
                            zlim <- c(-1,1)*as.numeric(quantile(abs(colors),p=gradient.range.quantile))
                            if(diff(zlim)==0) {
                                zlim <- c(-1,1)*as.numeric(max(abs(colors)))
                            }
                        }
                    }
                                        # restrict the values
                    colors[colors<zlim[1]] <- zlim[1]; colors[colors>zlim[2]] <- zlim[2];

                    if(!quiet) cat(' with zlim:',zlim,'\n')
                    colors <- (colors-zlim[1])/(zlim[2]-zlim[1])
                    cols <- gradientPalette[colors[match(rownames(emb),names(colors))]*(length(gradientPalette)-1)+1]
                } else {
                    stop("colors argument must be a cell-named vector of either character colors or numeric values to be mapped to a gradient")
                }
            }
        } else {
            if(!is.null(groups)) {
                if(min.group.size>1) { groups[groups %in% levels(groups)[unlist(tapply(groups,groups,length))<min.group.size]] <- NA; groups <- droplevels(groups); }
                groups <- as.factor(groups)[rownames(emb)]
                if(!quiet) cat("using provided groups as a factor\n")
                factor.mapping=TRUE;
                ## set up a rainbow color on the factor
                factor.colors <- fac2col(groups,s=s,v=v,shuffle=shuffle.colors,min.group.size=min.group.size,unclassified.cell.color=unclassified.cell.color,level.colors=group.level.colors,return.details=T)
                cols <- factor.colors$colors;
            }
        }
        names(cols) <- rownames(emb)
    }

    if(do.par) {
        par(mar = c(0.5,0.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
    }
    plot(emb,col=adjustcolor(cols,alpha=alpha),cex=cex,pch=19,axes=F, ...); box();
    if(mark.clusters) {
        if(!is.null(groups)) {
            cent.pos <- do.call(rbind,tapply(1:nrow(emb),groups,function(ii) apply(emb[ii,,drop=F],2,median)))
            cent.pos <- na.omit(cent.pos);
            text(cent.pos[,1],cent.pos[,2],labels=rownames(cent.pos),cex=mark.cluster.cex)
        }
    }
    if(show.legend) {
        if(factor.mapping) {
            legend(x=legend.x,pch=rep(19,length(levels(groups))),bty='n',col=factor.colors$palette,legend=names(factor.colors$palette))
        }
    }
}
## a utility function to translate factor into colors
fac2col <- function(x,s=1,v=1,shuffle=FALSE,min.group.size=1,return.details=F,unclassified.cell.color='gray50',level.colors=NULL) {
    x <- as.factor(x);
    if(min.group.size>1) {
        x <- factor(x,exclude=levels(x)[unlist(tapply(rep(1,length(x)),x,length))<min.group.size])
        x <- droplevels(x)
    }
    if(is.null(level.colors)) {
        col <- rainbow(length(levels(x)),s=s,v=v);
    } else {
        col <- level.colors[1:length(levels(x))];
    }
    names(col) <- levels(x);

    if(shuffle) col <- sample(col);

    y <- col[as.integer(x)]; names(y) <- names(x);
    y[is.na(y)] <- unclassified.cell.color;
    if(return.details) {
        return(list(colors=y,palette=col))
    } else {
        return(y);
    }
}
# quick utility to check if given character vector is colors
# thanks, Josh O'Brien: http://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
areColors <- function(x) {
      is.character(x) & sapply(x, function(X) {tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)})
}







##################### Runners

## test functions using simulated data
test.sim <- function() {
    mat <- simulate.data()
    com <- get.knn.membership(mat, k=10, PCA=TRUE, nPC=10)
    model <- model.lda(mat, com)
    emb <- tsne.lda(mat, model, com, perplexity=30, show.legend=TRUE, mark.clusters=TRUE)
}

## test functions using 10X data
test.10x <- function() {
    require(Matrix)
    path <- '../data/'
    files <- list.files(path)
    sn <- 'frozen'
    names <- files <- files[grepl(sn, files)]
    names
    files <- paste0(path, files)
    names <- gsub('_cd.RData', '', names)
    names

    cds <- lapply(1:length(files), function(i) {
        f = files[i]
        n = names[i]
        load(f)

        summary <- log10(Matrix::colSums(cd))
        hist(summary, main=n)

        colnames(cd) <- paste0(n, '_', colnames(cd))
        as(cd, 'dgCMatrix')

    })
    names(cds) <- names
    length(cds)

    sample <- unlist(lapply(1:length(cds), function(i) {
        cd <- cds[[i]]
        n <- names(cds)[i]
        sample <- rep(n, ncol(cd))
        names(sample) <- colnames(cd)
        return(sample)
    }))
    sample <- factor(sample)
    length(sample)
    table(sample)


    ##################### Filter
    k <- 30
    vi <- Reduce(intersect, lapply(cds, function(cd) {
                                vi <- rowSums(cd)>k
                                rownames(cd)[vi]
                            }))
    length(vi)
    head(vi)

    cds.filter <- lapply(cds, function(cd) cd[vi,])
    names(cds.filter)
    lapply(cds.filter, dim)


    ##################### All together
    comb.mat <- do.call(cbind, cds.filter)
    comb.mat <- t(t(comb.mat)/colSums(comb.mat))
    comb.mat <- log10(comb.mat*1e6+1)
    vg <- names(sort(apply(comb.mat, 1, var), decreasing=TRUE))[1:500]
    pcs <- fast.pca(t(comb.mat), 100)
    pca.model <- pcs$v
    rownames(pca.model) <- rownames(comb.mat)
    colnames(pca.model) <- paste0('PC', seq_len(100))

    m <- t(pcs$l)
    colnames(m) <- colnames(comb.mat)
    rownames(m) <- paste0('PC', seq_len(100))
    pca.mat <- m
    comb.com <- get.knn.membership(pca.mat, k=150)
    lda.model <- model.lda(comb.mat, comb.com, nfeatures=1000)

    lda.emb <- tsne.lda(comb.mat, lda.model, comb.com, perplexity=150, show.legend=TRUE, mark.clusters=TRUE)

    d <- dist(pca.mat)
    pca.emb <- Rtsne.multicore::Rtsne.multicore(d, is_distance=TRUE, perplexity=150, num_threads=10)$Y
    rownames(pca.emb) <- labels(d)

    par(mfrow=c(2,2))
    comb.emb <- lda.emb
    plotEmbedding(comb.emb, groups=sample, show.legend=TRUE, mark.clusters=TRUE, alpha=0.1)
    plotEmbedding(comb.emb, groups=comb.com, show.legend=TRUE, mark.clusters=TRUE)
    comb.emb <- pca.emb
    plotEmbedding(comb.emb, groups=sample, show.legend=TRUE, mark.clusters=TRUE, alpha=0.1)
    plotEmbedding(comb.emb, groups=comb.com, show.legend=TRUE, mark.clusters=TRUE)



    ##################### PCA model



    ##################### LDA models
    models <- lapply(cds.filter, function(mat) {
        ##mat <- cds.filter[[1]]
        mat <- t(t(mat)/colSums(mat))
        mat <- log10(mat*1e6+1)
        com <- get.knn.membership(mat, k=50, PCA=TRUE)
        model <- model.lda(mat, com, nfeatures=1000)
        emb <- tsne.lda(mat, model, com, perplexity=50, show.legend=TRUE, mark.clusters=TRUE)
        return(model)
    })

    ################### Apply individual models
    reduction <- do.call(cbind, lapply(models, function(model) {
                                    predict(model, data.frame(t(as.matrix(comb.mat))))$x
                                }))
    d <- dist(reduction)
    perplexity = 150
    ncores = 10
    verbose=TRUE
    ind.emb <- Rtsne.multicore::Rtsne.multicore(d, is_distance=TRUE, perplexity=perplexity, num_threads=ncores, verbose=verbose)$Y
    rownames(ind.emb) <- labels(d)

    plotEmbedding(ind.emb, groups=sample, show.legend=TRUE, mark.clusters=TRUE, alpha=0.1)
    plotEmbedding(ind.emb, groups=comb.com, show.legend=TRUE, mark.clusters=TRUE)


    ##################### Markers
    gene.sets <- list()
    gene.sets[['Mono/Macro/Neutro']] <- c('CD14', 'ITGAM', 'FCGRA1', 'FCGR3A', 'IL8', 'CD68', 'IL1B', 'AIF1', 'CD74', 'S100A8', 'S100A9', 'S100A12', 'CCL2', 'CCL3', 'CCL4', 'THBS1', 'CD36', 'FUT4')
    gene.sets[['Conventional DC']] <- c('THBD', 'HLA-DRA', 'ITGAX', 'CD1C')
    gene.sets[['Plasmacytoid DC']] <- c('IL3RA', 'HLA-DRA')
    gene.sets[['B cells']] <- c('MS4A1', 'CD19', 'CD79A', 'CD79B', 'CD74', 'CD37', 'BTK', 'SYK', 'LYN', 'BLK')
    gene.sets[['NK cells']] <- c('NCAM1', 'CD7', 'XCL1', 'NCR1', 'NCR3', 'GZMB', 'FCGR3A', 'SLAMF7', 'KIR3DL1', 'KIR3DL2', 'KIR3DL2', 'KIR3DL3')
    gene.sets[['EM T cells']] <- c('CD3', 'CD4', 'CD8', 'CD44')
    gene.sets[['CM T cells']] <- c('CD3', 'CD4', 'CD8', 'CD44', 'SELL', 'CCR7')
    gene.sets[['naïve T cells']] <- c('CD3', 'CD4', 'CD8', 'CCR7', 'SELL')
    gene.sets[['activated T cells']] <- c('HLA-DRA', 'HLA-DRB', 'CD38')
    gene.sets[['Inhibitory molecules']] <- c('CTLA4', 'BTLA', 'IDO1', 'LAG3', 'PDCD1', 'HAVCR2', 'KIR2DL2', 'KIR2DL3', 'KIR3DL1', 'KIR3DL2', 'KIR3DX1')
    gene.sets[['Stimulatory molecules']] <- c('CD27', 'CD28', 'ICOS', 'CD40', 'CD80', 'CD86', 'TNFRSF9', 'TNFRSF4', 'TNFRSF18')
    gene.sets[['Tregs']] <- c('CD3', 'CD4', 'IL2R', 'IL7R', 'FOXP3', 'IL10', 'TGFB', 'GZMB', 'CCR4')

    par(mfrow=c(3,4))
    lapply(1:length(gene.sets), function(i) {
        genes <- gene.sets[[i]]
        g <- intersect(genes, rownames(comb.mat))
        if(length(g)==0) {
            return(NA)
        }
        else if(length(g)>1) {
            col <- scale(colSums(comb.mat[g, ]))[,1]
        } else {
            col <- scale(comb.mat[g, ])[,1]
        }
        col[is.nan(col)] <- 0
        range(col)
        col[col < -1] <- -1
        col[col > 1] <- 1
        plotEmbedding(ind.emb,
                      col=col,
                      mark.clusters = FALSE,
                      show.legend=FALSE,
                      main=names(gene.sets)[i], alpha=0.1)
    })

}

test.reference.10x <- function() {
    load("~/Projects/Noelia/data/b_cells.RData")
    vi <- colSums(cd>0)>1000; table(vi)
    cd <- cd[, vi]
    bcells.cd <- cd

    load("~/Projects/Noelia/data/cd4_t_helper.RData")
    vi <- colSums(cd>0)>1000; table(vi)
    cd <- cd[, vi]
    thelper.cd <- cd

    load("~/Projects/Noelia/data/cd14_monocytes.RData")
    vi <- colSums(cd>0)>1000; table(vi)
    cd <- cd[, vi]
    monocytes.cd <- cd

    load("~/Projects/Noelia/data/cd34.RData")
    vi <- colSums(cd>0)>1000; table(vi)
    cd <- cd[, vi]
    cd34.cd <- cd

    load("~/Projects/Noelia/data/cd56_nk.RData")
    vi <- colSums(cd>0)>1000; table(vi)
    cd <- cd[, vi]
    nk.cd <- cd

    load("~/Projects/Noelia/data/cytotoxic_t.RData")
    vi <- colSums(cd>0)>1000; table(vi)
    cd <- cd[, vi]
    cytotoxict.cd <- cd

    load("~/Projects/Noelia/data/memory_t.RData")
    vi <- colSums(cd>0)>1000; table(vi)
    cd <- cd[, vi]
    memoryt.cd <- cd

    load("~/Projects/Noelia/data/naive_cytotoxic.RData")
    vi <- colSums(cd>0)>1000; table(vi)
    cd <- cd[, vi]
    naivecytotoxic.cd <- cd

    load("~/Projects/Noelia/data/naive_t.RData")
    vi <- colSums(cd>0)>1000; table(vi)
    cd <- cd[, vi]
    naivet.cd <- cd

    load("~/Projects/Noelia/data/regulatory_t.RData")
    vi <- colSums(cd>0)>1000; table(vi)
    cd <- cd[, vi]
    regulatoryt.cd <- cd

    reference.cd <- cbind(
        bcells.cd,
        thelper.cd,
        monocytes.cd,
        nk.cd,
        cytotoxict.cd,
        memoryt.cd,
        naivecytotoxic.cd,
        naivet.cd,
        regulatoryt.cd
    )
    dim(reference.cd)
    ct <- c(
        rep("bcells", ncol(bcells.cd)),
        rep("thelper", ncol(thelper.cd)),
        rep("monocytes", ncol(monocytes.cd)),
        rep("nk", ncol(nk.cd)),
        rep("cytotoxict", ncol(cytotoxict.cd)),
        rep("memoryt", ncol(memoryt.cd)),
        rep("naivecytotoxic", ncol(naivecytotoxic.cd)),
        rep("naivet", ncol(naivet.cd)),
        rep("regulatoryt", ncol(regulatoryt.cd))
    )
    colnames(reference.cd) <- make.unique(paste0(ct, colnames(reference.cd)))
    names(ct) <- colnames(reference.cd)
    reference.cd <- Matrix(reference.cd, sparse=TRUE)
    head(reference.cd)

    rownames(reference.cd) <- make.unique(rownames(reference.cd))
    #vi <- rowSums(reference.cd)>50
    #reference.cd <- reference.cd[vi,]
    #dim(reference.cd)
    #genes.have <- intersect(rownames(reference.cd), rownames(comb.mat))
    genes.have <- rownames(comb.mat)
    length(genes.have)
    reference.cd <- reference.cd[genes.have,]

    ## run mudan
    reference.mat <- reference.cd
    reference.mat <- t(t(reference.mat)/colSums(reference.mat))
    reference.mat <- log10(reference.mat*1e6+1)
    reference.model <- model.lda(reference.mat, ct, nfeatures=1000)
    reference.emb <- tsne.lda(reference.mat, reference.model, ct, perplexity=50, show.legend=TRUE, mark.clusters=TRUE)

    plotEmbedding(reference.emb,
                  groups=ct,
                  mark.clusters = TRUE,
                  show.legend=TRUE)

    genes <- c('CD4', 'CD8A', 'ITGAM', 'CD72')
    par(mfrow=c(2,2))
    lapply(genes, function(g) {
        col <- scale(reference.mat[g, ])[,1]
        col[is.nan(col)] <- 0
        range(col)
        col[col < -1] <- -1
        col[col > 1] <- 1
        plotEmbedding(reference.emb,
                      col=col,
                      mark.clusters = FALSE,
                      show.legend=FALSE,
                      main=g, alpha=0.1)
    })

    ## apply reference.model to all datasets
    foo <- cbind(as.matrix(comb.mat[genes.have,]), as.matrix(reference.mat[genes.have,]))
    #test.reduction <- predict(reference.model, data.frame(t(foo)))$x
    models.all <- models
    models.all[['reference']] <- reference.model
    names(models.all)
    test.reduction <- do.call(cbind, lapply(models.all, function(model) {
                                         predict(model, data.frame(t(foo)))$x
                                     }))
    dim(test.reduction)
    #test.reduction <- predict(reference.model, data.frame(t(as.matrix(comb.mat))))$x

    d <- dist(test.reduction)
    perplexity = 150
    ncores = 10
    verbose=TRUE
    test.emb <- Rtsne.multicore::Rtsne.multicore(d, is_distance=TRUE, perplexity=perplexity, num_threads=ncores, verbose=verbose)$Y
    rownames(test.emb) <- labels(d)

    par(mfrow=c(1,2))
    plotEmbedding(test.emb, groups=sample, show.legend=TRUE, mark.clusters=TRUE, alpha=0.1)
    plotEmbedding(test.emb, groups=ct, show.legend=TRUE, mark.clusters=TRUE, alpha=0.1)

}

## Noelia's data
test.noelia <- function() {

    require(Matrix)
    path <- '~/Projects/Noelia/data/'
    files <- list.files(path)
    sn <- 'CRC0019|CRC004|CRC-001|CRC003'
    names <- files <- files[grepl(sn, files)]
                                        #names <- files <- files[!grepl('CLL', files)]
    names
    files <- paste0(path, files)
    names <- gsub('_cd.RData', '', names)
    names

    cds <- lapply(1:length(files), function(i) {
        f = files[i]
        n = names[i]
        load(f)

        vi <- colSums(cd)>1000
        print(table(vi))
        cd <- cd[, vi]

        ## summary <- log10(Matrix::colSums(cd))
        ## hist(summary, main=n)
        ## threshold0 <- quantile(summary, 0.15, na.rm=TRUE)
        ## threshold1 <- quantile(summary, 0.95, na.rm=TRUE)
        ## print(10^threshold0)
        ## print(10^threshold1)
        ## abline(v=threshold0, col='red')
        ## abline(v=threshold1, col='red')

        ## vi <- colSums(cd)>(10^threshold0)
        ## print(table(vi))
        ## cd <- cd[, vi]

        ## vi <- colSums(cd)<(10^threshold1)
        ## print(table(vi))
        ## cd <- cd[, vi]

        colnames(cd) <- paste0(n, '_', colnames(cd))
        as(cd, 'dgCMatrix')

    })
    names(cds) <- names
    length(cds)

    cds[['reference']] <- reference.cd
    names(cds)

    sample <- unlist(lapply(1:length(cds), function(i) {
        cd <- cds[[i]]
        n <- names(cds)[i]
        sample <- rep(n, ncol(cd))
        names(sample) <- colnames(cd)
        return(sample)
    }))
    length(sample)
    table(sample)

    k <- 100
    vi <- Reduce(intersect, lapply(cds, function(cd) {
                                vi <- rowSums(cd)>k
                                rownames(cd)[vi]
                            }))
    length(vi)
    head(vi)

    cds.filter <- lapply(cds, function(cd) cd[vi,])
    names(cds.filter)
    lapply(cds.filter, dim)

    reference.mat <- cds.filter[['reference']]
    reference.mat <- t(t(reference.mat)/colSums(reference.mat))
    reference.mat <- log10(reference.mat*1e6+1)
    reference.model <- model.lda(reference.mat, ct, nfeatures=1000)
    reference.emb <- tsne.lda(reference.mat, reference.model, ct, perplexity=50, show.legend=TRUE, mark.clusters=TRUE)

    plotEmbedding(reference.emb,
                  groups=ct,
                  mark.clusters = TRUE,
                  show.legend=TRUE)



    all.mat <- cbind(cds.filter[[1]], cds.filter[[2]], cds.filter[[3]], cds.filter[[4]], cds.filter[['reference']])
    #all.mat <- do.call(cbind, cds.filter)
    all.mat <- t(t(all.mat)/colSums(all.mat))
    all.mat <- log10(all.mat*1e6+1)
    test.reduction <- predict(reference.model, data.frame(t(as.matrix(all.mat))))$x
    d <- dist(test.reduction)
    perplexity = 150
    ncores = 10
    verbose=TRUE
    test.emb <- Rtsne.multicore::Rtsne.multicore(d, is_distance=TRUE, perplexity=perplexity, num_threads=ncores, verbose=verbose)$Y
    rownames(test.emb) <- labels(d)

    par(mfrow=c(1,2))
    plotEmbedding(test.emb, groups=sample, show.legend=TRUE, mark.clusters=TRUE, alpha=0.1)
    plotEmbedding(test.emb, groups=ct, show.legend=TRUE, mark.clusters=TRUE, alpha=0.1)

    par(mfrow=c(3,4))
    lapply(1:length(gene.sets), function(i) {
        genes <- gene.sets[[i]]
        g <- intersect(genes, rownames(all.mat))
        if(length(g)==0) {
            return(NA)
        }
        else if(length(g)>1) {
            col <- scale(colSums(all.mat[g, ]))[,1]
        } else {
            col <- scale(all.mat[g, ])[,1]
        }
        col[is.nan(col)] <- 0
        range(col)
        col[col < -1] <- -1
        col[col > 1] <- 1
        plotEmbedding(test.emb,
                      col=col,
                      mark.clusters = FALSE,
                      show.legend=FALSE,
                      main=names(gene.sets)[i], alpha=0.1)
    })

}



drawfigureConos2=function(p2,appname,jcl3.coarse=NULL,cell_ano_sampleType=NULL,cell_ano_sample=NULL,saveRDS=NULL){

  if (!is.null(jcl3.coarse)){
    pdf1=paste(appname,'.cells.conos.png',sep='') 
    a1=con$plotGraph(groups =jcl3.coarse,show.legend=TRUE,title=appname)
#    ggsave(pdf1,a1,width = 7,height=7)
  }
  
  
  if (!is.null(cell_ano_sampleType)){ 
    pdf1=paste(appname,'.sampleType.conos.png',sep='') 
    a2=con$plotGraph(groups =cell_ano_sampleType,show.legend=TRUE,title=appname)
#    ggsave(pdf1,a1,width = 7,height=7)
  }
  
  
  
  if (!is.null(cell_ano_sample)){
    pdf1=paste(appname,'.sample.conos.png',sep='') 
    a3=con$plotGraph(groups =cell_ano_sample,show.legend=TRUE,title=appname)
#    ggsave(pdf1,a1,width = 9,height=7)
  }
  
  
  if (!is.null(saveRDS)){
    f1=paste(appname,'_conos.rds',sep='')
    saveRDS(p2,f1)
  } 
  
    b=  cowplot::plot_grid(plotlist=list(a1,a3,a2), ncol=3, nrow=1)

	return(b)
  
}






badG=function(x){
  
  nms=rownames(x)
  ig_genes = c(grep("^IGJ", nms, v=T), 
               grep("^IGH",nms,v=T),
               grep("^IGK", nms, v=T), 
               grep("^IGL", nms, v=T))
  
  bad_genes = unique(c(grep("^RP[LKS]", nms, v=T),c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes)))
  index=nms %in% c(bad_genes,ig_genes)
  x=x[!index,]
  print(table(index))
  return(x)
}
