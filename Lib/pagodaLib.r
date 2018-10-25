library(reshape2)
library(xtable)
library(pagoda2)
library(nbHelpers)
library(reshape2)
library(Matrix)
library(clusterMatch)
library(igraph)
library(ggplot2)
#library(nbHelpers)                  
library(parallel)
#library(fastSave)
library(pagoda2)
#library(nbHelpers)
library(reshape2)
library(pheatmap)
library(conos)
#library('betareg')
library(DESeq2)
library(Matrix)
library(conos)
library(GOstats)
library(RSQLite)
library(biomaRt)
suppressPackageStartupMessages(library(org.Hs.eg.db))



#  run pagoda apps 
runPagoda=function(cd,appname='none', n.cores = 1, batch = NULL, n.odgenes = 3000, nPcs = 100, 
                          k = 30, perplexity = 50, log.scale = TRUE, trim = 10, keep.genes = NULL, 
                          min.cells.per.gene = 30, get.largevis = TRUE, get.tsne = TRUE, 
                          make.geneknn = TRUE) {

    rownames(cd) <- make.unique(rownames(cd))
    p2 <- Pagoda2$new(cd, n.cores = n.cores, batch = batch, keep.genes = keep.genes, 
                      trim = trim, log.scale = log.scale, min.cells.per.gene = min.cells.per.gene)
    
    
    pv=paste(appname,'.adjustVariance.pdf',sep='')
    pdf(pv)
    p2$adjustVariance(plot = T, gam.k = 10)
    dev.off()
    
    p2$calculatePcaReduction(nPcs = nPcs, n.odgenes = n.odgenes, 
                             maxit = 1000)
    p2$makeKnnGraph(k = k, type = "PCA", center = TRUE, weight.type = "none", 
                    n.cores = n.cores, distance = "cosine")
    p2$getKnnClusters(method = igraph::infomap.community, type = "PCA", 
                      name = "infomap")
    p2$getKnnClusters(method = igraph::multilevel.community, 
                      type = "PCA", name = "multilevel")
    
    p2$getKnnClusters(method = igraph::walktrap.community, type = 'PCA', name='walktrap')
    
    p2$getDifferentialGenes(type='PCA',verbose=T,clusterType='multilevel')
    
    
    p2$getEmbedding(type = 'PCA', embeddingType = 'tSNE', perplexity = 50)
#    M <- 30
#    p2$getEmbedding(type = 'PCA', embeddingType = 'largeVis', M = M, perplexity = perplexity, gamma = 1 / M, alpha = 1 )
    
    
    pdf1=paste(appname,'.tsn.pdf',sep='') 
    pdf(pdf1)
    p2$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='clusters (tSNE)')
    dev.off()
    
    
    
    return(p2)
}



##    merge data from  list oblects
mergeDat=function(raw.mats){
  
  genelists <- lapply(raw.mats, function(x) rownames(x))
  str(genelists)
  commongenes <- Reduce(intersect,genelists)
  
  
  matrices2 <- mapply(function(m, name) {
    colnames(m) <- paste(name, colnames(m), sep='_');
    m[commongenes,]
  },
  raw.mats,
  names(raw.mats))
  cellNum=unlist(lapply(matrices2, function(x) ncol(x)))
  print('#commongenes')
  print(length(commongenes))
  
  print('#cell numbers')
  print(cellNum)
  
  return(matrices2)
}

#matrices=mergeDat(raw.mats)

creatLableList=function(matrices,res,key){
#  key='allName'
  celllists <- lapply(matrices, function(x) colnames(x))
  celss=unlist(celllists)
  names(celss)=celss
  
  for ( i in seq(length(matrices))){
    sampID=as.character(res[i,1])
    value=as.character(res[i,key])
    celss[celllists[[sampID]]]=value
  }
  print(table(celss))
  return(celss)
}

#ceLable=creatLableList(matrices,res,'allName')


makeWebLable=function(appName,p2,jfac2,n.cores=5){

  jfac2=as.factor(jfac2)
  cells.in.app <- rownames(p2$counts,jfac2)
  ## make hierarchical aspects
  hdea <- p2$getHierarchicalDiffExpressionAspects(type='PCA',clusterName='multilevel',z.threshold=3, n.cores = 5)
  ## make metadata   
  metadata.forweb <- list();
  metadata.forweb$multilevel <- p2.metadata.from.factor(p2$clusters$PCA$multilevel,displayname='Multilevel')
  ## Set colors manually
  #  myinfo=p2$clusters$PCA$infomap
  #  p1 <- rainbow(length(levels(myinfo)), s=1, v=1)
  #  names(p1) <- levels(myinfo)
  #  metadata.forweb$infomap <- p2.metadata.from.factor(myinfo[cells.in.app], pal=p1,displayname='Infomap')
  metadata.forweb$infomap <- p2.metadata.from.factor(p2$clusters$PCA$infomap, displayname = 'Infomap', start=0, end=0.5, s = 1, v=0.7)
  
  
  p1 <- rainbow(length(levels(jfac2)), s=1, v=1)
  names(p1) <- levels(jfac2)
  metadata.forweb$Label <- p2.metadata.from.factor(jfac2[cells.in.app], pal=p1,displayname='Label')
  ## get de sets
  
  #    p1 <- rainbow(length(levels(S_cells_ano)), s=1, v=1)
  #    names(p1) <- levels(S_cells_ano)
  
  #    metadata.forweb$jointsamp <- p2.metadata.from.factor(S_cells_ano[cells.in.app], pal=p1,displayname='jointsamp')
  
  deSets <- get.de.geneset(p2, groups=p2$clusters$PCA$multilevel, prefix='de_')
  ## Collect the genesets
  genesets = c(deSets, hierDiffToGenesets(hdea))
  
  #  genesets = hierDiffToGenesets(hdea)
  appmetadata = list(apptitle=appName)
  p2$makeGeneKnnGraph();
  ## make app
  wp <- make.p2.app(p2, additionalMetadata = metadata.forweb, geneSets = genesets,
                    dendrogramCellGroups=p2$clusters$PCA[[1]],show.clusters=F,
                    appmetadata=appmetadata)
  wp$serializeToStaticFast(paste0(appName,'.bin'))
  
  fout=paste(appName,'.rds',sep='')
  saveRDS(p2,fout)
}




makeWeb=function(p2,appName,n.cores=5){
  
  hdea <- p2$getHierarchicalDiffExpressionAspects(type = "PCA", 
                                                  clusterName = "multilevel", z.threshold = 3, n.cores = 5)
  metadata.forweb <- list()
  metadata.forweb$multilevel <- p2.metadata.from.factor(p2$clusters$PCA$multilevel, 
                                                        displayname = "Multilevel")
  
  metadata.forweb$infomap <- p2.metadata.from.factor(p2$clusters$PCA$infomap, displayname = 'Infomap', start=0, end=0.5, s = 1, v=0.7)
  #    metadata.forweb$walktrap <- p2.metadata.from.factor(p2$clusters$PCA$walktrap, displayname = 'Walktrap', s = 0.5)
  extraWebMetadata = NULL
  metadata.forweb <- c(metadata.forweb, extraWebMetadata)
  
  deSets <- get.de.geneset(p2, groups = p2$clusters$PCA[[1]], prefix = 'de_')
  
  
  genesets <- hierDiffToGenesets(hdea)
  
  genesets <- deSets
  
  
  appmetadata = list(apptitle = appName)
  cat("Making KNN graph...\\n")
  p2$makeGeneKnnGraph(n.cores = n.cores)
  myPagoda2WebObject=make.p2.app(p2, additionalMetadata = metadata.forweb,
                                 geneSets = genesets, 
                                 dendrogramCellGroups = p2$clusters$PCA$multilevel,
                                 show.clusters = F, 
                                 appmetadata = appmetadata)
  
  
  # Save serialised web object, RDS app and session image
  #myPagoda2WebObject$serialiseToStatic(text.file.directory = './tmp', binary.filename = paste0(appName,'.bin'))
  myPagoda2WebObject$serializeToStaticFast(binary.filename = paste0(appName,'.bin'))
  
  #    invisible(p2)
  
  
}



GOanalysis=function(markers,n){

  ENTREZID=unlist(mget(markers, org.Hs.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID=ENTREZID[!is.na(ENTREZID)]
  
  
  for(function_type in c("BP", "CC", "MF")){
    
    param <- new("GOHyperGParams", geneIds=ENTREZID,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology=function_type,pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp <- hyperGTest(param)
    sumTable <- summary(hyp)
    
    
    
    david=sumTable[1:20,]
    david$Pvalue=-log(david[,2])
    termNumber=nrow(david)
    
    
    library(ggplot2)
    p1 <- ggplot(data=david, aes(x=Pvalue, y=Term, size=Count, colour = factor(david$Count)))
    p1 <- p1 + geom_point()
    p1 <- p1 + guides(color = FALSE)
    p1 <- p1 + theme(panel.grid.major = element_line(colour='blue'),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank())
    p1 <- p1 + xlab(paste("-log10(Pvalue)", sep="")) + ylab("")
    p1 <- p1 + labs(title=paste("DAVID:", function_type, sep=""))
    p1 <- p1 + theme(axis.text.x=element_text(size=10, face="plain", colour ='black'))
    p1 <- p1 + theme(axis.text.y=element_text(size=6, face="plain", colour ='black'))
    #p1=p1+theme(axis.title.y=element_text(size=10,face="plain",colour ='black'))
    #p1=p1+theme(legend.position="bottom")
    p1 <- p1 + xlim(min(david$Pvalue), max(david$Pvalue))
    #p1=p1+ylim(-10,+15)
    print(p1)
    ggsave(file=paste(n,'_' ,function_type, ".png", sep=""), scale=0.8, dpi=600, width = 7, height=1+0.25*termNumber) 
  } 
  
}







#  draw figure for single genes
#plotWithGroupsGene(p2ens$p2objs,'PLAUR',jcl3.coarse, file='Myoloid_emb.jcl3.PLAUR.png',verbose=T,panel.size=300)
plotWithGroupsGene=function(p2objs,gene,groups = NULL, filename = NULL, panel.size = 300, mark.cluster.cex = 0.8, 
                            mar = c(0.5, 0.5, 0.5, 0.5), mgp = c(2, 0.65, 0), cex = 0.85, 
                            type = "PCA", embeddingType = "tSNE", mark.clusters = TRUE, 
                            verbose = FALSE) 
{
  require(Cairo)
  panel.dims <- getParMfrow(length(p2objs))
  if (verbose) 
    cat("Panel dimensions are ", panel.dims[1], " by ", panel.dims[2], 
        "\\n")
  if (!is.null(filename)) 
    CairoPNG(file = filename, height = panel.dims[1] * panel.size, 
             width = panel.dims[2] * panel.size)
  par(mfrow = c(panel.dims[1], panel.dims[2]), mar = mar, mgp = mgp, 
      cex = cex)
  lapply(names(p2objs), function(dn) {
    d <- p2objs[[dn]]
    g1 <- as.factor(groups)
    colors <- NULL
    if (!any(names(g1) %in% rownames(d$counts))) {
      g1 <- NULL
      cell.names <- rownames(d$counts)
      colors <- rep("grey70", length(cell.names))
      names(colors) <- cell.names
    }
    d$plotEmbedding(type = type, embeddingType = embeddingType, 
                    alpha = 0.2, min.group.size = 0, mark.clusters = mark.clusters, 
                    mark.cluster.cex = mark.cluster.cex, do.par = F, 
                    colors = d$counts[,gene])
    
    
    
    legend(x = "topleft", bty = "n", legend = dn)
  })
  if (!is.null(filename)) 
    dev.off()
  invisible(NULL)
}




