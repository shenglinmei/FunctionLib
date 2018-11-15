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





runP2Genes=function(p2,gs,nPcs=12,k=15,n.cores=10 ,alpha = 0.3 ){
  p2$adjustVariance(plot = F, gam.k = 10,alpha =alpha)
  g2=intersect(p2$misc$odgenes,gs)
  print(length(g2))
  p2$calculatePcaReduction(nPcs = 12, odgenes = g2, 
                           maxit = 1000)
  p2$makeKnnGraph(k = k, type = "PCA", center = TRUE, weight.type = "none", 
                  n.cores = n.cores, distance = "cosine")                             
  p2$getKnnClusters(method = igraph::multilevel.community, 
                    type = "PCA", name = "multilevel")
  p2$getEmbedding(type = 'PCA', embeddingType = 'tSNE', perplexity = 50)
  return(p2)
}


draw_conos_sgd=function(con,appname,sgd_batches,id,cell_ano=NULL){
  con$embedGraph(sgd_batches = sgd_batches)
  f1=paste(appname,'_',id,'_conos_clustering.png',sep='')
  a=con$plotGraph(groups=cell_ano_sampleType)
  ggsave(f1,a,width = 7,height=7)
  
  if(is.null(cell_ano)){
    f1=paste(appname,'_',id,'_conos_clustering2.png',sep='')
    a2= con$plotGraph(groups=cell_ano)
    ggsave(f1,a2,width = 7,height=7)
    f1=paste(appname,'_',id,'_conos_clustering2.rds',sep='')
    saveRDS(con,f1)
  } 
  
}

#draw_fig(con,appname2,1e+09,109)
#draw_fig(con,appname2,3e+08,308)
#draw_fig(con,appname2,6e+08,608)
#draw_fig(con,appname2,2e+09,209)
#draw_fig(con,appname2,4e+09,409)
#draw_fig(con,appname2,9e+09,909)








# drawfigureP2(p2,appname,jcl3.coarse=NULL,cell_ano_sampleType=NULL,cell_ano_sample=NULL,saveRDS=NULL)

drawfigureP2=function(p2,appname,jcl3.coarse=NULL,cell_ano_sampleType=NULL,cell_ano_sample=NULL,saveRDS=NULL,alpha=0.1){
  pdf1=paste(appname,'.clustering.tsn.png',sep='') 
  png(pdf1,height=600,width=600)
  p2$plotEmbedding(groups=p2$clusters$PCA$multilevel,type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='clusters (tSNE)')
  dev.off()
  
  
  if (!is.null(jcl3.coarse)){
    pdf1=paste(appname,'.cells.tsn.png',sep='') 
    png(pdf1,height=600,width=600)
    p2$plotEmbedding(type='PCA',groups =jcl3.coarse,embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='clusters (tSNE)')
    dev.off()
  }
  
  
  if (!is.null(cell_ano_sampleType)){ 
    pdf1=paste(appname,'.cells.tsn.sampleType.png',sep='') 
    png(pdf1,height=600,width=600)
    p2$plotEmbedding(type='PCA',groups =cell_ano_sampleType,embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='clusters (tSNE)')
    dev.off()
  }
  
  
  
  if (!is.null(cell_ano_sample)){
    pdf1=paste(appname,'.cells.tsn.sample.png',sep='') 
    png(pdf1,height=600,width=600)
    p2$plotEmbedding(type='PCA',groups =cell_ano_sample,embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=alpha,main='clusters (tSNE)')
    dev.off()
  }
  
  
  if (!is.null(saveRDS)){
    f1=paste(appname,'_p2combined.rds',sep='')
    saveRDS(p2,f1)
  } 
  
}


#drawfigureP2(p2,appname,jcl3.coarse=cell_ano_cell,cell_ano_sampleType=cell_ano_sampleType,cell_ano_sample=cell_ano_sample,saveRDS=NULL)
  



geneExpP2=function(p2,gs,lcol,lrow,appanme){
#  gs3=c('CD8A','CD8B','CD4','FCGR3A','FGFBP2','CD3E','CD3D','SH2D1A','TNFRSF4','CTLA4','TIGIT','GZMA','GZMB','GZMM','GZMK','GZMH','LAG3')
  filename=paste(appname,'exp.png',sep='')
  mar = c(0.5, 0.5, 0.5, 0.5)
  mgp = c(2, 0.65, 0)
  panel.size=300
  cex=1.1
  mark.cluster.cex = 0.8
  mark.clusters = TRUE
  
  png(file = filename, height = lrow * panel.size, 
      width = lcol * panel.size)
  par(mfrow = c(lrow, lcol), mar = mar, mgp = mgp, 
      cex = cex)
  
  for (gene in gs){
    p2$plotEmbedding(type='PCA',embeddingType='tSNE',colors=p2$counts[,gene],shuffle.colors=F,mark.cluster.cex=1,alpha=0.2,main=gene)
  }
  dev.off()
  
}

#gs=c('CD38','IGLL1','SOX4','CD9','CD52')
#geneExpP2(p2,gs,3,2,appanme)


geneExpConos=function(con,p2,gs,lrow,lcol,appanme,alpha =0.3){
  lis=list()
  for (gene in gs){
    print(gene)
    t=p2$counts[,gene]
    a=con$plotGraph(colors =t,title=gene,alpha =0.3)
    lis[[gene]]=a
  }
  b=  cowplot::plot_grid(plotlist=lis, ncol=lrow, nrow=lcol)
  fout=paste(appname,'expConos.png',sep='')
  ggsave(fout,b,width = 2.5*lcol,height=1.2*lrow)
}

#gs=c('CD38','CD8A','CD8B')
#geneExpConos(con,p2,gs,3,2,appanme)





#drawfigureConos=function(p2,appname,jcl3.coarse=NULL,cell_ano_sampleType=NULL,cell_ano_sample=NULL,saveRDS=NULL)
drawfigureConos=function(p2,appname,jcl3.coarse=NULL,cell_ano_sampleType=NULL,cell_ano_sample=NULL,saveRDS=NULL){
  pdf1=paste(appname,'.clustering.tsn.conos.png',sep='') 
  a1=con$plotGraph()
  ggsave(pdf1,a1,width = 7,height=7)
  
  
  if (!is.null(jcl3.coarse)){
    pdf1=paste(appname,'.cells.conos.png',sep='') 
    a1=con$plotGraph(groups =jcl3.coarse,show.legend=TRUE,title=appname)
    ggsave(pdf1,a1,width = 7,height=7)
  }
  
  
  if (!is.null(cell_ano_sampleType)){ 
    pdf1=paste(appname,'.sampleType.conos.png',sep='') 
    a1=con$plotGraph(groups =cell_ano_sampleType,show.legend=TRUE,title=appname)
    ggsave(pdf1,a1,width = 7,height=7)
  }
  
  
  
  if (!is.null(cell_ano_sample)){
    pdf1=paste(appname,'.sample.conos.png',sep='') 
    a1=con$plotGraph(groups =cell_ano_sample,show.legend=TRUE,title=appname)
    ggsave(pdf1,a1,width = 9,height=7)
  }
  
  
  if (!is.null(saveRDS)){
    f1=paste(appname,'_conos.rds',sep='')
    saveRDS(p2,f1)
  } 
  
}

#drawfigureConos(con,appname,jcl3.coarse=cell_ano_cell,cell_ano_sampleType=cell_ano_sampleType,cell_ano_sample=cell_ano_sample,saveRDS=NULL)






# DEcaculate(p2,appname,conosCluster,removeGene=TRUE)
# caculate differnential expressed gene based on cluster group
DEcaculate=function(p2,appname,conosCluster,removeGene=NULL){
  
  de1 <- p2$getDifferentialGenes(groups=conosCluster)
  
  
  f1=paste(appname,'_diffGene.rds',sep='')
  saveRDS(de1,f1)
  for(n in names(de1)){
    x=de1[[n]]
    z <- x[order(-x$Z),]
    if(!is.null(removeGene)){
      index=grepl('^RP[LKS]',rownames(z))
      z=z[!index,]
    }
    markers=rownames(z)[1:100]
    x <- as.matrix(p2$counts[names(conosCluster),markers])
    ## trim outliers
    x <- apply(x, 2, function(xp) {
      qs <- quantile(xp,c(0.01,0.98))
      xp[xp<qs[1]] <- qs[1]
      xp[xp>qs[2]] <- qs[2]
      xp
    })
    x <- x[,(apply(x,2,sd) != 0)]
    x <- t(scale(x))
    ## sample 2000 cells for plotting
    #x <- x[,sample(ncol(x),2000)]
    o <- order(as.numeric(as.character(conosCluster[colnames(x)])))
    annot <- data.frame(cluster=conosCluster[colnames(x)],row.names = colnames(x))
    
    annot$dtype='other'
    annot[as.character(annot[,1])==n,'dtype']=n
    annot$dtype=as.factor(annot$dtype)
    
    
    pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
    ## draw heatmap
    
    fout=paste(appname,'_',n,'_marker.heatmap.new.png',sep='')
    rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
    pheatmap(x[,o],labels_col=FALSE,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,
             cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=5,height=6,   #3.5*0.02*length(markers),
             breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
    
    
    iid= paste(appname,'_cluster_',n,sep='')
    GOanalysis(markers,iid) 
  }

}

#conosCluster=con$clusters$multi$groups
#DEcaculate(p2,appname,conosCluster,removeGene=TRUE)






# MakeWebConos(p2,appanme,conosTSN,conosCluster,jcl3.coarse2)
#conosTSN=con$embedding
#conosCluster=con$clusters$multi$groups

MakeWebConos<-function(p2,appname,conosTSN,conosCluster,jcl3.coarse2,cell_ano_sample=NULL,cell_ano_sampleType=NULL){
  n=appname
  p2$embeddings$PCA$conosEb=t(conosTSN)
  jfac2=as.factor(conosCluster)
  
  S_cells_ano=as.factor(jcl3.coarse2)  
  
  cells.in.app <- rownames(p2$counts,jfac2)
  ## make hierarchical aspects
  hdea <- p2$getHierarchicalDiffExpressionAspects(type='PCA',clusterName='multilevel',z.threshold=3)
  ## make metadata   
  metadata.forweb <- list();
  metadata.forweb$multilevel <- p2.metadata.from.factor(p2$clusters$PCA$multilevel,displayname='Multilevel')
  ## Set colors manually
  ##p1 <- colorRamps::primary.colors(n = nlevels(jfac2))
  p1 <- rainbow(length(levels(jfac2)), s=1, v=1)
  names(p1) <- levels(jfac2)
  metadata.forweb$conos <- p2.metadata.from.factor(jfac2[cells.in.app], pal=p1,displayname='conos')
  ## get de sets
  
  
  p1 <- rainbow(length(levels(S_cells_ano)), s=1, v=1)
  names(p1) <- levels(S_cells_ano)
  metadata.forweb$jointsamp <- p2.metadata.from.factor(S_cells_ano[cells.in.app], pal=p1,displayname='jointsamp')

  if(!is.null(cell_ano_sample)){
    cell_ano_sample=as.factor(cell_ano_sample)
    p1 <- rainbow(length(levels(cell_ano_sample)), s=1, v=1)
    names(p1) <- levels(cell_ano_sample)
    metadata.forweb$sample <- p2.metadata.from.factor(cell_ano_sample[cells.in.app], pal=p1,displayname='sample')
  } 
  
  if(!is.null(cell_ano_sample)){
    cell_ano_sampleType=as.factor(cell_ano_sampleType)
    p1 <- rainbow(length(levels(cell_ano_sampleType)), s=1, v=1)
    names(p1) <- levels(cell_ano_sampleType)
    metadata.forweb$sampleType <- p2.metadata.from.factor(cell_ano_sampleType[cells.in.app], pal=p1,displayname='sampleType')
  } 
  
  deSets <- get.de.geneset(p2, groups=jfac2, prefix='conosDE_')
  ## Collect the genesets
  genesets = c(deSets, hierDiffToGenesets(hdea))
  appmetadata = list(apptitle=n)
  p2$makeGeneKnnGraph();
  ## make app
  wp <- make.p2.app(p2, additionalMetadata = metadata.forweb, geneSets = genesets,
                    dendrogramCellGroups=p2$clusters$PCA[[1]],show.clusters=F,
                    appmetadata=appmetadata)
  wp$serializeToStaticFast(paste0(n,'.bin'))
  
}  

#conosTSN=con$embedding
#conosCluster=con$clusters$multi$groups
#MakeWebConos(p2,appname,conosTSN,conosCluster,cell_ano_cell,cell_ano_sample=cell_ano_sample,cell_ano_sampleType=cell_ano_sampleType)
  



creatLableList2=function(jcl3.coarse){
  #  key='allName'
  cell_ano_cell=as.character(jcl3.coarse)
  names(cell_ano_cell)=names(jcl3.coarse)
  cell_ano_sample=names(cell_ano_cell)
  names(cell_ano_sample)=names(cell_ano_cell)
  cell_ano_sample=apply(data.frame(cell_ano_sample),1,function(x) strsplit(x[1],'_')[[1]][1])
  cell_ano_sampleType=apply(data.frame(cell_ano_sample),1,function(x) strsplit(x[1],'-')[[1]][2])
  cellano=list('cell'=cell_ano_cell,'sample'=cell_ano_sample,'sampltType'=cell_ano_sampleType)
  print(table(cell_ano_cell))
  print(table(cell_ano_sampleType))
  print(table(cell_ano_sample))
  
  return(cellano)
}


mergeDat2=function(raw.mats){
  
  genelists <- lapply(raw.mats, function(x) rownames(x))
  str(genelists)
  commongenes <- Reduce(intersect,genelists)
  
  
  matrices2 <- mapply(function(m, name) {
 #   colnames(m) <- paste(name, colnames(m), sep='_');
    m[commongenes,]
  },
  raw.mats,
  names(raw.mats))
  cellNum=unlist(lapply(matrices2, function(x) ncol(x)))
  print('#commongenes')
  print(length(commongenes))
  
  print('#cell numbers')
  print(cellNum)
  bigM2 <- Reduce(cbind, matrices2)
  
  return(bigM2)
}




