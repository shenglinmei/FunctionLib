



#  output largeVis embdding 
embedGraph2=function (graph,method = "largeVis",ndim=20, M = 1, gamma = 1, alpha = 0.1, 
          perplexity = NA, sgd_batches = 1e+08, seed = 1, verbose = TRUE) 
{
  if (method != "largeVis") {
    stop("currently, only largeVis embeddings are supported")
  }
  wij <- as_adj(graph, attr = "weight")
  if (!is.na(perplexity)) {
    wij <- conos:::buildWijMatrix(wij, perplexity = perplexity, 
                                  threads = n.cores)
  }
  coords <- conos:::projectKNNs(wij = wij, dim = ndim, verbose = verbose, 
                                sgd_batches = sgd_batches, gamma = gamma, M = M, seed = seed, 
                                alpha = alpha, rho = 1, threads = 8)
  colnames(coords) <- V(graph)$name
  embedding <<- coords
  return(embedding)
}


#con$clusters$multi$groups
#con$clusters$walktrap$groups
#con$clusters$infomap$groups
#con$clusters$multi$groups


conCmunity=function(con){
 # t=embedGraph2(con$graph)
 # emb <-Rtsne.multicore::Rtsne.multicore(t(t),  perplexity=50, num_threads=10)$Y
 # rownames(emb) <- colnames(t)
  con$findCommunities(method=infomap.community)
  con$findCommunities(method=walktrap.community)
  con$findCommunities(method=label.propagation.community)
 # con$tSNE=emb
  return(con)
  }


conEmb=function(con){
  t=embedGraph2(con$graph)
  emb <-Rtsne.multicore::Rtsne.multicore(t(t),  perplexity=50, num_threads=10)$Y
  rownames(emb) <- colnames(t)
  
  pdf('emn.pdf')
  plotEmbedding(emb,groups=con$clusters$multi$groups, mark.clusters=TRUE, alpha=0.1)
  dev.off()
  return(emb)
}
  


runConos_Par=function(con,k,dtype,appname,ncomps=50,n.odgenes=1000,ano1,ano2){
#  k=15
#  dtype='CPCA'
#  ncomps=20
#  n.odgenes=200
#  ano1=anoCell
#  ano2=anoCell2
  con$buildGraph(k=k,space = dtype,ncomps = ncomps, n.odgenes = n.odgenes)
  con$findCommunities()
  
  
  appname2=paste(c(appname,dtype,k,ncomps,n.odgenes),collapse='_')
  f1=paste(appname2,'.conos.png',sep='')
  a=con$plotGraph()
  ggsave(f1,a,width = 7,height=7)
  
  draw_conos_sgd_par(con,appname2,1e+09,109,ano1,ano2)
  draw_conos_sgd_par(con,appname2,4e+09,409,ano1,ano2)
  draw_conos_sgd_par(con,appname2,5e+08,508,ano1,ano2)

  r=list('con'=con,'appname'=appname2)
  return(r)

}



draw_conos_sgd_par=function(con,appname,sgd_batches,id,cell_ano1,cell_ano2){
  con$embedGraph(sgd_batches = sgd_batches)
  a1=con$plotGraph()
  a2= con$plotGraph(groups=cell_ano1)
  a3= con$plotGraph(groups=cell_ano2)
 
  b=  cowplot::plot_grid(plotlist=list(a1,a2,a3), ncol=3, nrow=1)

  
  f1=paste(appname,'_',id,'_conos_clustering.png',sep='')
  ggsave(f1,b,width = 16,height=6)

}


#anoCell=con$clusters$multi$groups
#runConos_Par_mutiplePar(con,appname,anoCell,anoCell)

runConos_Par_mutiplePar=function(con,appname,anoCell,anoCellt,ncore=10){
  con$n.cores=ncore
  
  #Lcon=runConos_Par(con,30,'CPCA',appname,10,200,anoCell,anoCell)

  Lcon=runConos_Par(con,30,'CPCA',appname,50,1000,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
  
  
  Lcon=runConos_Par(con,30,'CPCA',appname,30,500,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
  
  
  Lcon=runConos_Par(con,30,'PCA',appname,50,1000,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
  
  Lcon=runConos_Par(con,30,'PCA',appname,30,500,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
  
  Lcon=runConos_Par(con,30,'genes',appname,50,1000,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
  
  Lcon=runConos_Par(con,30,'genes',appname,30,500,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
  
  
  Lcon=runConos_Par(con,15,'CPCA',appname,50,1000,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
  
  
  Lcon=runConos_Par(con,15,'CPCA',appname,30,500,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
  
  
  Lcon=runConos_Par(con,15,'PCA',appname,50,1000,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
  
  Lcon=runConos_Par(con,15,'PCA',appname,30,500,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
  
  Lcon=runConos_Par(con,15,'genes',appname,50,1000,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
  
  Lcon=runConos_Par(con,15,'genes',appname,30,500,anoCell,anoCell)
  drawfigureConos(Lcon[['con']],Lcon[['appname']],jcl3.coarse=anoCell,cell_ano_sampleType=anoCellt,cell_ano_sample=anoCell,saveRDS=TRUE)
}




conos_Gene_exp<- function(con,gs,appname,alpha=0.3,size=0.4){
  
  explis=lapply(con$samples,function(x) t(x$counts))
  exp=mergeDat2(explis)
  gs2=intersect(gs,rownames(exp))
  lrow=ceiling(length(gs2)/3)
  
  lis=list()
  for (gene in gs2){
    print(gene)
    t=exp[gene,]
    a=con$plotGraph(colors =t,title=gene,alpha =alpha,size=size)
    lis[[gene]]=a
  }
  b=  cowplot::plot_grid(plotlist=lis, ncol=3, nrow=lrow)
  fout=paste(appname,'.expConos.png',sep='')
  ggsave(fout,b,width = 6,height=2*lrow)


}




convert.factors.to.strings.in.dataframe <- function(dataframe)
{
  class.data  <- sapply(dataframe, class)
  factor.vars <- class.data[class.data == "factor"]
  for (colname in names(factor.vars))
  {
    dataframe[,colname] <- as.character(dataframe[,colname])
  }
  return (dataframe)
}



mappingStat10X <-function(path,appname){
  key=c('Estimated.Number.of.Cells','Mean.Reads.per.Cell','Median.Genes.per.Cell','Number.of.Reads','Total.Genes.Detected','Median.UMI.Counts.per.Cell')
  res=NULL
  
  path=apply(data.frame(path),1,function(x) strsplit(x,'outs')[[1]][1])
  
  nn=apply(data.frame(path),1,function(x) {
    tmp=strsplit(x,'/')[[1]]
    tmp[length(tmp)]}
           )
  path=apply(data.frame(path),1,function(x) paste(x,'/outs/metrics_summary.csv',sep=''))
  for( f in path){
    
    dat=read.csv(f,sep=',',header=T)[,key]
    res=rbind(res,convert.factors.to.strings.in.dataframe(dat))

  }
  rownames(res)=nn
}










