#  tracking embdding 


emd=readRDS('/home/meisl/Workplace/BMME/a.data/test_cell/cell/embdding.rds')

emd=con$embedding

cname2=colnames(emd)

n.anoCell=allano[['cell']][cname2]
n.anoSample=allano[['sample']][cname2]
n.anoSampleType=allano[['sampltType']][cname2]





lise=list()
	for (dt in c('Tumor','Whole','Involved','Noninvolved')){

index=n.anoSampleType==dt
tembdding=t(emd[,index])
t.anoSampleType=n.anoSampleType[index]
t.cell=n.anoCell[index]

index=t.cell %in% c('progenitors','prog_B_lymphoid','erythroid','PDC','neutrophils')
tembdding2=tembdding[index,]
t.cell2=t.cell[index]
pdf('t.pdf')
conos::embeddingPlot(tembdding2, groups=t.cell2,title=dt,alpha=0.3)
dev.off()


a6=conos::embeddingPlot(tembdding, groups=t.cell,title=dt,alpha=0.1)

lise[[dt]]=a6


}


b=  cowplot::plot_grid(plotlist=lise, ncol=4, nrow=1)

pdf1='t.png'
ggsave(pdf1,b,width = 16,height=4)






# GSEA




library('liger')
library(qusage)


msigdb.c2 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c7.all.v6.2.symbols.gmt')
msigdb.c5 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c5.all.v6.2.symbols.gmt')
msigdb.c6 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c6.all.v6.2.symbols.gmt')
msigdb.c7 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c7.all.v6.2.symbols.gmt')

sign.sets <- list(c2=msigdb.c2,c5=msigdb.c5,c6=msigdb.c6,c7=msigdb.c7)


fc <- fc[!is.na(fc)]
r=lapply(namedNames(sign.sets), function(n) {
  cat(paste0('\tProcessing set ',n,'\n'))
  signset <- sign.sets[[n]]
  gsea.res <- bulk.gsea(fc, signset, mc.cores=12)
  gsea.res[order(gsea.res$q.val),]
})


fo=paste(dtype,'.gsea.rds',sep='')
saveRDS(r,fo)


# Differential expressed gene add Fc
table(anoSampleType)

de1 <- p2$getDifferentialGenes(groups=anoSampleType,z.threshold =0)
de1=lapply(de1,function(x) x[order(x$Z,decreasing=T),])
ncell=table(anoSampleType)

res=list()
    if (cell %in% names(ncell)){
        tmp=de1[[cell]]
        tmp[is.na(tmp[,1]),1]=0
        tmp$rank=rank(1/(10000+tmp[,1]))
        
        n1=names(anoSampleType[anoSampleType==cell])
        n2=names(anoSampleType[anoSampleType!=cell])
        fc=rowMeans(t(p2$counts[n1,]))-rowMeans(t(p2$counts[n2,]))
        tmp$fc=fc[rownames(tmp)]
        res[[cell]]=tmp
  
      }
    
  
  f1=paste(appname,'_diffGene.all.rds',sep='')
  saveRDS(res,f1)  


# node to tSNE


source('/home/meisl/bin/FunctionLib/Lib/pagodaLib.r')

saveGraph <- function(g, path, directed.format=F) {
  file.conn <- file(path)
  edge.list <- igraph::get.edgelist(g) %>% as.factor() %>% as.integer() %>% 
    matrix(ncol=2) %>% cbind(igraph::E(g)$weight)

  if (directed.format) {
    edge.list <- rbind(edge.list, edge.list[,c(2, 1, 3)])
  }
  
  edge.list <- apply(edge.list, 1, paste, collapse=" ")
  
  writeLines(edge.list, file.conn)
  
  close(file.conn)
}

length(V(p2$graphs$PCA))

saveGraph(p2$graphs$PCA, "10x_9k_graph.lst")
saveGraph(con$graph, "con.lst")




python src/main.py --input ~/Workplace/BMME_new/conos/Whole_all/p21_2.lst --output emb/p2_test.emd



library(tidyverse)

p2=readRDS('Whole_all_p2combined.rds')


n2v_emb <- "/home/meisl/tools/node2vec/emb/p2.lst.emd" %>%
  read.table(skip=1) %>% arrange(V1) %>% select(-V1) %>% as.matrix()

n2v_tsne <- Rtsne::Rtsne(n2v_emb, verbose=T)$Y
rownames(n2v_tsne) <- paste(1:nrow(n2v_tsne))
rownames(n2v_tsne) <-names(p2$clusters$PCA$multilevel)


pdf('p2.pdf')
  conos::embeddingPlot(n2v_tsne, groups=anoCell)
  dev.off()


pdf('p2_samp.pdf')
  conos::embeddingPlot(n2v_tsne, groups=anoSample)
  dev.off()


pdf('p2_sampType.pdf')
  conos::embeddingPlot(n2v_tsne, groups=anoSample)
  dev.off()





 

n2v_emb <- "/home/meisl/tools/node2vec/emb/con_Whole_Tumor_new_t2_PCA.emd" %>%
  read.table(skip=1) %>% arrange(V1) %>% select(-V1) %>% as.matrix()

n2v_tsne <- Rtsne::Rtsne(n2v_emb, verbose=T)$Y
rownames(n2v_tsne) <- paste(1:nrow(n2v_tsne))
rownames(n2v_tsne) <-names(con$clusters$multi$groups)


pdf('con_PCA.pdf')
  conos::embeddingPlot(n2v_tsne, groups=anoCell,alpha=0.03,raster=T)
  dev.off()




 #con.gv.emb <- t(readRDS("/d0-mendel/home/viktor_petukhov/Copenhagen/GraphImprovement/embeddings/snap/examples/node2vec/emb/tm_tsne.rds"))

  pdf(file='g2v.pdf',width=8,height=8)
  con$embedding <- t(n2v_tsne);
  con$plotGraph(alpha=0.03,groups=anoCell,raster=T)  #groups=con$clusters$multi$groups,
  dev.off();




  pdf(file='g2v_samp.pdf',width=8,height=8)
  con$embedding <- t(n2v_tsne);
  con$plotGraph(alpha=0.03,groups=anoSample,raster=T)  #groups=con$clusters$multi$groups,
  dev.off();


  pdf(file='g2v_sampType.pdf',width=8,height=8)
  con$embedding <- t(n2v_tsne);
  con$plotGraph(alpha=0.03,groups=anoSampleType,raster=T)  #groups=con$clusters$multi$groups,
  dev.off();


