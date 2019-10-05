
GOanalysis.term2=function(markers,n,bg){
  
  ENTREZID=unlist(mget(markers, org.Hs.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID=ENTREZID[!is.na(ENTREZID)]
  
  Ebg=unlist(mget(bg, org.Hs.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID2=Ebg[!is.na(Ebg)]
  
  
  allr=list()
  
  for(function_type in c("BP", "CC", "MF")){
    
    param <- new("GOHyperGParams", geneIds=ENTREZID,
                 universeGeneIds=ENTREZID2,
                 annotation="org.Hs.eg.db", ontology=function_type,pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp <- hyperGTest(param)
    sumTable <- summary(hyp)
    
    
    
    david=sumTable[1:20,]
    david$Pvalue=-log(david[,2])
    termNumber=nrow(david)
    
    david$genes=apply(david,1,function(x) { paste(names(ENTREZID[ENTREZID %in% get(as.character(x[1]),org.Hs.egGO2ALLEGS)]),collapse = ',' ) } )
    
    allr[[function_type]]=david
    write.table(david,paste(n,'_' ,function_type, ".xls", sep=""),sep='\t',col.names=T,row.names=F,quote=F)
    
    
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
  
  saveRDS(allr,paste(n,'.GOterm.rds',sep=''))
  return(allr)
}








DEcaculate3.term=function(p2,appname,conosCluster,removeGene=NULL,cutoff=3,num=100,GO=NULL,bdg=NULL){
  allg=NULL  
  de1 <- p2$getDifferentialGenes(groups=conosCluster,z.threshold = cutoff)
  
  folderN=paste(appname,'.DE',sep='')
  pwd=getwd()
  pwd2=paste(pwd,'/',folderN,'/',sep='')
  
  system(paste('mkdir ',folderN))
  
  setwd(pwd2)
  
  
  
  
  f1=paste(appname,'_diffGene.rds',sep='')
  saveRDS(de1,f1)
  for(n in names(de1)){
    x=de1[[n]]
    z <- x[order(-x$Z),]
    z <- z[z$Z>0,]
   # z <- z[z$fe>0.4,]
    if(!is.null(removeGene)){
      index=grepl('^RP[LKS]',rownames(z))
      z=z[!index,]
    }
    markers=rownames(z)[1:num]
    markers=markers[!is.na(markers)]
    allg=cbind(allg,markers)
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
    o <- order(conosCluster[colnames(x)])
    annot <- data.frame(cluster=conosCluster[colnames(x)],row.names = colnames(x))
    
    annot$dtype='other'
    annot[as.character(annot[,1])==n,'dtype']=n
    annot$dtype=as.factor(annot$dtype)
    
    
    pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
    ## draw heatmap
    
    fout=paste(appname,'_',n,'_marker.heatmap.png',sep='')
    rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
    pheatmap(x[,o],labels_col=FALSE,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,
             cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=5,height=8,   #3.5*0.02*length(markers),
             breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
    
    
    iid= paste(appname,'_cluster_',n,sep='')
    if(!is.null(GO)){
      
      figGO=GOanalysis.term2(markers,iid,bdg) 
      figGO=figGO[['BP']]
      
      for (jj in seq(nrow(figGO))){
        
        gs=figGO[jj,'genes']
        gs=strsplit(gs,',')[[1]]
        
        if (length(gs) >= 4 ){
          
          fout=paste(appname,'_',n,'_',figGO[jj,'Term'],'_BP.heatmap.png',sep='')
          rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
          pheatmap(x[gs,o],labels_col=FALSE,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,
                   cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=5,height=5*0.02*length(gs),
                   breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
          
        }}
      
      
    }
  }
  colnames(allg)=names(de1)
  
  
  write.table(allg,paste(appname,'_Diff_gene.xls',sep=''),sep='\t',col.names=T,row.names=F,quote=F)
  setwd(pwd)
}

tmp=Toch(group4)
tmp[tmp=='Mono-2']='Mono-1'

de1 <- p2$getDifferentialGenes(groups=tmp,z.threshold = 0)

de1 <- p2$getDifferentialGenes(groups=group4[group4 %in% c('Mono-1','Mono-2') ],z.threshold = 5)


names(de1)

tmp2=de1$`Mono-1`
tmp2=tmp2[tmp2$Z>0,]
dim(tmp2)
bdg=rownames(tmp2)


exp=rowMeans(as.matrix(p2$counts[names(group4[group4 %in% c('Mono-1','Mono-2')]),]))
summary(exp)
bdg=names(exp[exp>0.017])
length(bdg)



setwd('/d0/home/meisl/Workplace/d.figure/new_F2/version2/')
DEcaculate3.term(p2,'Mono_top100.bdg3',group4[group4 %in% c('Mono-1','Mono-2') ],removeGene=TRUE,cutoff=2,num=200,GO=TRUE,bdg=bdg)


de1 <- p2$getDifferentialGenes(groups=group4[group4 %in% c('Mono-1','Mono-2') ],z.threshold = 0)
tmp2=de1$`Mono-2`
tmp2=tmp2[order(tmp2$Z,decreasing=T),]

fc=tmp2$Z
names(fc)=rownames(tmp2)
#fc=fc[!is.na(fc)]

#fc=fc[is.finite(fc)]


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

gsea=runGSEA(fc)

saveRDS(gsea,'momo1.rds')




library(fgsea)
fc=fc[fc!=0]

fgseaRes <- fgsea(msigdb.c5, fc, minSize=10, maxSize = 500, nperm=1000)

head(fgseaRes[order(padj, -abs(NES)), ], n=10)


plotEnrichment(msigdb.c5[["GO_KERATIN_FILAMENT"]], fc)

barplot(sort(fc, decreasing = T))



matrices2 <- do.call(c,mapply(function(m, name) {
 rep(name,length(m))
},
msigdb.h,
names(msigdb.h)))


matrices3 <- do.call(c,mapply(function(m, name) {
  m
},
msigdb.h,
names(msigdb.h)))


d1=data.frame(matrices2,matrices3)
d2=unique(data.frame(matrices2,matrices2))
USER_DATA <- build_Anno(d1,d2)

r=enricher_internal(gene = markers, pvalueCutoff = 0.1, 
                     pAdjustMethod = "BH", universe = universe, minGSSize = 10, 
                     maxGSSize = 500, qvalueCutoff = 0.1, USER_DATA = USER_DATA)










library(clusterProfiler)
markers=names(fc)[1:100]

ENTREZID=unlist(mget(markers, org.Hs.egSYMBOL2EG, ifnotfound=NA))
ENTREZID=ENTREZID[!is.na(ENTREZID)]

xx <- enrichMKEGG(ENTREZID, organism='hsa', minGSSize=1)


david = enrichDAVID(gene = ENTREZID, idType="ENTREZ_GENE_ID", 
                    listType="Gene", annotation="KEGG_PATHWAY")

library(RDAVIDWebService)

data(geneList)
gene = names(geneList)[abs(geneList) > 2]
david = enrichDAVID(gene = gene, idType="ENTREZ_GENE_ID", 
                     annotation="KEGG_PATHWAY")



david<-DAVIDWebService$new(email="ma.tongji@gmail.com")

enricher_internal <- DOSE:::enricher_internal
build_Anno <- DOSE:::build_Anno
EXTID2TERMID=DOSE:::EXTID2TERMID

tt=data.frame( 'b'=c('f1','f1','f2','f2'),'a'=c('AR','ESR1','CD4','CD8A'))
tt2=data.frame( 'b'=c('f1','f1','f2','f2'),'a'=c('f1','f1','f2','f2'))

USER_DATA <- build_Anno(tt,tt2)






# disease2gene=gda[, c("diseaseId", "geneId")]
# disease2name=gda[, c("diseaseId", "diseaseName")]


r=enricher_internal(gene = c('AR','ESR1'), pvalueCutoff = 0.1, 
                  pAdjustMethod = "BH", universe = universe, minGSSize = 10, 
                  maxGSSize = 500, qvalueCutoff = 0.1, USER_DATA = USER_DATA)


assign("PATHID2EXTID", msigdb.h, envir = Anno_clusterProfiler_Env)
assign("EXTID2PATHID", EXTID2PATHID, envir = Anno_clusterProfiler_Env)








x = enricher(markers, TERM2GENE=msigdb.c5[1:5], TERM2NAME=names(msigdb.c5)[1:5])

enricher_internal2=function (gene, pvalueCutoff, pAdjustMethod = "BH", universe = NULL, 
          minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, USER_DATA) 
{
  EXTID2TERMID=DOSE:::EXTID2TERMID
  
  gene=markers
  
  gene <- as.character(unique(gene))
  qExtID2TermID <- EXTID2TERMID(gene, USER_DATA)
  qTermID <- unlist(qExtID2TermID)
  print(qExtID2TermID)
  print(qTermID)
  if (is.null(qTermID)) {
    message("--> No gene can be mapped....")
    p2e <- get("PATHID2EXTID", envir = USER_DATA)
    sg <- unlist(p2e[1:10])
    sg <- sample(sg, min(length(sg), 6))
    message("--> Expected input gene ID: ", paste0(sg, collapse = ","))
    message("--> return NULL...")
    return(NULL)
  }
  qExtID2TermID.df <- data.frame(extID = rep(names(qExtID2TermID), 
                                             times = lapply(qExtID2TermID, length)), termID = qTermID)
  qExtID2TermID.df <- unique(qExtID2TermID.df)
  qTermID2ExtID <- with(qExtID2TermID.df, split(as.character(extID), 
                                                as.character(termID)))
  extID <- DOSE:::ALLEXTID(USER_DATA)
  if (missing(universe)) 
    universe <- NULL
  if (!is.null(universe)) {
    extID <- intersect(extID, universe)
  }
  qTermID2ExtID <- lapply(qTermID2ExtID, intersect, extID)
  qTermID <- unique(names(qTermID2ExtID))
  termID2ExtID <- DOSE:::TERMID2EXTID(qTermID, USER_DATA)
  termID2ExtID <- lapply(termID2ExtID, intersect, extID)
  geneSets <- termID2ExtID
  idx <- DOSE:::get_geneSet_index(termID2ExtID, minGSSize, maxGSSize)
  if (sum(idx) == 0) {
    msg <- paste("No gene set have size >", minGSSize, "...")
    message(msg)
    message("--> return NULL...")
    return(NULL)
  }
  termID2ExtID <- termID2ExtID[idx]
  qTermID2ExtID <- qTermID2ExtID[idx]
  qTermID <- unique(names(qTermID2ExtID))
  k <- sapply(qTermID2ExtID, length)
  k <- k[qTermID]
  M <- sapply(termID2ExtID, length)
  M <- M[qTermID]
  N <- rep(length(extID), length(M))
  n <- rep(length(qExtID2TermID), length(M))
  args.df <- data.frame(numWdrawn = k - 1, numW = M, numB = N - 
                          M, numDrawn = n)
  pvalues <- apply(args.df, 1, function(n) phyper(n[1], n[2], 
                                                  n[3], n[4], lower.tail = FALSE))
  GeneRatio <- apply(data.frame(a = k, b = n), 1, function(x) paste(x[1], 
                                                                    "/", x[2], sep = "", collapse = ""))
  BgRatio <- apply(data.frame(a = M, b = N), 1, function(x) paste(x[1], 
                                                                  "/", x[2], sep = "", collapse = ""))
  Over <- data.frame(ID = as.character(qTermID), GeneRatio = GeneRatio, 
                     BgRatio = BgRatio, pvalue = pvalues, stringsAsFactors = FALSE)
  p.adj <- p.adjust(Over$pvalue, method = pAdjustMethod)
  qobj <- tryCatch(qvalue(p = Over$pvalue, lambda = 0.05, pi0.method = "bootstrap"), 
                   error = function(e) NULL)
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  }else {
    qvalues <- NA
  }
  geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse = "/"))
  geneID <- geneID[qTermID]
  Over <- data.frame(Over, p.adjust = p.adj, qvalue = qvalues, 
                     geneID = geneID, Count = k, stringsAsFactors = FALSE)
  Description <- DOSE:::TERM2NAME(qTermID, USER_DATA)
  if (length(qTermID) != length(Description)) {
    idx <- qTermID %in% names(Description)
    Over <- Over[idx, ]
  }
  Over$Description <- Description
  nc <- ncol(Over)
  Over <- Over[, c(1, nc, 2:(nc - 1))]
  Over <- Over[order(pvalues), ]
  Over$ID <- as.character(Over$ID)
  Over$Description <- as.character(Over$Description)
  row.names(Over) <- as.character(Over$ID)
  x <- new("enrichResult", result = Over, pvalueCutoff = pvalueCutoff, 
           pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff, 
           gene = as.character(gene), universe = extID, geneSets = geneSets, 
           organism = "UNKNOWN", keytype = "UNKNOWN", ontology = "UNKNOWN", 
           readable = FALSE)
  return(x)
}

r=enricher_internal2(gene = markers, pvalueCutoff = 0.1, 
                    pAdjustMethod = "BH", universe = universe, minGSSize = 10, 
                    maxGSSize = 500, qvalueCutoff = 0.1, USER_DATA = USER_DATA)



enricher_internal <- DOSE:::EXTID2TERMID(c('AR'), USER_DATA)


