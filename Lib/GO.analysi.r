

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
  
  fo=paste(appname,'DAVID.rds',sep='')
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



# fin  different comparesion 
# nname=c('Tcyto','NK','TNK','Treg','T_help')
# path is result of DESeq result; 

gene_list=list()
fclist=list()

KEGG_lis=list()
BP_lis=list()
for ( i in seq(length(path))){
  tmp=readRDS(path[i])[[fin]][[1]][['res']]
  tmp=tmp[!is.na(tmp$log2FoldChange),]
  fc=tmp$log2FoldChange
  names(fc)=rownames(tmp)
  fc=fc[order(fc,decreasing=TRUE)]
  gene_list[[nname[i]]]=names(fc[1:50])
  fclist[[nname[i]]]=fc
  n1=paste(fin,'.',nname[i],sep='')
  
  
  KEGG_lis[[nname[i]]] <- tryCatch({
    DAVID_merge(n1,'KEGG_PATHWAY')
  },error=function(cond) {
    print('error')
    return(NULL)
  }
  )
 
  BP_lis[[nname[i]]] <- tryCatch({
    DAVID_merge(n1,'GOTERM_BP_DIRECT')
  },error=function(cond) {
    print('error')
    return(NULL)
  }
  )  
  
  
  pKEGG= listTomatrix_Value(KEGG_lis)
  pBP= listTomatrix_Value(BP_lis)
  
  pKEGG=as.matrix(pKEGG)
  pBP=as.matrix(pBP)
  
  rgb.palette <- colorRampPalette(c("white","red"), space = "rgb" )
  mi=min(pKEGG)
  ma=max(pKEGG)
  pdf=paste(nname[i],'.KEGG.merge.pdf',sep='')
  a=pheatmap(pKEGG,filename=pdf,color=rgb.palette(100),scale='none',border_color='NA',cluster_rows=T,cluster_cols=T,
             show_colnames=T,fontsize_col=6,fontsize_row=8,height=0.3*nrow(pKEGG),
             breaks = c(mi,seq(2,ma,length.out = 99)))
  
  
  rgb.palette <- colorRampPalette(c("white","red"), space = "rgb" )
  mi=min(pBP)
  ma=max(pBP)
  pdf=paste(nname[i],'.BP.merge.pdf',sep='')
  a=pheatmap(pBP,filename=pdf,color=rgb.palette(100),scale='none',border_color='NA',cluster_rows=T,cluster_cols=T,
             show_colnames=T,fontsize_col=6,fontsize_row=8,height=0.3*nrow(pBP),
             breaks = c(mi,seq(2,ma,length.out = 99)))
  
}




