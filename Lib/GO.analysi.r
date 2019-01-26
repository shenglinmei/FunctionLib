


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



DAVID_topG=function(fc,appname,dorder=FALSE){
  if (dorder){
    print('order by name')}else{
    fc <- fc[!is.na(fc)]
    fc=fc[order(fc,decreasing=T)] }
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



# bigM2 is big data matrix 
# all ano
#  group4 and atype; choosed 

prepare_rawList=function(bigM2,group4,atype,allano,ncut=15){
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
  
  t.bigM2=t.bigM2[rowSums(t.bigM2)>ncut,]
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
  
  r=list('dat_DEseq'=dat_DEseq,'dat_cpm'=dat_cpm,'raw'=raw.mats,'tab1'=tab1)
  
  fo=paste(appname,'.normlized.count.rds',sep='')
  saveRDS(r,fo)
  
  return(r)
}






runGexp=function(gs,appname,dat){
  dat_DEseq=dat[['dat_DEseq']]
  dat_cpm=dat[['dat_cpm']]
  #tab1=dat[['tab1']]
  
  tab1=apply(data.frame(colnames(dat_cpm)),1,function(x) strsplit(x,'-')[[1]][2])
  
  gs=intersect(rownames(dat_DEseq),gs)
  
  lrow=ceiling(length(gs)/5)
  
  lis=list()
  lis2=list()
  for (gene in gs) {
    library(ggplot2)
    
    
    nn=as.factor(tab1)
    nn <- ordered(nn, levels = c("Whole","Noninvolved", "Involved","Tumor"))
    table(nn==tab1)
    
    dat=data.frame('exp_DEseq'=dat_DEseq[gene,],'exp_cpm'=dat_cpm[gene,],'dtype'=nn)
    
    p <- ggplot(dat, aes(x=dtype,fill=dtype,y=exp_cpm)) + geom_boxplot(outlier.shape = NA)  + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.text.x=element_blank()) +xlab("Type") + ylab("Log Exp")
    p=p+ geom_point(aes(fill = dtype), size = 2, shape = 21, position = position_jitterdodge()) + ggtitle(gene) 
    
    lis[[gene]]=p
    
    
    
    p1 <- ggplot(dat, aes(x=dtype,fill=dtype,y=exp_DEseq)) + geom_boxplot(outlier.shape = NA)  + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.text.x=element_blank()) +xlab("Type") + ylab("Exp DESeq ")
    p1=p1+ geom_point(aes(fill = dtype), size = 2, shape = 21, position = position_jitterdodge()) + ggtitle(gene) 
    
    lis2[[gene]]=p1
    
  }
  
  # 
  #p <- ggplot(dat, aes(x=dtype,fill=dtype,y=exp_DEseq)) + geom_boxplot(outlier.shape = NA)  + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.text.x=element_blank()) +xlab("Type") + ylab("Exp")
  
  
  #p<- p+geom_point(aes(fill = atype), size = 0.6, shape = 21, position = position_jitterdodge())
  
  #ggsave(f2,plot=p,width=12,height=8)
  
  
  b=  cowplot::plot_grid(plotlist=lis, ncol=5, nrow=lrow)
  
  
  ggsave(paste(appname,'.png',sep=''),b,width = 5*5,height=3*lrow)
  
  
  
  
  
  b=  cowplot::plot_grid(plotlist=lis2, ncol=5, nrow=lrow)
  
  
  ggsave(paste(appname,'.DEseq.png',sep=''),b,width = 5*5,height=3*lrow)
  
  
  
}



GO_singleVec=function(fc,appname){
  
 # gsea=runGSEA(fc)
#  fo=paste(appname,'.gsea.new.rds',sep='')
#  saveRDS(gsea,fo)
  
  
  rlis=DAVID_topG(fc,appname)
  
  fo=paste(appname,'.DAVID.rds',sep='')
  saveRDS(rlis,fo)
  
  
  tmp=DAVID_merge(appname,'GOTERM_BP_DIRECT')
  
  tmp2=DAVID_merge(appname,'KEGG_PATHWAY')
  
}





