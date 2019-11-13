prepare_rawList3=function(expA,bigM,group4,atype,anoSample,fractionf,ncut=2){
  
  
  #group4=ano
  #atype='Mono-1'
  #anoSample=samplef
  #ncut=2
  #fractionf=fractionf
  
  gname=names(group4[group4==atype])
  print(length(gname))
  
  tab=table(anoSample[gname])
  nnames=names(tab[tab>ncut])
  
  raw.mats=list()
  
  n.anoSample=anoSample[gname]
  for (n in nnames){
    cd=bigM[,names(n.anoSample[n.anoSample==n])]
    n1=n
    raw.mats[[n1]]=cd
  }
  
  nn=unlist(lapply(raw.mats, function(x) ncol(x)))
  print(atype)
  print(nn)
  
  cname=do.call(c,lapply(raw.mats,function(x) colnames(x)))
  expMatrix=t(p2$counts[cname,])
  
  
  count1=do.call(cbind,lapply(raw.mats,function(x) rowSums(as.matrix(x))))
  
  cm.tmp <- as.matrix(count1)
  rownames(cm.tmp) <- rownames(cm)
  ## calculate cpm
  cpm <- sweep(cm.tmp, 2, apply(cm.tmp,2, sum), FUN='/')
  dat_cpm <- log10(cpm * 1e6 + 1)
  rownames(dat_cpm)=rownames(count1)
  
  
  raw.exp=lapply(raw.mats,function(x) expA[,colnames(x)])
  unlist(lapply(raw.exp, function(x) ncol(x)))
  exp=do.call(cbind,lapply(raw.exp,function(x) rowMeans(as.matrix(x))))
  
  fraction=apply(data.frame(colnames(exp)),1, function(x) strsplit(x[1],'-')[[1]][2])
  
  res=list('cpm'=dat_cpm,
           'Mexp'=exp,
           'CellNum'=nn,
           'raw.mats'=raw.mats,
           'fraction'=fraction,
           'cSample'=anoSample[cname],
           'cfraction'=fractionf[cname],
           'expMatrix'=expMatrix,
           'zscore'=de1[[atype]])
  
  return(res)
  
}


sn=function(x) { names(x) <- x; return(x); }




caculateCor=function(tRL,m1,m2,g1,g2,inter,method='pearson'){
  corr=apply(tRL,1,function(x) {
    
    g1=as.character(unlist(x['gs1']))
    g2=as.character(unlist(x['gs2']))
    
    cor(m1[g1,inter],m2[g2,inter],method=method)
  })
  
  
  if (length(inter)>4){
    
    pvalue=apply(tRL,1,function(x) {
      
      g1=as.character(unlist(x['gs1']))
      g2=as.character(unlist(x['gs2']))
      
      cor.test(m1[g1,inter],m2[g2,inter],method=method)$p.value
    })
  }else{
    pvalue=0
  }
  
  return(list('corr'=corr,'pvalue'=pvalue))
}



getInfor=function(cell1,cell2,dat,gs1,gs2){
  print('Start')
  print(cell1)
  print(cell2)
  
  e1n=colnames(dat[[cell1]]$expMatrix[gs1,])
  e2n=colnames(dat[[cell2]]$expMatrix[gs2,])
  
  
  m1=dat[[cell1]]$cpm
  m2=dat[[cell2]]$cpm
  inter=intersect(colnames(m1),colnames(m2))
  
  
  e1n=names(anoSample[e1n][anoSample[e1n] %in% inter] )
  e2n=names(anoSample[e2n][anoSample[e2n] %in% inter] )
  
  
  e1=bigM[gs1,e1n]
  e2=bigM[gs2,e2n] 
  
  
  e1r=apply(e1,1,function(x) length(x[x>0]))/ncol(e1)
  e2r=apply(e2,1,function(x) length(x[x>0]))/ncol(e2)
  
  e2r[1:12]
  e1r[1:12]
  
  
  
  
  
  
  Lexp=rowMeans(dat[[cell1]]$Mexp[gs1,inter])
  Rexp=rowMeans(dat[[cell2]]$Mexp[gs2,inter])
  
  Lcpm=rowMeans(dat[[cell1]]$cpm[gs1,inter])
  Rcpm=rowMeans(dat[[cell2]]$cpm[gs2,inter])
  
  
  Lzscore=dat[[cell1]]$zscore[gs1,'Z']
  Rzscore=dat[[cell2]]$zscore[gs2,'Z']
  
  tRL=data.frame('gs1'=gs1,'gs2'=gs2)
  
  corr_pearson=caculateCor(tRL,m1,m2,g1,g2,inter,method='pearson')
  corr_spearman=caculateCor(tRL,m1,m2,g1,g2,inter,method='spearman')
  
  tmp=data.frame('Gene.a'=gs1,'Gene.b'=gs2,'MeanExpression.a'=Lexp,'MeanExpression.b'=Rexp,
                 'exp.ratio.a'=e1r, 'exp.ratio.b'=e2r,
                                 'zscore.a'=Lzscore,'zscore.b'=Rzscore,
                 'M.CPM.a'=Lcpm,'M.CPM.b'=Rcpm,
                 'corr.pearson'=corr_pearson$corr,'pvalue.pearson'=corr_pearson$pvalue,
                 'corr.spearman'=corr_spearman$corr,'pvalue.spearman'=corr_spearman$pvalue,'l'=length(inter))
  
  
  
  return(tmp)
}





Outputf=function(lis,fl1,fl2,outN,cut.ratio.a=0.1,cut.ratio.b=0.1){

    #outN='Res2'
   # cut.ratio.a=0.01
    #cut.ratio.b=0.01
    
    folderN=outN
    pwd=getwd()
    pwd2=paste(pwd,'/',folderN,'/',sep='')
    
    system(paste('mkdir ',folderN))
    
    setwd(pwd2)
  
  lc1=length(fl1)
  lc2=length(fl2)
  
  stat=matrix(rep(NA,lc1*lc2),lc1,lc2)
  rownames(stat)=fl1
  colnames(stat)=fl2
  
  stat2=stat
  LL=stat
  RR=stat
  
  
  flis=list()
  flis2=list()
  
  for(i in fl1){
    for (j in fl2){
      
      cell1=i
      cell2=j
      
      if (cell1!=cell2){
      
      fname=paste(cell1,cell2,sep='_')
      tmp=lis[[fname]]
      print(cell1)
      print(cell2)
      
      
      
      t1=tmp
      print(t1[1,'l'])
      
      # exp.ratio.a> 0.1 & exp.ratio.b > 0.1 &
      if (t1[1,'l']>5){
        
        
        t2=subset(t1 , MeanExpression.a >0 & MeanExpression.b >0 & 'M.CPM.a'>0.01 & 'M.CPM.b'>0.01)
        t2=subset(t2 ,exp.ratio.a> 0.01 & exp.ratio.b > 0.01 )
        
        t2=subset(t2 ,exp.ratio.a> cut.ratio.a | exp.ratio.b > cut.ratio.b )
        
        t2=t2[!is.na(t2$pvalue.spearman),]
        
        dim(t2)
        
        
        t3=subset(t2 , !is.na(corr.spearman) & !is.na(pvalue.spearman) )
        t3=subset(t3 , pvalue.spearman < 0.05)
        t3=subset(t3 , corr.spearman >0)
        
        
        
        #  t3=t3[,c(1,2,3,4,7,8,9,10,11)]
        # t3=t2
        num=nrow(t2)
        flis[[fname]]=t3
        flis2[[fname]]=t2
        num2=nrow(t3)
        
        Lr=length(unique(t2[,1]))/length(unique(RL[,1]))
        Rr=length(unique(t2[,2]))/length(unique(RL[,2]))
      }else{
        num=NA
        num2=NA
        Lr=NA
        Rr=NA
      }
      
      stat[i,j]=num
      stat2[i,j]=num2
      LL[i,j]=Lr
      RR[i,j]=Rr
      
    }
    }
  }
  
  
  fout=function(tmp){
    for( i in seq(ncol(tmp)))
      if( is.numeric(tmp[1,i])){
        tmp[,i]=round(tmp[,i],4)
      }
    return(tmp)
  }
  
  
  
  rgb.palette <- colorRampPalette(c("white","red"), space = "rgb" )
  a1=pheatmap(stat,show_rownames = T,width=10,height=8,display_numbers = stat,
              cluster_rows = F, cluster_cols=F,color=rgb.palette(100),filename=paste(outN,'.num.pdf',sep=''),number_color='black')
  
  #breaks = c(seq(0,quantile(stat,0.1),length.out = 2),seq(quantile(stat,0.1)+1,max(stat),length.out = 98)))
  
  
  rgb.palette <- colorRampPalette(c("white","red"), space = "rgb" )
  a1=pheatmap(stat2,show_rownames = T,width=10,height=8,display_numbers = stat2,
              cluster_rows = F, cluster_cols=F,color=rgb.palette(100),filename=paste(outN,'.num2.pdf',sep=''),number_color='black')
  
  
  
  LL=round(LL,3)
  rgb.palette <- colorRampPalette(c("white","red"), space = "rgb" )
  a2=pheatmap(LL,show_rownames = T,width=5,height=4,display_numbers = LL,
              cluster_rows = F, cluster_cols=F,color=rgb.palette(100),filename=paste(outN,'.LigandRatio.pdf',sep=''),number_color='black')
  
  
  
  
  RR=round(RR,3)
  rgb.palette <- colorRampPalette(c("white","red"), space = "rgb" )
  a3=pheatmap(LL,show_rownames = T,width=5,height=4,display_numbers = RR,
              cluster_rows = F, cluster_cols=F,color=rgb.palette(100),filename=paste(outN,'.ReceptorRatio.pdf',sep=''),number_color='black')
  
  
  
  flis2=flis

  for( fname in names(flis2)){
    print(dim((flis2[[fname]])))
    write.csv(fout(flis2[[fname]]),paste(fname,'.csv',sep=''),col.names=T,row.names=F,quote=F,sep=',' )
  }
  
  
  
  
  
  cell.types=fl1
  
  comp.prefix <- fl2
  
  names(cell.types) <- cell.types
  names(comp.prefix) <- comp.prefix
  
  
  
  corr.prefix='Filter'
  names(corr.prefix)=corr.prefix
  
  
  ## prefix for local files
  prefix.local = './'
  ## prefix for links to data
  ## url of the viewer
  viewerprefix <- paste('http://pklab.med.harvard.edu/shenglin/test/LigandReceptor/',outN,'/index?file=',sep='')
  
  x <- melt(lapply(comp.prefix, function(comp) {
    lapply(cell.types, function(celltype) {
      lapply(corr.prefix, function(corr) {
        file <- as.character(paste0( celltype,'_',comp ,'.csv',sep=''))
        file
      })
    })
  }))
  x$exists <- file.exists(paste0(prefix.local,'/',as.character(x$value)))
  
  library(xtable)
  
  x$num=apply(x,1,function(x) stat2[match(x['L2'],rownames(stat2)),match(x['L1'],colnames(stat2))])
  
  
  ## Table for non-corrected
  x.uncorr <- subset(x, x$L3 == 'Filter')
  x.uncorr$link <- paste0(viewerprefix, x.uncorr$value)
  x.uncorr$html <- ifelse(x.uncorr$exists,paste0('<a href="',x.uncorr$link,'">',x.uncorr$num, '</a>'),'N/A')
  x.uncorr.html <- print(xtable(dcast(x.uncorr, L2 ~ L1, value.var='html')),type='html',sanitize.text.function=I)
  
  
  write(x.uncorr.html, 'table.html')
  
  setwd(pwd)
  
  tmp=list('a'=a1,'a1'=a2,'a2'=a3)
  return(tmp)
}






#tmp1=Outputf(lis,fl1,fl2,'Myeloid.To.Tcell',cut.ratio.a=0.01,cut.ratio.b=0.01)
#tmp2=Outputf(lis,fl2,fl1,'Tcell.To.Myeloid',cut.ratio.a=0.01,cut.ratio.b=0.01)






Outputf=function(lis,fl1,fl2,outN,cut.ratio.a=0.01,cut.ratio.b=0.01){
  
  #outN='Res2'
  # cut.ratio.a=0.01
  #cut.ratio.b=0.01
  
  
  
  folderN=outN
  pwd=getwd()
  pwd2=paste(pwd,'/',folderN,'/',sep='')
  
  system(paste('mkdir ',folderN))
  
  setwd(pwd2)
  
  lc1=length(fl1)
  lc2=length(fl2)
  
  stat=matrix(rep(NA,lc1*lc2),lc1,lc2)
  rownames(stat)=fl1
  colnames(stat)=fl2
  
  stat2=stat
  LL=stat
  RR=stat
  
  
  flis=list()
  flis2=list()
  
  for(i in fl1){
    for (j in fl2){
      
      cell1=i
      cell2=j
      
      if (cell1!=cell2){
        
        fname=paste(cell1,cell2,sep='_')
        tmp=lis[[fname]]
        print(cell1)
        print(cell2)
        
        
        
        t1=tmp
        print(t1[1,'l'])
        
        # exp.ratio.a> 0.1 & exp.ratio.b > 0.1 &
        if (t1[1,'l']>5){
          
          
          t2=subset(t1 , MeanExpression.a >0 & MeanExpression.b >0 & 'M.CPM.a'>0.01 & 'M.CPM.b'>0.01)
          t2=subset(t2 ,exp.ratio.a> 0.01 & exp.ratio.b > 0.01 )
          
          t2=subset(t2 ,exp.ratio.a> cut.ratio.a | exp.ratio.b > cut.ratio.b )
          
          t2=t2[!is.na(t2$pvalue.spearman),]
          
          dim(t2)
          
          
          t3=subset(t2 , !is.na(corr.spearman) & !is.na(pvalue.spearman) )
          t3=subset(t3 , pvalue.spearman < 0.05)
          t3=subset(t3 , corr.spearman >0)
          
          #  t3=t3[,c(1,2,3,4,7,8,9,10,11)]
          # t3=t2
          num=nrow(t2)
          flis[[fname]]=t3
          flis2[[fname]]=t2
          num2=nrow(t3)
          
          Lr=length(unique(t2[,1]))/length(unique(RL[,1]))
          Rr=length(unique(t2[,2]))/length(unique(RL[,2]))
        }else{
          num=NA
          num2=NA
          Lr=NA
          Rr=NA
        }
        
        stat[i,j]=num
        stat2[i,j]=num2
        LL[i,j]=Lr
        RR[i,j]=Rr
        
      }
    }
  }
  
  
  fout=function(tmp){
    for( i in seq(ncol(tmp)))
      if( is.numeric(tmp[1,i])){
        tmp[,i]=round(tmp[,i],4)
      }
    return(tmp)
  }
  
  
  
  rgb.palette <- colorRampPalette(c("white","red"), space = "rgb" )
  a1=pheatmap(stat,show_rownames = T,width=10,height=8,display_numbers = stat,
              cluster_rows = F, cluster_cols=F,color=rgb.palette(100),filename=paste(outN,'.num.pdf',sep=''),number_color='black')
  
  #breaks = c(seq(0,quantile(stat,0.1),length.out = 2),seq(quantile(stat,0.1)+1,max(stat),length.out = 98)))
  
  
  rgb.palette <- colorRampPalette(c("white","red"), space = "rgb" )
  a1=pheatmap(stat2,show_rownames = T,width=10,height=8,display_numbers = stat2,
              cluster_rows = F, cluster_cols=F,color=rgb.palette(100),filename=paste(outN,'.num2.pdf',sep=''),number_color='black')
  
  
  
  LL=round(LL,3)
  rgb.palette <- colorRampPalette(c("white","red"), space = "rgb" )
  a2=pheatmap(LL,show_rownames = T,width=5,height=4,display_numbers = LL,
              cluster_rows = F, cluster_cols=F,color=rgb.palette(100),filename=paste(outN,'.LigandRatio.pdf',sep=''),number_color='black')
  
  
  
  
  RR=round(RR,3)
  rgb.palette <- colorRampPalette(c("white","red"), space = "rgb" )
  a3=pheatmap(LL,show_rownames = T,width=5,height=4,display_numbers = RR,
              cluster_rows = F, cluster_cols=F,color=rgb.palette(100),filename=paste(outN,'.ReceptorRatio.pdf',sep=''),number_color='black')
  
  
  
  
  stat2=stat
  
  for( fname in names(flis2)){
    print(dim((flis2[[fname]])))
    write.csv(fout(flis2[[fname]]),paste(fname,'.csv',sep=''),col.names=T,row.names=F,quote=F,sep=',' )
  }
  
  
  
  
  
  cell.types=fl1
  
  comp.prefix <- fl2
  
  names(cell.types) <- cell.types
  names(comp.prefix) <- comp.prefix
  
  
  
  corr.prefix='Filter'
  names(corr.prefix)=corr.prefix
  
  
  ## prefix for local files
  prefix.local = './'
  ## prefix for links to data
  ## url of the viewer
  viewerprefix <- paste('http://pklab.med.harvard.edu/shenglin/test/LigandReceptor/',outN,'/index?file=',sep='')
  
  x <- melt(lapply(comp.prefix, function(comp) {
    lapply(cell.types, function(celltype) {
      lapply(corr.prefix, function(corr) {
        file <- as.character(paste0( celltype,'_',comp ,'.csv',sep=''))
        file
      })
    })
  }))
  x$exists <- file.exists(paste0(prefix.local,'/',as.character(x$value)))
  
  library(xtable)
  
  x$num=apply(x,1,function(x) stat2[match(x['L2'],rownames(stat2)),match(x['L1'],colnames(stat2))])
  
  
  ## Table for non-corrected
  x.uncorr <- subset(x, x$L3 == 'Filter')
  x.uncorr$link <- paste0(viewerprefix, x.uncorr$value)
  x.uncorr$html <- ifelse(x.uncorr$exists,paste0('<a href="',x.uncorr$link,'">',x.uncorr$num, '</a>'),'N/A')
  x.uncorr.html <- print(xtable(dcast(x.uncorr, L2 ~ L1, value.var='html')),type='html',sanitize.text.function=I)
  
  
  write(x.uncorr.html, 'table.html')
  
  setwd(pwd)
  
  tmp=list('a'=a1,'a1'=a2,'a2'=a3)
  return(tmp)
}






