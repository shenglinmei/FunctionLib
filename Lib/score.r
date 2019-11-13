


n(p2myeloid,markers,appname){
  markers=intersect(markers,colnames(p2$counts))
  
  
  x <- as.matrix(p2myeloid$counts[names(conosCluster),markers])
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
  
  
  x=x[,o]
  #ss=apply(x,2,function(x) sum(x))
  #x=x[,abs(ss)>0.1]
  
  annot <- data.frame(CellType=conosCluster[colnames(x)],row.names = colnames(x))
  
  
  annot2=annot
  
  #annot2 = data.frame(ID = as.factor(as.character(annot[,1])))
  rownames(annot2)=colnames(x)
  
  
  cellA=annot2[,1]
  names(cellA)=rownames(annot2)
  
  o <- order(cellA)
  
  x=x[,o]
  
  pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
  ## draw heatmap
  
  fout=paste(appname,'.DiffG.png')
  rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
  heat=pheatmap(x,cluster_cols=FALSE,annotation_col = annot2,show_colnames = F,annotation_legend = TRUE, #,gaps_col =  1,
                cluster_rows = T,color=rgb.palette(100),filename=fout,fontsize_row =5,width=8,height=4*0.02*length(markers),
                breaks = c(seq(min(x),-0.01,length.out = 50),0.01,seq(0.1,2,length.out = 48),max(x)))
  
}



Boxplot_dat2=function(m2,p2,fraction,anoSample,group4){
  
  m2=intersect(m2,colnames(p2$counts))
  cname=names(group4)
  exp=as.matrix(p2$counts[cname,m2])
  
  #exp=apply(exp,2,function(x) x/max(x))
  
  m2score=rowMeans(exp)
  m2score_scale=scale(m2score)
  names(m2score_scale)=names(m2score)
  print(summary(m2score))
  dat=data.frame('m2score'=m2score,'cell'=group4,'Type'=fraction[cname],sample=anoSample[cname])
  nn=as.factor(dat$Type)
  
  
  dat$Type2=apply(dat,1,function(x) paste(x['cell'],x['sample']))
  tmp=tapply(dat$m2score,dat$Type2,mean)
  index=match(names(tmp),dat$Type2)
  tmp_type=dat$Type2[index]
  tmp_type[1:5]
  dat2=data.frame('Type2'=tmp_type,'m2score'=tmp,'Type'=dat$Type[index],'cell'=dat$cell[index],
                  'sample'=dat$sample[index])
  
  dat2$name=paste(rownames(dat2),as.character(dat2[,'cell']))
  print(table(dat2$Type))
  print(table(dat2$cell))
  return(dat)
  
}


drawBoxplot=function(appname,tmp,name2,clpalette.fraction=NULL,limHeight=1.5,dsize=3,dtype=NULL,dcell=NULL){
  # limHeight=1.45
  
  if (!is.null(dtype)){
    tmp=tmp[tmp$Type==dtype,]
    
  }
  if (!is.null(dcell)){
    tmp=tmp[tmp$cell==dcell,]
    tmp$cell2=tmp$cell
    tmp$cell=tmp$Type
  }
  
  t=tapply(tmp$m2score , tmp$cell,median)
  nn=names(t)[which.max(t)]
  
  # filter significant pairs 
  sig=compare_means(m2score ~ cell,  data = tmp)
  write.table(sig,paste(name2,'.pvalue.xls',sep=''),col.names=T,row.names=F,quote=F,sep='\t')
  sig=sig[sig$p.signif!='ns',]
  sig=sig[(sig$group1==nn | sig$group2==nn),]
  print(sig)
  siglis=split(sig, seq(nrow(sig)))
  pair=lapply(siglis,function(x) as.character(x[,2:3]))
  
  
  p1=ggboxplot(tmp, x = "cell", y = "m2score",fill ="cell",xlab = "",ylab=name2,width=0.6,
               color = "black",outlier.shape = NA)+
    stat_compare_means(comparisons = pair,label = "p.signif",hide.ns=TRUE,size=dsize,tip.lengt=0.01,bracket.size =0.3,label.y.npc = "bottom")  #+ # Add pairwise comparisons p-value
  p1=p1+ylim(c(min(tmp$m2score)*.7,max(tmp$m2score)*limHeight))
  p1=p1+ geom_point(data = tmp,aes(fill = cell), size = 0.5, shape = 21, position = position_jitterdodge()) # +scale_fill_manual(values=cols)
  p1=p1+ theme(legend.position="none")
  p1=p1+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p1
  
  
  ggsave(paste(appname,'',name2,'.score.pvalue.pdf',sep=''),p1,w=5,h=5)
  
  saveRDS(tmp,paste(appname,'',name2,'.dat.rds',sep=''))
  
  
  p1 <- ggplot(tmp, aes(x=cell,fill=cell,y=m2score)) + geom_boxplot(outlier.shape = -1,width=0.5,position=position_dodge(width=0.1)) +theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("")  + ylab(name2)
  p1=p1+ geom_point(data = tmp,aes(fill = cell), size = 0.3, shape = 21, position = position_jitterdodge()) # +scale_fill_manual(values=cols)
  p1=p1+ theme(legend.position="none")
  p1
  
  ggsave(paste(appname,'',name2,'.score.pdf',sep=''),p1,w=5,h=4)
  
  return(p1)
}









draw_heatmap=function(group4,p2,cols,name2,m2,fraction,dtype){
  ano=group4[fraction[names(group4)]==dtype ]
  
  marker=intersect(m2,colnames(p2$counts))
  
  dat=p2$counts[names(ano),marker]
  
  atype=fraction[names(ano)]    #  fraction
  
  ano=as.factor(ano)
  cells=levels(ano)
  
  res=NULL
  for( i in cells){
    index=ano==i
    tmp=dat[index,]
    ttype=atype[index]
    
    # if ('Involved' %in% ttype){
    #   ttype=as.factor(ttype)
    #   ttype <- ordered(ttype, levels = c( "Whole", "Noninvolved","Involved","Tumor"))
    
    # }
    
    tm=apply(tmp,2,function(x) tapply(x,ttype,mean))
    #rownames(tm)=paste(i,'.',rownames(tm),sep='')
    
    res=rbind(res,tm)
  }
  
  rownames(res)=cells
  
  x=res
  x <- apply(x, 2, function(xp) {
    qs <- quantile(xp,c(0.05,0.95))
    xp[xp<qs[1]] <- qs[1]
    xp[xp>qs[2]] <- qs[2]
    xp
  })
  x <- x[,(apply(x,2,sd) != 0)]
  x <- t(scale(x))
  
  
  annot4 = data.frame(CellType = unique(ano))
  rownames(annot4)=unique(ano)
  
  
  #aka3 = list(ID = cols) ,annotation_colors = aka3[1]
  fout=paste(paste(name2,'.gene.hetmap.pdf',sep=''),sep='')
  
  rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
  aa=pheatmap(x,annotation_col = annot4,show_rownames = T,show_colnames = T, fontsize_row =5,width=3,
              cluster_rows = TRUE, cluster_cols=FALSE,color=rgb.palette(100),filename=fout,height=3,
              breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
  
  
  #n1=aa$tree_col$order
  n2=aa$tree_row$order
  xx=x[n2,]
  
  
  rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
  aa=pheatmap(xx[,cells],annotation_legend = FALSE,annotation_col = annot4,annotation_colors = aka3[1],show_rownames = T,show_colnames = T, fontsize_row =7.2,width=3,
              cluster_rows = FALSE, cluster_cols=FALSE,color=rgb.palette(100),filename=fout,height=3+0.1*length(m2),
              breaks = c(seq(min(xx),-0.01,length.out = 50),seq(0.01,max(xx),length.out = 50)))
  
  
  fout=paste(paste(name2,'.gene.hetmap.pdf',sep=''),sep='')
  rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
  aa=pheatmap(xx[,cells],annotation_legend = FALSE,annotation_col = annot4,annotation_colors = aka3[1],show_rownames = T,show_colnames = T, fontsize_row =(5.3+nrow(xx)*0.1),width=2.7,
              cluster_rows = FALSE, cluster_cols=FALSE,color=rgb.palette(100),filename=fout,height=4.3,
              breaks = c(seq(min(xx),-0.01,length.out = 50),seq(0.01,max(xx),length.out = 50)))
  
  saveRDS(xx,paste(name2,'.dat.heatmap.rds',sep=''))
  return(xx)
}




getProportion=function(anoCell,anoSample,cname){
  ano2=data.frame('Cell'=anoCell[cname],'SampleType'=anoSample[cname])
  
  # Annotation vs sample
  tmp2 <- acast(ano2, Cell ~ SampleType, fun.aggregate=length)
  head(tmp2)
  # Normalise for the number of cells in each library
  tmp3 <- (sweep(tmp2, 2, colSums(tmp2), FUN='/'))
  tmp4 <- melt(tmp3)
  head(tmp4)
  names(tmp4) <- c('annot', 'sample','pc.of.sample')
  head(tmp4)
  return(tmp3)
}




getScore=function(p2,gs,nn_ano,kkey,atype,tsample){
  m2=intersect(gs,colnames(p2$counts))
  
  nn=nn_ano
  tmp=nn[nn %in% kkey]
  table(tmp)
  cname=names(tmp)
  
  cname=intersect(cname,rownames(p2$counts))
  
  # exp=p2$counts[cname,m2]
  
  m2score=rowMeans(as.matrix(p2$counts[cname,m2]))
  m2score_scale=scale(m2score)
  names(m2score_scale)=names(m2score)
  
  length(m2score_scale)
  
  
  dat=data.frame('m2score'=m2score_scale,'cell'=nn_ano[cname],'Type'=as.character(atype[cname]),sample=tsample[cname])
  
  
  
  dat$Type2=apply(dat,1,function(x) paste(x[2],x[3]))
  
  dat$Type3=apply(dat,1,function(x) paste(x[2],x[3],x[4]))
  
  tmp=tapply(dat$m2score,dat$Type3,mean)
  index=match(names(tmp),dat$Type3)
  tmp_type=dat$Type2[index]
  tmp_type[1:5]
  
  dat2=data.frame('Type2'=tmp_type,'m2score'=tmp,'Type'=dat$Type[index],'cell'=dat$cell[index],'sample'=dat$sample[index])
  print(table(dat2$Type2))
  return(dat2)
}




dat_heatmap=function(p2,m2,nn_ano,atype){
  
  m2=intersect(m2,colnames(p2$counts))
  
  nn=nn_ano
  dat=p2$counts[names(nn),m2]
  
  res=NULL
  for( i in unique(nn)){
    indexNmae=names(nn[nn==i])
    tmp=dat[indexNmae,]
    ttype=atype[indexNmae]
    print(i)
    print(table(ttype))
    #  if ('Involved' %in% ttype){
    ttype=as.factor(ttype)
    # ttype <- ordered(ttype, levels = c( "Whole", "Noninvolved","Involved","Tumor"))
    
    # }
    
    tm=apply(tmp,2,function(x) tapply(x,ttype,mean))
    
    rownames(tm)=paste(i,rownames(tm),sep='.')
    res=rbind(res,tm)
  }
  
  
  x=res
  x <- apply(x, 2, function(xp) {
    qs <- quantile(xp,c(0.05,0.95))
    xp[xp<qs[1]] <- qs[1]
    xp[xp>qs[2]] <- qs[2]
    xp
  })
  x <- x[,(apply(x,2,sd) != 0)]
  x <- t(scale(x))
  
  
  dim(x)
  return(x)
}





















#t=dat_heatmap(p2,m2,anoCell.merge2[c(cN,cT)],anoSampleType)

#r=getScore(p2,m2,anoCell.merge2[c(cN,cT)],c('Epitheial_Basal','Epitheial_Club'),anoSampleType,anoSample)


## M2 Score 
#gg='ARG1|ARG2|IL10|CD32|CD163|CD23|FCER2|CD200R1|PD-L2|PD-L1|MARCO|CSF1R|CD206|IL1RA|IL14R|CCL4|CCL13|CCL20|CCl17|CCL18|CCl22|CCL24|LYVE1|VEGFA|VEGFB|VEGFC|VEGFD|EGF|CTSA|CTSB|CTSC|CTSD|TGFB1|TGFB2|TGFB3|MMP14|MMP19|MMP9|CLEC7A|WNT7B|FASL|TNSF12|TNSF8CD276|VTCN1|MSR1|FN1|IRF4'
#name2='M2 score'
#m2=strsplit(gg,split='|', fixed=TRUE)[[1]]

#boxScore=Boxplot_dat2(m2,p2,anoSampleType,anoSample,anoCell.merge2[c(cN,cT)])

#drawBoxplot('Test',boxScore,name2,clpalette.fraction=NULL,dsize=4,dcell='Epitheial_Basal')

#t=draw_heatmap(anoCell.merge2[c(cN,cT)],p2,cols,name2,m2,anoSampleType,dtype='T')









### cytokine expression 
#cyto=readRDS('/home/meisl/bin/data/cytokine.rds')

#nfkb=read.csv('/home/meisl/Workplace/BMME/d.figure/NFKB.markerG.txt',sep='\t',header=F)

#cyto=as.character(nfkb[,1])


#cyto=intersect(cyto,colnames(p2$counts))

#exp=p2$counts[,cyto]

heatmapDat_sample=function(gs,p2,group4,anoSampleType,anoSample, dtype='N',cutoff=10){

 
  gs=intersect(gs,colnames(p2$counts))
  exp=p2$counts[,gs]
  
  ##
  
  cname=names(group4)
  
  dat=data.frame('cell'=group4,'Type'=anoSampleType[cname],sample=anoSample[cname])
  
  dat$Type2=apply(dat,1,function(x) paste(x['cell'],x['sample']))
  
  
  dat3=dat[as.character(dat[,'Type'])==dtype,]
  
  cname=rownames(dat3)
  
  exp=exp[cname,]
  
  dim(exp)
  dim(dat3)
  
  
  raw.mats=list()
  
  t.bigM2=t(exp)
  
  n.anoSample=dat3$Type2
  tab=table(n.anoSample)
  nname= names(tab[tab>cutoff])
  
  for (n in nname){
    cd=t.bigM2[,n.anoSample==n]
    n1=n
    raw.mats[[n1]]=cd
  }
  nn=unlist(lapply(raw.mats, function(x) ncol(x)))
  
  count1=do.call(cbind,lapply(raw.mats,function(x) rowMeans(as.matrix(x))))
  
  cn=apply(data.frame(colnames(count1)),1,function(x) strsplit(x,' ')[[1]][1])

  x=t(count1)
  x <- apply(x, 2, function(xp) {
    qs <- quantile(xp,c(0.05,0.95))
    xp[xp<qs[1]] <- qs[1]
    xp[xp>qs[2]] <- qs[2]
    xp
  })
  x <- x[,(apply(x,2,sd) != 0)]
  x <- t(scale(x))
  
  x1=x
  fout=paste('inflamatory.signature.heatmap2.pdf',sep='')

  rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
  pp=pheatmap(x1,,show_rownames = T,show_colnames = T, fontsize_row =6,width=8,
              cluster_rows = TRUE, cluster_cols=FALSE,color=rgb.palette(100),filename=fout,height=13,
              breaks = c(seq(min(x1),-0.01,length.out = 50),seq(0.01,max(x1),length.out = 50)))
  print(dim(x))
  return(x1)
}


#x=heatmapDat_sample(cyto,p2,group4,anoSampleType,anoSample, dtype='N',cutoff=10)
