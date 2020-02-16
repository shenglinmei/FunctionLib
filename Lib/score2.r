

scoreData=function(gs,exp,fraction,anoSample,groups,cn1=NULL){
    
    cname=names(groups)
    cname=intersect(cname,colnames(exp))
    
    if (!is.null(cn1)){
      cname=intersect(cname,cn1)
    }
    
    print(table(groups[cname]))
    print(table(fraction[cname]))
    
    gs=intersect(gs,rownames(exp))
    
    t.exp=t(as.matrix(exp[gs,cname]))
    #t.exp=apply(t.exp,2,function(x) x/max(x))
    m2score=rowMeans(t.exp)
    m2score_scale=scale(m2score)
    names(m2score_scale)=names(m2score)
    print(summary(m2score))
    dat=data.frame('m2score'=m2score,'cell'=groups[cname],'Type'=fraction[cname],sample=anoSample[cname])
    #nn=as.factor(dat$Type)
    dat$Type2=apply(dat,1,function(x) paste(x['cell'],x['sample']))
    tmp=tapply(dat$m2score,dat$Type2,mean)
    index=match(names(tmp),dat$Type2)
    tmp_type=dat$Type2[index]
    tmp_type[1:5]
    dat2=data.frame('Type2'=tmp_type,'m2score'=tmp,'Type'=dat$Type[index],'cell'=dat$cell[index],
                    'sample'=dat$sample[index])
    
    dat2$name=paste(rownames(dat2),as.character(dat2[,'cell']))
  return(dat2)
}
    


figbox=function(tmp,name2,fname){
  p1 <- ggplot(tmp, aes(x=Type,fill=Type,y=m2score)) + geom_boxplot(outlier.shape = -1,width=0.5,position=position_dodge(width=0.1)) +theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("")  + ylab(name2)
  p1=p1+ geom_point(data = tmp,aes(fill = Type), size = 0.3, shape = 21, position = position_jitterdodge()) # +scale_fill_manual(values=cols)
  p1=p1+ theme(legend.position="none")
  
  fout=paste(fname,'_',name2,'.pdf',sep='')
  
  ggsave(fout,p1)
  
}



samToFrac=function(anoSample,fraction){
  inter=intersect(names(anoSample),names(fraction))
  tab=unique(data.frame('sample'=anoSample[inter],'frac'=fraction[inter]))
  rownames(tab)=tab$sample

  tab2=tab[,2] 
  names(tab2)=rownames(tab)
  
  return(tab2)
}





proportion.box=function(anoCell,anoSample,fraction,appname,cname){
  
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
  
  p=ggplot(tmp4, aes(x=annot, fill=sample, y = pc.of.sample)) + geom_bar(stat='identity', position='fill') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(appname,'.porpotion.png',sep=''),plot=p)
  
  
  tab=samToFrac(anoSample,fraction)
  
  tmp4$dtype=tab[tmp4$sample]
  
  df=tmp4
  p <- ggplot(na.omit(df),aes(x=annot,y=pc.of.sample,dodge=dtype,fill=dtype))+geom_boxplot(notch=FALSE,outlier.shape=NA)  +  geom_point(position = position_jitterdodge(jitter.width=0.1),color=adjustcolor(1,alpha=0.3),aes(pch=dtype))+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5))  +xlab("") +ylab("fraction of total cells")
  
  
  png(file=paste(appname,'.fraction.png',sep=''),width=600,height=400)
  print(p)
  dev.off();
  
  return(df)
  
  
}


