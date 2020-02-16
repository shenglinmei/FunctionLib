
proporgateBy.cluster=function(cluster,ano){
  cname=intersect(names(cluster),names(ano))
  tab=table(ano[cname],cluster[cname])
  tmp3 <- (sweep(tab, 2, colSums(tab), FUN='/'))
  index=apply(tmp3,2,function(x) which.max(x))
  ratio=apply(tmp3,2,function(x) max(x))
  res=data.frame('cluster'=colnames(tmp3),'cells'=rownames(tmp3)[index],'ratio'=ratio)
  
  
  pdf(file='ratio.pdf',width=8,height=8)
  lcl <- cluster[cname];
  x <- jdf(as.factor(lcl),as.factor(ano[cname]))
  x$l2 <- factor(as.character(x$l2),levels=levels(x$l2)[order(unlist(tapply(1:nrow(x),x$l2,function(ii) { w <- x$ov[ii]^5; sum(as.integer(x$l1)[ii] * w)/sum(w) })),decreasing=F)])
  ggplot(x,aes(x=l2,y=l1))+geom_point(aes(size=sqrt(overlap),alpha=jaccard),col='blue')+scale_size(range = c(0, 6))+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1,size=rel(0.9)), axis.text.y=element_text(size=rel(0.9)))+labs(x='',y='clusters')
  dev.off()
  gc()
  
  
  cluster2=Toch(cluster)
  tmp=as.character(res$cluster)
  names(tmp)=res$cells
  
  cluster3=names(tmp)[match(cluster2,tmp)]
  names(cluster3)=names(cluster2)
  
  return(list('res'=res,'ano'=cluster3))
}



