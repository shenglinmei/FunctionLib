

myvenn=function(res,fout){
  library(VennDiagram)
  library(grid)
  
  dnames=names(res)
  nc=length(dnames)
  cols=c('red','blue','yellow','green','grey')
  aaa=venn.diagram(res,fout,category.names =dnames,filename = NULL,fill=cols[1:nc])
  jpeg(fout);
  grid.draw(aaa);
  dev.off();
  
}



