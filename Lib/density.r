

anoCell=cluster2[cname]
fraction=samp2[cname]
unique(fraction)

plotlis=list()
myeloid_emb=con$embedding
for( iterm in unique(fraction)){
  
  #iterm='Begin'
  data2=con$embedding[names(fraction[fraction==iterm]),]
  colnames(data2)=c('x','y')
  data2=as.data.frame(data2)
  dim(data2)
  
  # raster polygon  xlim(xlims) + ylim(ylims)+

  p=con$plotGraph(groups=anoCell[names(fraction[fraction==iterm])],plot.na=FALSE,alpha=0.1,show.legend=FALSE,mark.groups = TRUE)
  
  
  
  gg=ggplot(data2, aes(x=x, y=y) ) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)  +xlab('')+ylab('')+
    #scale_fill_distiller(palette=2, direction=0.1,expand = c(0, 0))+     # +ggtitle(iterm)
    
    scale_fill_viridis(option='B',alpha = 1,direction=1) +  #ggtitle(iterm)+
    geom_point(aes(x=x,y=y), col='#FCFDBFFF',size=0.00001,alpha=0.2)+
    #scale_x_continuous(expand = c(0, 0)) +
    #scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(
      legend.position='none',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      #legend.position='none',
      plot.margin = margin(0,0,0,0,"cm")) 
  
  gg=gg+theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))
  
  gg <- gg + theme(axis.title.x=element_blank(),
                   # axis.text.x=element_text(size=7,margin = margin(0, unit = "cm")),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_line(),
                   axis.title.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   # axis.text.y=element_text(size=7,margin = margin(0, unit = "cm")),
                   axis.text.y=element_blank(),
                   axis.ticks.length = unit(0, "cm"))
  
  
  ggsave(paste(iterm,'.pdf',sep=''),gg,width = 2.5,height=2.5)
  
  
  
  
  
  
  
  
  
  plotlis[[iterm]]=p
}


b=  cowplot::plot_grid(plotlist=plotlis, ncol=3, nrow=2)

b

ggsave(paste(appname,'.density.png',sep=''),b,width = 12,height=8)

