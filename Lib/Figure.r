#

library(ggpubr)

data("ToothGrowth")
head(ToothGrowth)


compare_means(len ~ supp, data = ToothGrowth)


p <- ggboxplot(ToothGrowth, x = "supp", y = "len",
               color = "supp", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()
# Change method
p + stat_compare_means(method = "t.test"
                       
                       
                       
 
 # Visualize: Specify the comparisons you want
 my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
 ggboxplot(ToothGrowth, x = "dose", y = "len",
           color = "dose", palette = "jco")+ 
   stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
   stat_compare_means(label.y = 50)     # Add global p-value                      
# 








diffPlot4 <- eventReactive(input$submit4, {
  gene4 <- sepGene(input$gene4)[1]
  updateTextInput(session, inputId = "gene4", value = paste(gene4, collapse = " "))
  validate(
    need(gene4 %in% geneAliasUpper$Symbol, "")      
  )
  epat <- cbind(meta, gene=get(load(paste0("data/expr/", geneIndex[gene4], ".Rdata")))[,gene4])
  epat <- expandSub(epat, skcm = F)
  d <- filter(epat, tissue %in% c("01","03","06","11"), !is.na(gene)) %>% 
    mutate(tissue=as.factor(as.character(tissue)))
  levels(d$tissue) <- c("Tumor","Tumor", "Metastasis","Normal")
  d <- mutate(d, caLabel=paste0(cancer, ".", tissue)) %>% 
    filter(grepl("Tumor", caLabel) | grepl("Normal", caLabel) | grepl("SKCM\\.Metastasis", caLabel) & !grepl("SKCM\\.Normal", caLabel))
  allLevels <- c("ACC.Tumor", "BLCA.Tumor", "BLCA.Normal", "BRCA.Tumor", "BRCA.Normal", "BRCA-Basal.Tumor", "BRCA-Her2.Tumor", "BRCA-Luminal.Tumor", "CESC.Tumor", "CHOL.Tumor", "CHOL.Normal", "COAD.Tumor", "COAD.Normal", "DLBC.Tumor", "ESCA.Tumor", "ESCA.Normal", "GBM.Tumor", "HNSC.Tumor", "HNSC.Normal", "HNSC-HPVpos.Tumor", "HNSC-HPVneg.Tumor", "KICH.Tumor", "KICH.Normal", "KIRC.Tumor", "KIRC.Normal", "KIRP.Tumor", "KIRP.Normal", "LAML.Tumor", "LGG.Tumor", "LIHC.Tumor", "LIHC.Normal", "LUAD.Tumor", "LUAD.Normal", "LUSC.Tumor", "LUSC.Normal", "MESO.Tumor", "OV.Tumor", "PAAD.Tumor", "PCPG.Tumor", "PRAD.Tumor", "PRAD.Normal", "READ.Tumor", "READ.Normal", "SARC.Tumor", "SKCM.Tumor", "SKCM.Metastasis", "STAD.Tumor", "STAD.Normal", "TGCT.Tumor", "THCA.Tumor", "THCA.Normal", "THYM.Tumor", "UCEC.Tumor", "UCEC.Normal", "UCS.Tumor", "UVM.Tumor")
  d$caLabel <- factor(d$caLabel, levels = allLevels)
  tmp <- table(d$caLabel) < 2
  levels(d$caLabel)[allLevels %in% names(tmp[tmp])] <- NA
  d <- filter(d, !is.na(caLabel))
  d$caColor <- d$caLabel
  tmp <- levels(d$caColor)
  tmp[grep("Tumor", tmp)] <- 1
  tmp[grep("Normal", tmp)] <- 2
  tmp[grep("Metastasis", tmp)] <- 3
  levels(d$caColor) <- tmp
  rect <- data.frame(x=c(grep("\\.[NM]", levels(d$caLabel)),grep("HPVneg", levels(d$caLabel)))) %>% 
    mutate(xmin=x-1.45, xmax=x+0.45)
  caPairs <- data.frame(x=rect$x-0.5, Tumor=levels(d$caLabel)[rect$x-1], Nor=levels(d$caLabel)[rect$x])
  caPairs$p <- apply(caPairs, 1, function(x) wilcox.test(gene~caLabel, data = subset(d, caLabel %in% x))$p.value)
  caPairs$labels <- cut(caPairs$p, breaks = c(0, 0.001, 0.01, 0.05, 0.1), labels = c("***","**","*","·"))
  p <- ggplot(d, aes(x=caLabel, y = gene))+
    labs(x="", y=paste(gene4,"Expression Level (log2 RSEM)"))+
    geom_rect(data=rect, aes(NULL,NULL,xmin=xmin, xmax=xmax), ymin=-Inf, ymax=Inf, fill="Light Gray")+
    geom_jitter(aes(color=caColor),position = position_jitter(width = .3), alpha=0.3)+
    geom_boxplot(aes(color=caColor),outlier.colour = NA, alpha=0.5)+
    scale_color_manual(values=c("red","blue","purple"))+
    geom_text(data = caPairs, aes(x=x, label=labels), y=Inf, vjust=2, size=7)+
    theme_bw(base_size = 17)+
    theme(legend.position="none", axis.text.x = element_text(angle=90, hjust=1, vjust=.5,color="black"))
  txt <- caPairs[,-c(1,5)]
  list(p=p, txt=txt)
})

#  color size 

ggplot(data,aes(x=x, y=y, color=z,size=z1)) +   #z1
	geom_point(pch=15)+
#	geom_jitter()+
  	theme(axis.title.x = element_text(face = "bold", size = 0.5),axis.title.y = element_text(face = "bold", size = 0.8),plot.background = element_blank(),axis.text.x=element_text(size=6,face = "bold",angle=90),axis.text.y=element_text(size=5.5,face= "bold"),
 	 panel.background = element_blank(),axis.ticks.y=element_blank(),axis.ticks.length=unit(0,'cm'))+
  	scale_x_continuous(breaks=1:Ncol,labels=colnames(dat))+
  	scale_y_continuous(breaks=1:Nrow,labels=rownames(dat))+
  	scale_size(name="RP",range=c(1,3))+
#  	labs(title = name[4])+
  	theme(plot.margin = unit(c(0,0.5,0.5,0.3),"cm"))+
  	scale_colour_gradient2(name='Cor', low = "blue",
    mid = "white", high = "red", midpoint = 0,
    space = "rgb", na.value = "grey50",
    guide = "colourbar")


dev.off()








#  dot plot with label 
library(ggplot2)
suppressPackageStartupMessages(library(ggrepel))

InsectSprays[,2]=as.character(InsectSprays[,2])

p=ggplot(InsectSprays,aes(x=count,y=count,text =paste("name:", spray)))+ theme_bw() 
p=p+geom_point(shape=".",alpha=1/1,size =1)
p=p+geom_jitter(size = 1)
#  p=p+labs(x=g1,y=g2,title=paste(round(cc,4),round(cc2,4),sep='/'))
p=p+labs(x='FC',y='-log(p-value)')

nn=c('t1','t2')
InsectSprays[2,2]='t1'

p=p+geom_vline(xintercept = 0,linetype = "dotted")
p=p+geom_hline(yintercept = 0,linetype = "dotted")

p=p+geom_text_repel(
  data = InsectSprays[InsectSprays$spray %in% nn,],
  aes(label = spray),
  size = 5,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines")
)
p





# density

set.seed(1234)
df <- data.frame(
  sex=factor(rep(c("F", "M"), each=200)),
  weight=round(c(rnorm(200, mean=55, sd=5),
                 rnorm(200, mean=65, sd=5)))
  )
head(df)









##
library(ggplot2)

p <- ggplot(InsectSprays, aes(x=spray, y=count, fill=spray)) +geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16, position=position_jitter(0.2))
p=p+theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.text.x=element_blank())
p





p <- ggplot(InsectSprays, aes(x=spray, y=count, fill=spray)) +geom_violin()+ geom_boxplot(width=0.1)
p



#box plot

nn=as.factor(nn)
nn <- ordered(nn, levels = c("Tumor", "Involved", "Noninvolved","Whole"))


dat=data.frame('Faction'=nn,'Percentage'=tmp3['B_lymphoid',])


p <- ggplot(dat, aes(x=Faction, y=Percentage, fill=Faction)) +geom_boxplot(outlier.shape = NA)

p=p+theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.text.x=element_blank())

ggsave('b1.png',p,height=4,width=6)







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



