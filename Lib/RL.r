

source('/home/meisl/bin/FunctionLib/Lib/pagodaLib.r')
source('/home/meisl/bin/FunctionLib/Lib/conosLib.r')



scon <- readRDS('../../ncdcon.rds')

p2=readRDS('../../allp2.rds')


to.coarse.ann <- function(ca) {
  ca=gsub(' ', '.', ca)
  ca=gsub('/', '_', ca)
  ca=gsub('\\+', '_', ca)
  ca=gsub('-', '', ca)
  
  return(ca)
}



ano=readRDS('../../annotation.new.1103.rds')

unique(ano)

ano=ano[ano!="Th17 IL17-"]
ano=ano[ano!="NKcell"]



#anoT=readRDS('../Tcell.ano1103.rds')
#anoM=readRDS('../myeloid.ano10.rds')
#anoT=anoT[anoT!="Th17 IL17-"]

#ano=ano[ano %in% c(unique(anoT),unique(anoM),'Fibroblast','Bridge')]
#ano[ano=="CD8+ CTL naïve "]="CD8+ CTL naive"
#ano[ano=="CD4 naïve"]="CD4 naive"


anoT=to.coarse.ann(ano)
anoM=to.coarse.ann(ano)






anoT[1:4]
anoM[1:4]

anoA=anoT

de1 <- p2$getDifferentialGenes(groups=anoA,z.threshold = -1000)
names(de1)





anoSample=apply(data.frame(names(anoA)),1,function(x) strsplit(x,'_')[[1]][1])
names(anoSample)=names(anoA)


fractionf=anoSample
fractionf[!is.na(fractionf)]='NB'

bigM=t(p2$misc$rawCounts)
expA=t(p2$counts)

table(names(anoA) %in% colnames(bigM))


dat=lapply(sn(unique(anoA)), function(x) prepare_rawList3(expA,bigM,anoA,x,anoSample,fractionf,ncut=2))






RL=readRDS('/home/meisl/Workplace/Data/CellPhoneDB/CellDB.rds')
dim(RL)

RL=unique(RL[,c('Ligand','Receptor')])

dim(RL)


dim(expA)
dim(bigM)


index1= as.character(RL[,'Ligand']) %in% rownames(bigM)
RL=RL[index1,]

index1= as.character(RL[,'Receptor']) %in% rownames(bigM)
RL=RL[index1,]
dim(RL)






lis=list()
for(i in unique(anoM)){
  for (j in unique(anoT)){

    
  #  for(i in seq(3)){
  #   for (j in seq(3)){    
    
      cell1=i
      cell2=j
      if (cell1!=cell2){
      gs1=RL$Ligand
      gs2=RL$Receptor
      
      rLR=getInfor(cell1,cell2,dat,gs1,gs2)
     # rRL=getInfor(cell2,cell1,dat,gs1,gs2)
      
      lis[[paste(cell1,cell2,sep='_')]]=rLR
     # lis[[paste(cell2,cell1,sep='_')]]=rRL
    }
  }
}

table(anoA,anoSample)


#cell1="Mono1"
#cell2='Th17.IL17'


fl1=unique(anoM)
fl2=unique(anoT)

tmp1=Outputf(lis,fl1,fl2,'Result',cut.ratio.a=0.1,cut.ratio.b=0.1)

tmp1=Outputf2(lis,fl1,fl2,'Result.cut',cut.ratio.a=0.1,cut.ratio.b=0.1)


save.image(file='F1111.RData')




save(lis,fl2,fl1,file='res2.RData')





tmp2=Outputf(lis,fl2,fl1,'Tcell.To.Myeloid',cut.ratio.a=0.01,cut.ratio.b=0.01)
save(lis,fl2,fl1,file='Myeloid.Tcell.RData')








tmp1=Outputf(lis,fl1,fl2,'Myeloid.To.Tumor',cut.ratio.a=0.01,cut.ratio.b=0.01)
tmp2=Outputf(lis,fl2,fl1,'Tumor.To.Myeloid',cut.ratio.a=0.01,cut.ratio.b=0.01)
save(lis,fl2,fl1,file='Myeloid.Tumor.RData')





tmp1=Outputf(lis,fl1,fl2,'Tcell.To.Tumor',cut.ratio.a=0.01,cut.ratio.b=0.01)
tmp2=Outputf(lis,fl2,fl1,'Tumor.To.Tcell',cut.ratio.a=0.01,cut.ratio.b=0.01)
save(lis,fl2,fl1,file='Tcell.Tumor.RData')





