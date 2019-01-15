
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



	dat.mat=readRDS('/home/meisl/Workplace/BMME/a.data/raw.mats.rds')
	index=names(dat.mat) %in% c("BMET11-PB",'BMET11-TumorEpi')

	dat.mat=dat.mat[!index]


	bigM2=mergeDat2(dat.mat)
	rm(dat.mat)




    input_cell='T_helper'
    
    cname2=names(tmp)
	anoCell=allano[['cell']][cname2]
	anoSample=allano[['sample']][cname2]
	anoSampleType=allano[['sampltType']][cname2]

    
    
S_cells <- anoCell[anoCell %in% input_cell]
print('cell')
head(S_cells)
length(S_cells)

nname=names(S_cells)

anoSample=anoSample[nname]

index=grepl('BMET11-',names(anoSample))
anoSample=anoSample[!index]
index=grepl('BMET10-',names(anoSample))
anoSample=anoSample[!index]

index=grepl('BMET1-',names(anoSample))
anoSample=anoSample[!index]



	gname=names(anoSample)

	t.bigM2=bigM2[,gname]


	tab=table(allano[['sample']][gname])

	#tab=rowSums(nt[,i])

	nnames=names(tab[tab>15])

	cname2=colnames(t.bigM2)
	n.anoCell=allano[['cell']][cname2]
	n.anoSample=allano[['sample']][cname2]
	n.anoSampleType=allano[['sampltType']][cname2]

	raw.mats=list()

	p2.lis=list()
	for (n in nnames){
	  cd=t.bigM2[,n.anoSample==n]
	  # p2 <- Pagoda2$new(cd,min.cells.per.gene = 0,trim = 0)
	  # n1=paste('Group',i,'_',n,sep='')
	  #n1=paste('Group','_',n,sep='')
	  n1=n
	  # p2.lis[[n1]]=p2
	  raw.mats[[n1]]=cd
	}


	nn=unlist(lapply(raw.mats, function(x) ncol(x)))


nn



	raw.mats2=lapply(raw.mats,function(x) t(x))


	jcl3.coarse=as.character(n.anoCell)
	names(jcl3.coarse)=names(n.anoCell)

	#jcl3.coarse[jcl3.coarse!='B_lymphoid']='B_lymphoid'

	jcl3.coarse=as.factor(jcl3.coarse)



#raw.mats2=raw.mats2[!grepl('_BMET1-',names(raw.mats2))]
#ggg=unlist(lapply(raw.mats2,function(x) rownames(x)))
#jcl3.coarse=jcl3.coarse[ggg]


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
	        raw.mats3=lapply(raw.mats2,function(x) x[,names(fc.TvsW)])
	        
	        all.percl.TvsW.nocorr <- getPerCellTypeDECorrected_raw_list(conObj=raw.mats3,
	                                                                    groups=jcl3.coarse,
	                                                                    sampleGroups=sampleGroups,
	                                                                    n.cores=32,correction=fc.TvsW,reflevel='control')
	        
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


	saveRDS(allres,'T_help_noBMET1.rds')






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






save.image(file='Thelp.RData')


