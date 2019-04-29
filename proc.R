

library(nbHelpers)
library(parallel)
library(pagoda2)


ss1 <- read.csv('primary_samplesheet.csv',stringsAsFactors=F)
prefix <- ss1$path
names(prefix) <- ss1$sample
prefix <- prefix[!prefix == '']
suffix <- '/outs/filtered_gene_bc_matrices/hg19'
fullpaths <- paste0(prefix,suffix)
names(fullpaths) <- names(prefix)

raw.mats <- mclapply(fullpaths, function(p) {
    try({
        read10xMatrix(p)
    })
},mc.cores=32)

## Prefix the cells with the names
named.mats <- mapply(function(m,n) {
    colnames(m) <- paste0(n,'-',colnames(m))
    m
},raw.mats, names(raw.mats))



p2objs <- mclapply(named.mats, function(m) {
    try({
        basicP2proc(m,n.cores=2)
    })
},mc.cores=28)



p2objs <- p2objs[!unlist(lapply(p2objs,is.error))]

ls()
rm(fullpaths,named.mats,prefix,raw.mats,ss1,suffix)
gc()

## How many cells?
unlist(lapply(p2objs, function(p2) {nrow(p2$counts)}))

## Make a conos object and get a joint clustering
library(conos)
con <- Conos$new(p2objs,n.cores=30)
con$buildGraph(verbose=T)
con$findCommunities()

fastSave::save.image.lbzip2('savepoint1.RDataFS',n.cores=32)
## 

p <- con$plotPanel()
p
ggplot2::ggsave('jc_overview.png',p,width=10,height=10)

## JC
jc1 <- as.factor(con$clusters$multi$groups)
head(jc1)

p2wobjs <- lapply(namedNames(p2objs), function(n) {
    p2 <- p2objs[[n]]
    jc.app <- jc1[rownames(p2$counts)]
    additionalMetadata <- factorListToMetadata(list(JC=jc.app))
    wo <- basicP2web(p2,app.title=n,n.cores=4,extraWebMetadata=additionalMetadata)
    wo
})

## Save bins
lapply(names(p2wobjs), function(n) {
    wo <- p2wobjs[[n]]
    wo$serializeToStaticFast(paste0('apps/',n,'.bin'))
})

paste0('http://pklab.med.harvard.edu/nikolas/pagoda2/frontend/current/pagodaURL/index.html?fileURL=http://pklab.med.harvard.edu/nikolas/hirz/pcaProc20180814/',list.files('apps/'))

## save final state
fastSave::save.image.lbzip2('savepoint2.RDataFS',n.cores=28);




## Genereate N only, T only and N+T
fastSave::load.lbzip2('savepoint2.RDataFS',n.cores=28);

## reload raw.mats
library(nbHelpers)
library(parallel)
library(pagoda2)
ss1 <- read.csv('primary_samplesheet.csv',stringsAsFactors=F)
prefix <- ss1$path
names(prefix) <- ss1$sample
prefix <- prefix[!prefix == '']
suffix <- '/outs/filtered_gene_bc_matrices/hg19'
fullpaths <- paste0(prefix,suffix)
names(fullpaths) <- names(prefix)
raw.mats <- mclapply(fullpaths, function(p) {
    try({
        read10xMatrix(p)
    })
},mc.cores=32)
## Prefix the cells with the names
named.mats <- mapply(function(m,n) {
    colnames(m) <- paste0(n,'-',colnames(m))
    m
},raw.mats, names(raw.mats))
rm(raw.mats); gc()

exclude.samples <- c('SCG_PCA3A')

all.sample.names <- names(named.mats)
all.sample.names <- setdiff(all.sample.names, exclude.samples)
all.sample.names

T.sample.names <- all.sample.names[grepl('-T',all.sample.names)]
N.sample.names <- all.sample.names[grepl('-N',all.sample.names)]

## Tumor
T.bigmat <- do.call(cbind,named.mats[T.sample.names])
T.p2 <- basicP2proc(T.bigmat,n.cores=32)
snms <- unlist(strpart(rownames(T.p2$counts),'-[ATGC]{5,}',1))
names(snms) <- rownames(T.p2$counts)
snms <- as.factor(snms)
meta <- list(sample=snms, jc=jc1[names(snms)])
meta2 <- factorListToMetadata(meta)
T.p2w <- basicP2web(T.p2, 'PCA-Tumor', extraWebMetadata= factorListToMetadata(meta))
T.p2w$serializeToStaticFast('pca-tumor.bin')
rm(T.bigmat,T.p2,T.p2w); gc()

## Normal
N.bigmat <- do.call(cbind, named.mats[N.sample.names])
N.p2 <- basicP2proc(N.bigmat,n.cores=32)
snms <- unlist(strpart(rownames(N.p2$counts),'-[ATGC]{5,}',1))
names(snms) <- rownames(N.p2$counts)
snms <- as.factor(snms)
meta <- list(sample=snms, jc=jc1[names(snms)])
meta2 <- factorListToMetadata(meta)
N.p2w <- basicP2web(N.p2, 'PCA-Normal', extraWebMetadata=meta2)
N.p2w$serializeToStaticFast('pca-normal.bin')
rm(N.bigmat, N.p2,N.p2w)
gc()

## All
all.bigmat <- do.call(cbind, named.mats[c(N.sample.names, T.sample.names)])

set.seed(2018);cells.keep <- sample(colnames(all.bigmat),50000)

all.p2 <- basicP2proc(all.bigmat[,cells.keep], n.cores=32)

snms <- unlist(strpart(rownames(all.p2$counts),'-[ATGC]{5,}',1))
names(snms) <- rownames(all.p2$counts)
snms <- as.factor(snms)

sampletype <- ifelse(grepl('-T',as.character(snms)),'Tumor','Normal')
names(sampletype) <- names(snms)
sampletype <- as.factor(sampletype)

individual <- strpart(names(sampletype),'-',1)
names(individual) <- names(sampletype)
individual <- as.factor(individual)


meta <- list(sample=snms, jc=jc1[names(snms)],sampletype=sampletype,individual=individual)
meta2 <- factorListToMetadata(meta)


all.p2w <- basicP2web(all.p2, 'PCA-All-Subset', extraWebMetadata=meta2)
all.p2w$serializeToStaticFast('pca-all.bin')


rm(N.bigmat, N.p2,all.p2w)
gc()


