msigdb.h <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c2.cp.kegg.v6.2.symbols.gmt')
msigdb.c2 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c2.all.v6.2.symbols.gmt')
msigdb.c5 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c5.all.v6.2.symbols.gmt')
msigdb.c6 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c6.all.v6.2.symbols.gmt')
msigdb.c7 <- read.gmt('/d0-mendel/home//meisl/lib/msigdb/c7.all.v6.2.symbols.gmt')


msigdb.h=msigdb.c5

matrices2 <- do.call(c,mapply(function(m, name) {
  rep(name,length(m))
},
msigdb.h,
names(msigdb.h)))


matrices3 <- do.call(c,mapply(function(m, name) {
  m
},
msigdb.h,
names(msigdb.h)))


d1=data.frame(matrices2,matrices3)
d2=unique(data.frame(matrices2,matrices2))
USER_DATA <- build_Anno(d1,d2)

r=enricher_internal(gene = markers, pvalueCutoff = 0.1, 
                    pAdjustMethod = "BH", universe = universe, minGSSize = 10, 
                    maxGSSize = 500, qvalueCutoff = 0.1, USER_DATA = USER_DATA)

