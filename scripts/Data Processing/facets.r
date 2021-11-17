library(facets)
set.seed(1234)
args = commandArgs(trailingOnly=TRUE)
sample = args[1]
rcmat = readSnpMatrix(paste(sample,'.csv',sep=''))
xx = preProcSample(rcmat, ndepth=15)
oo = procSample(xx, min.nhet=5)
fit = emcncf(oo)

print("ploidy")
print(fit$ploidy)

print("purity")
print(fit$purity)

rdatafile = paste(sample,'.RData',sep='')
save(xx, oo, fit, file=rdatafile)

fitfile = paste(sample,'.cncf.txt',sep='')
write.table(fit$cncf, fitfile, sep="\t", row.names=F, col.names=T, quote=F)
pdf(paste(sample,".facets.pdf",sep=''))
plotSample(x=oo, emfit=fit, sname=paste(sample,"; ploidy=",round(fit$ploidy,2),"; purity=",round(fit$purity,2),sep=''))
dev.off()
