# export Eric's saved baseline bulk matrices (post log2+QN) to tsv gene means + rownames
d<-"/space/scratch/ericchu/r_cache/041_CH4_FINAL/output/"
for(n in c("rosmap","velmeshev")){x<-readRDS(paste0(d,"bk_exprmats_",n,".rds"))$exprmat
 o<-data.frame(gene=rownames(x),mean=rowMeans(x),median=apply(x,1,median),sd=apply(x,1,sd),n_samples=ncol(x))
 write.table(o,file.path("/home/paul/xcell-xbulk-revision/A_bulk_check",paste0("saved_",n,"_genestats.tsv")),sep="\t",quote=FALSE,row.names=FALSE)}
