# Bulk gene-count and expression-correlation check (R1/Methods claim: 21,518 / 25,326 genes, 21,239 shared, r=0.9)
# Rebuilds ROSMAP & Velmeshev bulk from Eric's source files exactly as 041_MAIN_01_exprmats.R, restricted to the saved
# baseline genes/samples (== cleaned set), and checks equality with the saved bk_exprmats_*.rds. Stages: raw(as given), CPM, log2(CPM+1), QN(=saved).
suppressMessages(library(data.table)); options(stringsAsFactors=FALSE)
out<-"/space/scratch/ericchu/r_cache/041_CH4_FINAL/output/"; dd<-"/space/scratch/ericchu/r_cache/041_CH4_FINAL/data/"
qn<-function(m){r<-apply(m,2,rank,ties.method="min");s<-apply(m,2,sort);mu<-apply(s,1,mean)
 matrix(mu[r],nrow=nrow(m),dimnames=dimnames(m))}   # check against saved below
sv<-lapply(c(rosmap="bk_exprmats_rosmap.rds",velmeshev="bk_exprmats_velmeshev.rds"),function(f)readRDS(paste0(out,f)))
# ROSMAP (FPKM file as used by the pipeline)
x<-fread(paste0(dd,"rosmap/ROSMAP_RNAseq_FPKM_gene.tsv"),data.table=FALSE); rownames(x)<-sub("\\..*$","",x[,2]); x<-as.matrix(x[,3:ncol(x)])
cn<-sub("_[^_]*$","",colnames(x)); cn<-regmatches(colnames(x),regexpr(".+_.+(?=_)",colnames(x),perl=TRUE)); colnames(x)<-cn
raw_ros<-x; cat("ROSMAP file genes",nrow(x),"samples",ncol(x),"\n")
v<-fread(paste0(dd,"velmeshev/ASD_BULK.txt"),data.table=FALSE); rownames(v)<-v[,1]; v<-as.matrix(v[,-1]); storage.mode(v)<-"numeric"
cat("Velmeshev file genes",nrow(v),"samples",ncol(v),"\n"); raw_vel<-v
res<-list(rosmap=raw_ros,velmeshev=raw_vel); M<-list()
for(n in names(res)){g<-rownames(sv[[n]]$exprmat); s<-colnames(sv[[n]]$exprmat)
 ct<-res[[n]][g,s,drop=FALSE]; cpm<-t(t(ct)/colSums(ct))*1e6; lg<-log2(cpm+1); q<-qn(lg)
 cat(n,"max|QN recompute - saved| =",max(abs(q-sv[[n]]$exprmat)),"\n")
 M[[n]]<-list(raw=ct,cpm=cpm,log2cpm=lg,qn=sv[[n]]$exprmat)}
saveRDS(M,"/tmp/claude-30003/-home-paul-XCellXGene/4c0a934e-110e-4375-9d9c-cb86509a7125/scratchpad/bulk_stages.rds")
