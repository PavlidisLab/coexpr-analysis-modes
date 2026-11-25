M<-readRDS("/tmp/claude-30003/-home-paul-XCellXGene/4c0a934e-110e-4375-9d9c-cb86509a7125/scratchpad/bulk_stages.rds")
a<-M$rosmap; b<-M$velmeshev; sh<-intersect(rownames(a$qn),rownames(b$qn))
cnt<-data.frame(item=c("rosmap_genes","velmeshev_genes","shared_genes","rosmap_only","velmeshev_only"),
 n=c(nrow(a$qn),nrow(b$qn),length(sh),nrow(a$qn)-length(sh),nrow(b$qn)-length(sh)))
write.table(cnt,"gene_counts.tsv",sep="\t",quote=FALSE,row.names=FALSE)
rows<-list(); add<-function(stage,reading,r,p=NA,n=length(sh))rows[[length(rows)+1]]<<-data.frame(stage=stage,reading=reading,value=r,n=n)
for(st in c("raw","cpm","log2cpm","qn")){
 A<-a[[st]][sh,]; B<-b[[st]][sh,]
 for(tr in if(st%in%c("raw","cpm")) c("as_is","log2p1_of_gene_mean","mean_of_log2p1") else "as_is"){
  if(tr=="as_is"){ma<-rowMeans(A);mb<-rowMeans(B)} else if(tr=="log2p1_of_gene_mean"){ma<-log2(rowMeans(A)+1);mb<-log2(rowMeans(B)+1)} else {ma<-rowMeans(log2(A+1));mb<-rowMeans(log2(B+1))}
  s<-paste0(st,"|",tr)
  add(s,"pearson_across_genes_of_mean_expr",cor(ma,mb)); add(s,"spearman_across_genes_of_mean_expr",cor(ma,mb,method="spearman"))
  add(s,"pearson_across_genes_of_median_expr",cor(apply(if(tr=="as_is")A else log2(A+1),1,median),apply(if(tr=="as_is")B else log2(B+1),1,median)))
 }}
# sample-level readings on log2cpm and qn
for(st in c("log2cpm","qn")){A<-a[[st]][sh,];B<-b[[st]][sh,]
 C<-cor(A,B); add(st,"mean_of_all_636x41_sample_pairs_pearson_across_genes",mean(C))
 add(st,"mean_over_rosmap_samples_of_r_vs_velmeshev_mean_profile",mean(cor(A,rowMeans(B))))
 add(st,"mean_over_velmeshev_samples_of_r_vs_rosmap_mean_profile",mean(cor(B,rowMeans(A))))
 Ca<-cor(A);Cb<-cor(B);add(st,"mean_within_rosmap_sample_pair_r",mean(Ca[upper.tri(Ca)]));add(st,"mean_within_velmeshev_sample_pair_r",mean(Cb[upper.tri(Cb)]))
 Cs<-cor(A,B,method="spearman");add(st,"mean_of_all_sample_pairs_spearman",mean(Cs))
 # gene-wise mean-centred? (correlation of per-gene expression across 'matched' samples is not defined: different samples)
 }
r<-do.call(rbind,rows); r$stage<-as.character(r$stage); write.table(r,"results.tsv",sep="\t",quote=FALSE,row.names=FALSE)
print(r,digits=4,row.names=FALSE); print(cnt)
