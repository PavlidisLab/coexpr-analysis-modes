# Export raw-count submatrices (group genes x group cells) for the 39 groups from Eric Chu's sc_exprmats_<dataset>.rds ($raw$ctmat restricted to each cell type's retained genes/cells).
# Binary CSC: <key>__i.i32 (row idx 0-based into <key>__genes.i32 order), __p.i32, __x.f64, __genes.i32 (global gene idx, C-sorted global list), __cellsbj.txt (subject id per column, column order).
suppressMessages({library(Matrix); library(parallel)})
O <- "/space/scratch/ericchu/r_cache/041_CH4_FINAL/output/"; OUT <- path.expand("~/xcell-xbulk-revision/eric/expr/")
G0 <- readLines(path.expand("~/xcell-xbulk-revision/eric/edges/global_genes.txt"))
ds_all <- c("lau","lim","nagy","pineda","ramos","rosmap","velmeshev")
have <- sub("^sc_coexmats_","",sub("\\.rds$","",list.files(O, pattern="^sc_coexmats_.*\\.rds$"))); have <- have[!grepl("consensus",have)]
one <- function(ds) {
  x <- readRDS(paste0(O,"sc_exprmats_",ds,".rds")); C0 <- as(x$raw$ctmat, "CsparseMatrix"); res <- list()
  for (ct in names(x$exprmats)) {
    key <- paste0(ds,"_",ct); if (!(key %in% have)) next
    E <- x$exprmats[[ct]]; genes <- rownames(E); cells <- colnames(E); S <- x$samples[[ct]]
    stopifnot(all(cells %in% S$sample)); sbj <- S$sample_bk[match(cells, S$sample)]
    C <- C0[genes, cells]; C <- as(C, "CsparseMatrix"); stopifnot(max(C@x) < 2^31, all(C@x == round(C@x)))
    gi <- match(genes, G0); stopifnot(!anyNA(gi))
    pre <- paste0(OUT, key)
    writeBin(as.integer(C@i), paste0(pre,"__i.i32"), size=4); writeBin(as.integer(C@p), paste0(pre,"__p.i32"), size=4); writeBin(as.double(C@x), paste0(pre,"__x.f64"), size=8)
    writeBin(as.integer(gi), paste0(pre,"__genes.i32"), size=4); writeLines(as.character(sbj), paste0(pre,"__cellsbj.txt"))
    # also CPM check vector: colsums of retained-gene counts (to confirm CPM definition) for first 5 cells
    cs <- Matrix::colSums(C)[1:5]; cpm_chk <- as.numeric(E[1:3, 1:5]); 
    res[[key]] <- list(n_genes=length(genes), n_cells=length(cells), nnz=length(C@x), n_sbj=length(unique(sbj)), cpm_def_ok = isTRUE(all.equal(as.numeric(E[1:3,1:5]), as.numeric(C[1:3,1:5] %*% diag(1e6/cs)), tolerance=1e-8)))
    cat(key, res[[key]]$n_genes, res[[key]]$n_cells, res[[key]]$nnz, "cpm_def_ok", res[[key]]$cpm_def_ok, "\n")
  }
  res
}
out <- mclapply(ds_all, one, mc.cores=7, mc.preschedule=FALSE)
tab <- do.call(rbind, lapply(unlist(out, recursive=FALSE), function(r) as.data.frame(r)))
tab$key <- names(unlist(out, recursive=FALSE)); write.table(tab, paste0(OUT,"export_qc.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
