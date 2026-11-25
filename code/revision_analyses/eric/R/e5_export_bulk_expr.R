O <- "/space/scratch/ericchu/r_cache/041_CH4_FINAL/output/"; OUT <- path.expand("~/xcell-xbulk-revision/eric/expr/")
G0 <- readLines(path.expand("~/xcell-xbulk-revision/eric/edges/global_genes.txt"))
for (nm in c("rosmap","velmeshev")) {
  x <- readRDS(paste0(O,"bk_exprmats_",nm,".rds")); cat(nm, names(x), "\n"); M <- x$exprmat; cat(dim(M), "sorted:", !is.unsorted(rownames(M)), "\n")
  gi <- match(rownames(M), G0); stopifnot(!anyNA(gi))
  writeBin(as.integer(gi), paste0(OUT,"bulk_",nm,"__genes.i32"), size=4); writeBin(as.double(M), paste0(OUT,"bulk_",nm,"__expr.f64"), size=8)   # column-major genes x samples
  writeLines(colnames(M), paste0(OUT,"bulk_",nm,"__samples.txt")); write(dim(M), paste0(OUT,"bulk_",nm,"__dim.txt"))
}
