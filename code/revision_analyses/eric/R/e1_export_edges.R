# Export positive-edge sets (value >= 1-frac, lower triangle) from Eric Chu's saved rank-normalised coexpression matrices
# (/space/scratch/ericchu/r_cache/041_CH4_FINAL/output/{sc,sbj,bk}_coexmats_*.rds). Read-only on the source. Pair id = max(gi,gj)*G + min(gi,gj) in a global C-sorted gene list.
suppressMessages(library(parallel))
O <- "/space/scratch/ericchu/r_cache/041_CH4_FINAL/output/"; OUT <- path.expand("~/xcell-xbulk-revision/eric/edges/")
args <- commandArgs(TRUE); which <- args[1]; ncores <- as.integer(args[2])
gm <- readRDS(paste0(O,"genes_metadata.rds")); G0 <- sort(unique(gm$gene_id), method="radix"); G <- length(G0)
writeLines(G0, paste0(OUT,"global_genes.txt")); cat("global genes", G, "\n")
fracs <- c(0.005,0.01,0.05,0.10)
files <- list.files(O, pattern=paste0("^", which, "_coexmats_.*\\.rds$"))
files <- files[!grepl("consensus", files)]
one <- function(f) {
  nm <- sub("\\.rds$","",f); t0 <- Sys.time()
  M <- readRDS(paste0(O,f)); genes <- rownames(M); stopifnot(identical(genes, colnames(M)))
  gi <- match(genes, G0); stopifnot(!anyNA(gi)); n <- length(genes)
  sha <- strsplit(system(paste("sha256sum", shQuote(paste0(O,f))), intern=TRUE), " ")[[1]][1]
  writeBin(as.integer(gi), paste0(OUT, nm, "__genes.i32"), size=4)
  res <- list(file=f, n_genes=n, sha256=sha, n_na=sum(is.na(M)), sorted_c=!is.unsorted(genes), mono_global=!is.unsorted(gi))
  for (fr in fracs) {
    k <- which(M >= 1-fr); r <- (k-1) %% n + 1; cc <- (k-1) %/% n + 1; keep <- r > cc
    a <- gi[r[keep]]; b <- gi[cc[keep]]; id <- pmax(a,b) * G + pmin(a,b)   # fits int32 if G*G < 2^31
    stopifnot(max(id) < 2^31)
    id <- sort(id); writeBin(as.integer(id), paste0(OUT, nm, "__", fr, ".i32"), size=4); res[[paste0("n_pos_", fr)]] <- length(id)
    rm(k, r, cc, keep, a, b, id)
  }
  res$secs <- round(as.numeric(difftime(Sys.time(), t0, units="secs"))); cat(nm, n, res$secs, "\n"); res
}
out <- mclapply(files, one, mc.cores = ncores, mc.preschedule = FALSE)
tab <- do.call(rbind, lapply(out, function(x) as.data.frame(x[sapply(x, length)==1], stringsAsFactors=FALSE)))
write.table(tab, paste0(OUT, "qc_", which, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
