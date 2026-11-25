
### ODDS ratio comparison for different sources edges 


###. prep data  ------------------------------------------------------------------
genes = read_csv("/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/Ling/lings_PEER4_gene_loadings_for_nairuz.csv") %>% head(1000)

ast_pos_genes <- genes %>% filter(cell_type == "astrocyte") %>% filter(sign_PEER_4 > 0) 
exc_pos_genes <- genes %>% filter(cell_type == "glutamatergic") %>% filter(sign_PEER_4 > 0)
#ast_neg <- genes %>% filter(cell_type == "astrocyte") %>% filter(sign_PEER_4 < 0) 
#exc_neg <- genes %>% filter(cell_type == "glutamatergic") %>% filter(sign_PEER_4 < 0)  

test = combined1$`ling-CTL` %>% filter(gene_list == "top_abs", 
                                       gene1 %in% ast_pos$gene & gene2 %in% ast_pos$gene#, 
                                       #cell_type %in% c("Ast", "Exc")
)




filter_edges_by_genes <- function(
    df,
    int_genes,
    int_gene_list = NULL,
    keep_self_loops = FALSE,
    undirected = TRUE
) {
  stopifnot(is.data.frame(df))
  
  out <- df
  
  if (!is.null(gene_list)) out <- dplyr::filter(out, gene_list == int_gene_list)
  out <- dplyr::filter(out, gene1 %in% int_genes, gene2 %in% int_genes)
  
  # optional cleanup
  if (!keep_self_loops) {
    out <- dplyr::filter(out, gene1 != gene2)
  }
  
  # if you want an undirected edge representation, enforce canonical ordering
  if (undirected) {
    out <- dplyr::mutate(out,
                         gene_low  = pmin(gene1, gene2),
                         gene_high = pmax(gene1, gene2)
    ) |>
      dplyr::select(-gene1, -gene2) |>
      dplyr::rename(gene1 = gene_low, gene2 = gene_high)
  }
  
  out
}


ast_df <- filter_edges_by_genes( df = combined1$`ling-CTL`, int_genes = ast_pos_genes$gene, int_gene_list = "top_abs")
exc_df <- filter_edges_by_genes( df = combined1$`ling-CTL`, int_genes = exc_pos_genes$gene, int_gene_list = "top_abs")

ast_pos_list = list (`ling-CTL` = ast_df) 
exc_pos_list = list (`ling-CTL` = exc_df) 




### functions ------------------------------------------------------------------
library(data.table)

compute_peak_or_by_source <- function(
    combined_list,
    dataset,
    cell_types = NULL,
    gene_list = NULL,
    score_col = "corrected_cor_z",
    threshold = 1.5,
    sources = NULL,
    direction = "pos",          # default for all sources: "pos" | "neg" | "abs"
    direction_by_source = NULL  # named list or named vector: e.g. c(xCell="pos", xSubject="neg")
) {
  # ---- input checks ----
  if (!dataset %in% names(combined_list)) {
    stop("`dataset` not found in `combined_list`. Available: ", paste(names(combined_list), collapse = ", "))
  }
  
  dt <- as.data.table(combined_list[[dataset]])
  
  if (!score_col %in% names(dt)) {
    stop("`score_col` not found in dataset. Available columns include: ", paste(names(dt), collapse = ", "))
  }
  
  # ---- filters ----
  if (!is.null(cell_types)) dt <- dt[cell_type %in% cell_types]
  if (!is.null(gene_list))  dt <- dt[gene_list %in% gene_list]
  if (nrow(dt) == 0) stop("No rows left after filtering (check dataset/cell_types/gene_list).")
  
  # ---- sources ----
  if (is.null(sources)) sources <- unique(dt$source)
  sources <- intersect(sources, unique(dt$source))
  if (length(sources) < 2) stop("Need at least 2 sources after filtering.")
  
  cells <- unique(dt$cell_type)
  
  # ---- edge id (order-invariant) ----
  dt[, edge := fifelse(gene1 < gene2,
                       paste(gene1, gene2, sep = "_"),
                       paste(gene2, gene1, sep = "_"))]
  
  # ---- direction map per source ----
  # normalize direction strings
  norm_dir <- function(x) {
    x <- tolower(x)
    if (!x %in% c("pos", "neg", "abs")) stop("Direction must be one of: 'pos', 'neg', 'abs'.")
    x
  }
  
  dir_map <- setNames(rep(norm_dir(direction), length(sources)), sources)
  
  if (!is.null(direction_by_source)) {
    # allow named vector OR named list
    if (is.list(direction_by_source)) {
      direction_by_source <- unlist(direction_by_source)
    }
    if (is.null(names(direction_by_source))) {
      stop("`direction_by_source` must be named by source (e.g., c(xCell='pos', xSubject='neg')).")
    }
    for (nm in intersect(names(direction_by_source), sources)) {
      dir_map[[nm]] <- norm_dir(direction_by_source[[nm]])
    }
  }
  
  # ---- peak function factory ----
  make_peak_fun <- function(dir, thr) {
    if (dir == "pos") return(function(x) x >= thr)
    if (dir == "neg") return(function(x) x <= -thr)
    # abs
    function(x) abs(x) >= thr
  }
  
  # ---- compute OR for one pair ----
  compute_or_pair <- function(cell, src1, src2) {
    peak1 <- make_peak_fun(dir_map[[src1]], threshold)
    peak2 <- make_peak_fun(dir_map[[src2]], threshold)
    
    d1 <- dt[cell_type == cell & source == src1,
             .(edge, is_peak = peak1(get(score_col)))]
    d2 <- dt[cell_type == cell & source == src2,
             .(edge, is_peak = peak2(get(score_col)))]
    
    merged <- merge(d1, d2, by = "edge", all = TRUE)
    merged[is.na(is_peak.x), is_peak.x := FALSE]
    merged[is.na(is_peak.y), is_peak.y := FALSE]
    
    a <- factor(merged$is_peak.x, levels = c(FALSE, TRUE))
    b <- factor(merged$is_peak.y, levels = c(FALSE, TRUE))
    tab <- table(a, b)
    
    ft <- fisher.test(tab)
    
    data.table(
      dataset   = dataset,
      gene_list = if (is.null(gene_list)) NA_character_ else paste(gene_list, collapse = ","),
      score_col = score_col,
      threshold = threshold,
      
      cell_type = cell,
      source_1  = src1,
      source_2  = src2,
      dir_1     = dir_map[[src1]],
      dir_2     = dir_map[[src2]],
      
      odds_ratio = unname(ft$estimate),
      conf_low   = ft$conf.int[1],
      conf_high  = ft$conf.int[2],
      p_value    = ft$p.value,
      
      n_edges  = nrow(merged),
      n_peak_1 = sum(merged$is_peak.x),
      n_peak_2 = sum(merged$is_peak.y),
      
      tab_FF = tab[1,1],
      tab_FT = tab[1,2],
      tab_TF = tab[2,1],
      tab_TT = tab[2,2]
    )
  }
  
  # ---- run all ----
  results <- rbindlist(
    lapply(cells, function(cell) {
      rbindlist(
        combn(sources, 2, simplify = FALSE, FUN = function(x) {
          compute_or_pair(cell, x[1], x[2])
        }),
        use.names = TRUE, fill = TRUE
      )
    }),
    use.names = TRUE, fill = TRUE
  )
  
  results[]
}








###. RUN ANALYSIS ------------------------------------------------------------------
###. positive peaks 
res_ast<- compute_peak_or_by_source(
  combined_list = ast_pos_list,
  dataset = "ling-CTL",
  cell_types = c("Ast", "Exc"),
  gene_list = "top_abs",
  score_col = "corrected_cor_z",
  threshold = 1,
  direction = "pos"
)

res_exc<- compute_peak_or_by_source(
  combined_list = exc_pos_list,
  dataset = "ling-CTL",
  cell_types = c("Ast", "Exc"),
  gene_list = "top_abs",
  score_col = "corrected_cor_z",
  threshold = 1,
  direction = "pos"
)


res <- compute_peak_or_by_source(
  combined_list = combined1,
  dataset = "ling-CTL",
  cell_types = c("Ast", "Exc"),
  gene_list = "neur_pos",
  score_col = "corrected_cor_z",
  threshold = 1,
  direction = "pos"
)
res

###. positive peaks vs negative 
res_pos_vs_neg <- compute_peak_or_by_source(
  combined_list = combined1,
  dataset = "ling-CTL",
  cell_types = c("Ast", "Exc"),
  gene_list = "astro_pos",
  score_col = "corrected_cor_z",
  threshold = 1.5,
  direction = "pos",
  direction_by_source = c(xCell = "neg")  # override just this one source
)
res_pos_vs_neg

###. absolute peaks 
res_abs <- compute_peak_or_by_source(
  combined_list = combined1,
  dataset = "ling-CTL",
  cell_types = c("Ast", "Exc"),
  gene_list = "astro_pos",
  threshold = 1.5,
  direction = "abs"
)


# 1600 vs 350 
ggplot(res_exc, aes(
  x = odds_ratio,
  y = paste(source_1, source_2, sep = " vs ")
)) +
  geom_point(size = 6) +
  geom_errorbarh(
    aes(xmin = conf_low, xmax = conf_high),
    width = 0.6
  ) +
  scale_x_log10() +
  facet_wrap(~ cell_type, scales = "free_y") +
  theme_classic(base_size = 12) +
  labs(
    x = "Odds ratio (log scale)",
    title = "Enrichment of peak-driving gene-gene edges across analysis modes"
  ) +
  theme(axis.text = element_text(size = 12))










