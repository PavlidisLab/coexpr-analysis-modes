
### this is to plot SNAP cor distribution plots [per study] 


dir <- "/cosmos/data/project-data/NW-rewiring/coExpr/ling-PEER4"
cells = c("Ast", "Inh", "Exc")
cells = c("Ast", "Exc")

studies <- c(#"Lau-2020-ALL", "Lim-2022-ALL", 
             "ling-CTL"#, 
             #"Nagy-2020-ALL",
             #"Pineda-2021-ALL", "ROSMAP-Mathys-2019-ALL", "Velmeshev-2019-ALL"
             )


#### process dataset ============================= ============================= =============================
library(data.table)
read_and_tag <- function(study, source) {
  message("Reading ", study, " — ", source)
  file <- file.path(dir, sprintf("%s_PEER4_gene_correlations_with_background_%s.csv", study, source))
  fread(file) %>%
    filter(cell_type %in% cells) %>%
    mutate(source = source)
}


combined <- map(studies, function(study) {
  map_dfr(c("xCell", "xCellxSubject", "xSubject"), ~ read_and_tag(study, .x))
}) %>%
  set_names(studies)


#saveRDS(combined, "/cosmos/data/project-data/NW-rewiring/coExpr/ling-PEER4/with_background_combined_corrected.rds")


combined1 <- lapply(combined, function(df) {
  setDT(df)
  
  # Compute background stats per source × cell_type
  bg_stats <- df[gene_list == "background_random1000",
                 .(bg_mean = mean(correlation, na.rm = TRUE),
                   bg_sd   = sd(correlation, na.rm = TRUE)),
                 by = .(source, cell_type)]
  
  # Merge stats back into main dataframe
  df <- merge(df, bg_stats, by = c("source", "cell_type"), all.x = TRUE)
  
  # Compute corrected correlations
  df[, corrected_cor_mean := correlation - bg_mean]
  df[, corrected_cor_z := (correlation - bg_mean) / bg_sd]
  
  # Rescale corrected correlations to [0, 1] within each source × cell_type
  #df[, corrected_cor_z_scaled :=
  #     if (all(is.na(corrected_cor_z))) NA_real_
  #   else (corrected_cor_z - min(corrected_cor_z, na.rm = TRUE)) /
  #     (max(corrected_cor_z, na.rm = TRUE) - min(corrected_cor_z, na.rm = TRUE)),
  #   by = .(source, cell_type)]
  
  # Optionally clean up
  df[, c("bg_mean", "bg_sd") := NULL]
  
  return(df)
})














#### plotting - Main ============================= ============================= =============================
plot_correlation_freqpoly_combined <- function(data,
                                               correlation_col = "correlation",
                                               celltype_col = "cell_type",
                                               gene_list_col = "gene_list",
                                               source_col = "source",
                                               binwidth = 0.05,
                                               ncol = 3,
                                               title = "Gene-gene correlation frequency across datasets") {
  # Check required columns
  required_cols <- c(correlation_col, celltype_col, gene_list_col, source_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }
  
  p <- ggplot(data, aes_string(x = correlation_col, color = celltype_col)) +
    geom_freqpoly(binwidth = binwidth, linewidth = 1) +
    facet_grid(as.formula(paste(gene_list_col, "~", source_col)), scales = "free_y") + #scales = "free_y", "fixed"
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      x = "Rank-Transformed Correlation",
      y = "Number of gene pairs",
      color = "Cell Type"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA)
    )
  
  if (!is.null(cell_type_palette)) {
    p <- p +
      scale_color_manual(values = cell_type_palette)
  }
  
  return(p)
}


plot_correlation_freqpoly_combined(combined$`ling-CTL`, correlation_col = "correlation")
plot_correlation_freqpoly_combined(combined1$`ling-CTL` %>% filter(!(gene_list == "background_random1000")), correlation_col = "corrected_cor_z", binwidth = 0.05)




all_study_plots <- purrr::map(
  names(combined1),
  ~ {
    df <- combined1[[.x]] %>% 
      filter(gene_list != "background_random1000")
    
    plot_correlation_freqpoly_combined(
      df,
      correlation_col = "corrected_cor_z",
      binwidth = 0.05
    ) +
      ggtitle(paste("Study:", .x))
  }
) %>% set_names(names(combined1))

# 1200 by 1000
all_study_plots[["Lau-2020-ALL"]]
all_study_plots[["Lim-2022-ALL"]]
all_study_plots[["ling-CTL"]]
all_study_plots[["Nagy-2020-ALL"]]
all_study_plots[["Pineda-2021-ALL"]]
all_study_plots[["ROSMAP-Mathys-2019-ALL"]]
all_study_plots[["Velmeshev-2019-ALL"]]






#### plotting - ZOOM in top abs ============================= ============================= =============================
genes = read_csv("/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/Ling/lings_PEER4_gene_loadings_for_nairuz.csv") %>% head(1000)


snap = combined$`ling-CTL` %>% filter(gene_list == "top_abs") %>% select(c("gene1", "gene2")) 
all_genes = unique(c(snap$gene1, snap$gene2))
snap_n = data.frame(gene = all_genes) %>% left_join(genes, by = "gene")


genes_ast = genes %>% filter(cell_type == "astrocyte")
genes_exc = genes %>% filter(cell_type == "glutamatergic")
genes_inh = genes %>% filter(cell_type == "gabaergic")

ast_pos <- genes_ast %>% filter(sign_PEER_4 > 0) %>% nrow()
ast_neg <- genes_ast %>% filter(sign_PEER_4 < 0) %>% nrow()
exc_pos <- genes_exc %>% filter(sign_PEER_4 > 0) %>% nrow()
exc_neg <- genes_exc %>% filter(sign_PEER_4 < 0) %>% nrow()
inh_pos <- genes_inh %>% filter(sign_PEER_4 > 0) %>% nrow()
inh_neg <- genes_inh %>% filter(sign_PEER_4 < 0) %>% nrow()

make_and_plot_top_sets <- function(
    data_list,
    study_name,
    cell_types_to_plot = NULL,  # e.g., c("Ast", "Exc") or c("Ast","Inh")
    show_pairs = TRUE
) {
  
  make_top_sets <- function(celltype, df_study) {
    genes_df <- get(paste0("genes_", celltype))
    df <- df_study %>% dplyr::filter(gene_list == "top_abs")
    
    list(
      pos = df %>% dplyr::filter(
        gene1 %in% genes_df$gene[genes_df$sign_PEER_4 > 0] &
          gene2 %in% genes_df$gene[genes_df$sign_PEER_4 > 0]
      ),
      neg = df %>% dplyr::filter(
        gene1 %in% genes_df$gene[genes_df$sign_PEER_4 < 0] &
          gene2 %in% genes_df$gene[genes_df$sign_PEER_4 < 0]
      )
    )
  }
  
  df_study <- data_list[[study_name]]
  
  top_ast <- make_top_sets("ast", df_study)
  top_exc <- make_top_sets("exc", df_study)
  # top_inh <- make_top_sets("inh", df_study)
  
  combined_top <- dplyr::bind_rows(
    top_ast$pos %>% dplyr::mutate(gene_list = "ast_pos"),
    top_ast$neg %>% dplyr::mutate(gene_list = "ast_neg"),
    top_exc$pos %>% dplyr::mutate(gene_list = "exc_pos"),
    top_exc$neg %>% dplyr::mutate(gene_list = "exc_neg"),
    # top_inh$pos %>% dplyr::mutate(gene_list = "inh_pos"),
    # top_inh$neg %>% dplyr::mutate(gene_list = "inh_neg"),
    df_study %>% dplyr::filter(gene_list == "ribo")
  )
  
  combined_top <- combined_top %>% filter(
    (cell_type == "Ast" & gene_list %in% c("ast_pos", "ast_neg")) |  
    (cell_type == "Exc" & gene_list %in% c("exc_pos", "exc_neg")) | 
    (cell_type %in% c("Ast", "Exc") & gene_list == "ribo")
  )
  
  
  # ---- NEW: choose which cell types (lines) to show ----
  if (is.null(cell_types_to_plot)) {
    cell_types_to_plot <- unique(combined_top$cell_type)
  }
  
  combined_top <- combined_top %>%
    dplyr::filter(cell_type %in% cell_types_to_plot)
  
  # ---- NEW: facet labels with n_genes (and optionally n_pairs) ----
  facet_counts <- combined_top %>%
    dplyr::mutate(
      pair_id = paste(pmin(gene1, gene2), pmax(gene1, gene2), sep = "__")
    ) %>%
    dplyr::group_by(source, gene_list) %>%  # assumes facets are by source x gene_list
    dplyr::summarise(
      n_genes = dplyr::n_distinct(c(gene1, gene2)),
      n_pairs = dplyr::n_distinct(pair_id),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      facet_label = if (show_pairs) {
        paste0(gene_list, " (n=", n_genes, ", pairs=", n_pairs, ")")
      } else {
        paste0(gene_list, " (n=", n_genes, ")")
      }
    )
  
  combined_top_labeled <- combined_top %>%
    dplyr::left_join(facet_counts, by = c("source", "gene_list")) %>%
    dplyr::mutate(gene_list = facet_label)
  
  plot_correlation_freqpoly_combined(
    combined_top_labeled,
    correlation_col = "corrected_cor_z",
    binwidth = 0.05
  ) + ggplot2::ggtitle(study_name)
}

# softer blue "#0072B2"
# 1200 by 1000
cell_type_palette <- c(
  Ast = "#E31A1C",  # red
  Exc = "blue",  # green
  Inh = "#4DAF4A"   # blue
)


make_and_plot_top_sets(combined1, "ling-CTL", cell_types_to_plot = c("Ast", "Exc"), show_pairs = FALSE)

make_and_plot_top_sets(combined1, "Lau-2020-ALL")
make_and_plot_top_sets(combined1, "Lim-2022-ALL")
make_and_plot_top_sets(combined1, "ling-CTL", show_pairs = FALSE)
make_and_plot_top_sets(combined1, "Nagy-2020-ALL")
make_and_plot_top_sets(combined1, "Pineda-2021-ALL")
make_and_plot_top_sets(combined1, "ROSMAP-Mathys-2019-ALL")
make_and_plot_top_sets(combined1, "Velmeshev-2019-ALL")












### same as above but OG code =====================================================
genes = read_csv("/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/Ling/lings_PEER4_gene_loadings_for_nairuz.csv") %>% head(1000)

genes_ast = genes %>% filter(cell_type == "astrocyte")
genes_exc = genes %>% filter(cell_type == "glutamatergic")
genes_inh = genes %>% filter(cell_type == "gabaergic")


make_top_sets <- function(celltype) {
  genes_df <- get(paste0("genes_", celltype))
  df <- combined1$`ling-CTL` %>% filter(gene_list == "top_abs")
  list(
    pos = df %>% filter(gene1 %in% genes_df$gene[genes_df$sign_PEER_4 > 0] &
                          gene2 %in% genes_df$gene[genes_df$sign_PEER_4 > 0]),
    neg = df %>% filter(gene1 %in% genes_df$gene[genes_df$sign_PEER_4 < 0] &
                          gene2 %in% genes_df$gene[genes_df$sign_PEER_4 < 0])
  )
}

top_ast <- make_top_sets("ast")
top_exc <- make_top_sets("exc")
top_inh <- make_top_sets("inh")

combined_top <- bind_rows(
  top_ast$pos %>% mutate(gene_list = "ast_pos"),
  top_ast$neg %>% mutate(gene_list = "ast_neg"),
  top_exc$pos %>% mutate(gene_list = "exc_pos"),
  top_exc$neg %>% mutate(gene_list = "exc_neg"),
  top_inh$pos %>% mutate(gene_list = "inh_pos"),
  top_inh$neg %>% mutate(gene_list = "inh_neg"), 
  combined1$`ling-CTL` %>% filter(gene_list == "ribo") 
)

plot_correlation_freqpoly_combined(combined_top, correlation_col = "corrected_cor_z", binwidth = 0.05)
plot_correlation_freqpoly_combined(combined_top %>% filter(gene_list %in% c("ast_pos", "ast_neg")), correlation_col = "corrected_cor_z", binwidth = 0.05)
plot_correlation_freqpoly_combined(combined_top, correlation_col = "correlation", binwidth = 0.05)



#### ZOOMED IN comparing genes at xCell ============================= ============================= =============================

library(dplyr)
library(tidyr)
library(purrr)

analyze_celltype <- function(data, celltype, gene_list_name = "ast_pos", cutoff = 1) {
  
  # Subsets --------------------------------------------------------
  pos_pos <- data %>%
    filter(
      gene_list == gene_list_name,
      source == "xCell",
      cell_type == celltype,
      corrected_cor_z > cutoff
    )
  
  pos_neg <- data %>%
    filter(
      gene_list == gene_list_name,
      source == "xCell",
      cell_type == celltype,
      corrected_cor_z < -cutoff
    )
  
  # Frequency tables ----------------------------------------------
  freq_tbl <- full_join(
    as.data.frame(table(pos_pos$gene1)) %>% rename(gene = Var1, freq_pos = Freq),
    as.data.frame(table(pos_neg$gene1)) %>% rename(gene = Var1, freq_neg = Freq),
    by = "gene"
  ) %>%
    mutate(
      freq_pos = replace_na(freq_pos, 0),
      freq_neg = replace_na(freq_neg, 0),
      diff_freq = freq_pos - freq_neg
    )
  
  # Unique / shared genes -----------------------------------------
  genes_pos <- unique(c(pos_pos$gene1, pos_pos$gene2))
  genes_neg <- unique(c(pos_neg$gene1, pos_neg$gene2))
  
  unique_pos  <- setdiff(genes_pos, genes_neg)
  unique_neg  <- setdiff(genes_neg, genes_pos)
  shared      <- intersect(genes_pos, genes_neg)
  
  # Gene counts ----------------------------------------------------
  count_pos <- as.data.frame(table(c(pos_pos$gene1, pos_pos$gene2))) %>%
    rename(gene = Var1, count_pos = Freq)
  
  count_neg <- as.data.frame(table(c(pos_neg$gene1, pos_neg$gene2))) %>%
    rename(gene = Var1, count_neg = Freq)
  
  gene_counts <- full_join(count_pos, count_neg, by = "gene") %>%
    replace_na(list(count_pos = 0, count_neg = 0)) %>%
    mutate(diff_count = count_pos - count_neg)
  
  # Output all components in a single list -------------------------
  list(
    pos_pos = pos_pos,
    pos_neg = pos_neg,
    freq_table = freq_tbl,
    unique_pos = unique_pos,
    unique_neg = unique_neg,
    shared = shared,
    gene_counts = gene_counts
  )
}

### gene partner extraction 
get_gene_partners <- function(data, gene) {
  filtered <- data %>%
    filter(gene1 == gene | gene2 == gene)
  
  unique(c(
    filtered$gene1[filtered$gene2 == gene],
    filtered$gene2[filtered$gene1 == gene]
  ))
}

ast_results <- analyze_celltype(combined_top, "Ast")
exc_results <- analyze_celltype(combined_top, "Exc")
ast_results$freq_table
exc_results$unique_neg
exc_results$gene_counts





