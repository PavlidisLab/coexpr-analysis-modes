

### this is to measure reproducability 
library(data.table)
library(dplyr)

#### functions ------------------------------------------------------------------
read_and_tag <- function(study, source) {
  message("Reading ", study, " — ", source)
  file <- file.path(dir, sprintf("%s_PEER4_gene_correlations_with_background_%s.csv", study, source))
  fread(file) %>%
    filter(cell_type %in% cells) %>%
    mutate(source = source)
}
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

compare_studies <- function(df) {
  
  # Format wide study-edge matrix
  edge_matrix <- df %>%
    mutate(edge = paste(gene1, gene2, sep = "_")) %>%
    select(edge, study, corrected_cor_z) %>%
    pivot_wider(names_from = study, values_from = corrected_cor_z)
  
  studies <- colnames(edge_matrix)[-1]
  
  # Must have ≥2 studies
  if (length(studies) < 2) return(NULL)
  
  study_pairs <- combn(studies, 2, simplify = FALSE)
  
  map_dfr(study_pairs, function(p) {
    s1 <- p[1]
    s2 <- p[2]
    
    vals <- edge_matrix %>% select(all_of(c(s1, s2)))
    vals <- vals %>% filter(!is.na(.data[[s1]]) & !is.na(.data[[s2]]))
    
    pearson_cor <- suppressWarnings(cor(vals[[s1]], vals[[s2]], method = "pearson"))
    spearman_cor <- suppressWarnings(cor(vals[[s1]], vals[[s2]], method = "spearman"))
    
    pos_s1 <- which(vals[[s1]] > 0)
    pos_s2 <- which(vals[[s2]] > 0)
    jaccard_pos <- length(intersect(pos_s1, pos_s2)) / length(union(pos_s1, pos_s2))
    
    sign_preservation <- mean(sign(vals[[s1]]) == sign(vals[[s2]]))
    
    # --- Top-k overlap ---
    k = 50
    ord_s1 <- order(abs(vals[[s1]]), decreasing = TRUE)
    ord_s2 <- order(abs(vals[[s2]]), decreasing = TRUE)
    
    topk_s1 <- ord_s1[1:min(k, length(ord_s1))]
    topk_s2 <- ord_s2[1:min(k, length(ord_s2))]
    
    topk_overlap <- length(intersect(topk_s1, topk_s2)) / min(length(topk_s1), length(topk_s2))
    
    
    tibble(
      study1 = s1,
      study2 = s2,
      pearson_cor,
      spearman_cor,
      jaccard_pos,
      sign_preservation, 
      topk_overlap
    )
  })
}
run_one_combination <- function(df, source, cell_type, gene_list) {
  
  sub <- df %>%
    filter(
      source == !!source,
      cell_type == !!cell_type,
      gene_list == !!gene_list
    )
  
  res <- compare_studies(sub)
  if (is.null(res)) return(NULL)
  
  res %>%
    mutate(
      source = source,
      cell_type = cell_type,
      gene_list = gene_list,
      .before = study1
    )
}
run_all_comparisons <- function(df) {
  
  sources <- sort(unique(df$source))
  cell_types <- sort(unique(df$cell_type))
  gene_lists <- sort(unique(df$gene_list))
  
  combos <- expand_grid(
    source = sources,
    cell_type = cell_types,
    gene_list = gene_lists
  )
  
  map_dfr(seq_len(nrow(combos)), function(i) {
    run_one_combination(
      df,
      combos$source[i],
      combos$cell_type[i],
      combos$gene_list[i]
    )
  })
}

add_preserved_geneset <- function(df, base_geneset, new_geneset) {
  
  base <- df %>%
    filter(gene_list == base_geneset, source %in% sources_all) %>%
    mutate(
      g_low  = pmin(gene1, gene2),
      g_high = pmax(gene1, gene2),
      edge_id = paste(g_low, g_high, sep = "__")
    )
  
  preserved_ids <- base %>%
    distinct(cell_type, source, edge_id) %>%
    count(cell_type, edge_id, name = "n_sources") %>%
    filter(n_sources == length(sources_all)) %>%
    select(cell_type, edge_id)
  
  preserved_rows <- base %>%
    semi_join(preserved_ids, by = c("cell_type", "edge_id")) %>%
    select(-g_low, -g_high, -edge_id) %>%
    mutate(gene_list = new_geneset)
  
  bind_rows(df, preserved_rows)
}

normalize_gene_list <- function(x) {
  recode(
    x,
    astro_pos = "ast_pos",
    neur_pos  = "exc_pos",
    top_abs   = "top_abs", # keep if used
    ribo      = "ribo",
    background_random1000 = "background_random1000",
    .default = x
  )
} #Helper to standardize gene_list names 
project_to_template <- function(df_study, study_name) {
  
  df_keyed <- df_study %>% filter(gene_list %in% c("top_abs", "ribo"))
    #mutate(
      #gene_list = normalize_gene_list(gene_list),
      #g_low  = pmin(gene1, gene2),
      #g_high = pmax(gene1, gene2),
      #edge_id = paste(g_low, g_high, sep="__")
    #) %>%
    # if you can have duplicates per key, collapse deterministically:
    #group_by(source, cell_type, gene_list, edge_id) %>%
    #summarise(
      #correlation = dplyr::first(correlation),
      #corrected_cor_mean = dplyr::first(corrected_cor_mean),
      #corrected_cor_z = dplyr::first(corrected_cor_z),
      #.groups = "drop"
    #)
  
  out <- template %>%
    left_join(df_keyed, by = c("source", "cell_type", "gene1", "gene2")) %>%
    mutate(study = study_name) %>%
    select(source, cell_type, gene1, gene2, correlation, gene_list.x , study,
           corrected_cor_mean, corrected_cor_z)
  out$gene_list  = out$gene_list.x
  
  out
} # 3) map one study df -> "same dimensions as template"

#### gene sets 1.0 ------------------------------------------------------------------
genes = read_csv("/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/Ling/lings_PEER4_gene_loadings_for_nairuz.csv") %>% head(1000)
genes_ast = genes %>% filter(cell_type == "astrocyte")
genes_exc = genes %>% filter(cell_type == "glutamatergic")



#### all data studies ------------------------------------------------------------------
dir <- "/cosmos/data/project-data/NW-rewiring/coExpr/ling-PEER4"
cells = c("Ast", "Exc")

studies <- c("Lau-2020-ALL", "Lim-2022-ALL", 
  "ling-CTL", 
  "Nagy-2020-ALL",
  "Pineda-2021-ALL", "ROSMAP-Mathys-2019-ALL", "Velmeshev-2019-ALL", "Ramos-2022-ALL"
)

combined <- map(studies, function(study) {
  map_dfr(c("xCell", "xCellxSubject", "xSubject"), ~ read_and_tag(study, .x))
}) %>%
  set_names(studies)

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
rm(combined)


#### gene sets 2.0 ------------------------------------------------------------------
df_study <- combined1$`ling-CTL`

top_ast <- make_top_sets("ast", df_study)
top_exc <- make_top_sets("exc", df_study)

combined_top <- dplyr::bind_rows(
  top_ast$pos %>% dplyr::mutate(gene_list = "ast_pos"),
  top_ast$neg %>% dplyr::mutate(gene_list = "ast_neg"),
  top_exc$pos %>% dplyr::mutate(gene_list = "exc_pos"),
  top_exc$neg %>% dplyr::mutate(gene_list = "exc_neg"),
  df_study %>% dplyr::filter(gene_list == "ribo")
)

threshold <- 1
threshold_sets <- combined_top %>%
  filter(gene_list %in% c("ast_pos", "exc_pos"),
         corrected_cor_z >= threshold) %>%
  mutate(gene_list = case_when(
    gene_list == "ast_pos" ~ "threshold_ast_pos",
    gene_list == "exc_pos" ~ "threshold_exc_pos",
    .default = gene_list
  ))

combined_gene_sets <- bind_rows(combined_top, threshold_sets) %>%   filter(!(gene_list %in% c("ast_neg", "exc_neg")))


sources_all <- c("xCell", "xSubject", "xCellxSubject")
combined_gene_sets2 <- combined_gene_sets %>%
  add_preserved_geneset("threshold_ast_pos", "threshold_ast_pos_preserved_all") %>%
  add_preserved_geneset("threshold_exc_pos", "threshold_exc_pos_preserved_all")
rm(combined_gene_sets) 
rm(combined_top) 
rm(df_study) 



#### filter studies ------------------------------------------------------------------
template <- combined_gene_sets2 %>%
  mutate(
    g_low  = pmin(gene1, gene2),
    g_high = pmax(gene1, gene2),
    edge_id = paste(g_low, g_high, sep="__")
  ) %>%
  select(source, cell_type, gene_list, edge_id, gene1, gene2)  # keep gene1/gene2 for display
#rm(combined_gene_sets2)

combined1_projected <- lapply(seq_along(names(combined1)), function(i) {
  st <- names(combined1)[i]
  cat(sprintf("[%d/%d] Processing study: %s\n", i, length(combined1), st))
  project_to_template(combined1[[st]], st)
})
names(combined1_projected) <- names(combined1)
stopifnot(all(vapply(combined1_projected, nrow, integer(1)) == nrow(combined_gene_sets2))) # Sanity checks



#### compare studies ------------------------------------------------------------------
df <- dplyr::bind_rows(combined1_projected, .id = "study") %>% select(-"gene_list.x") %>% unique()
results_studies <- run_all_comparisons(df)

gene_lists_summary = df %>% filter(study == "ling-CTL")
table(gene_lists_summary$gene_list, gene_lists_summary$source, gene_lists_summary$cell_type)

#### visualization ------------------------------------------------------------------
d1_fixed <- results_studies %>%
  mutate(
    swap = (study2 == "ling-CTL" & study1 != "ling-CTL"),
    study1_new = if_else(swap, study2, study1),
    study2_new = if_else(swap, study1, study2)
  ) %>%
  transmute(
    source, cell_type, gene_list,
    study1 = study1_new,
    study2 = study2_new,
    pearson_cor, spearman_cor, jaccard_pos, sign_preservation, topk_overlap
  ) %>%
  filter(study1 == "ling-CTL")

gene_list_order <- c(
  "ribo",
  "ast_pos",
  "exc_pos",
  "threshold_ast_pos",
  "threshold_exc_pos",
  "threshold_ast_pos_preserved_all",
  "threshold_exc_pos_preserved_all" 
)

# Convert gene_list to factor with your custom order
d1_fixed <- d1_fixed %>%
  mutate(gene_list = factor(gene_list, levels = gene_list_order))

d1_fixed <- d1_fixed %>% filter(gene_list %in% c("ribo", "ast_pos", "exc_pos", "threshold_ast_pos_preserved_all", "threshold_exc_pos_preserved_all" ))
# 1800 by 800 
ggplot(d1_fixed, aes(study1, study2, fill = spearman_cor)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  facet_grid(cell_type ~ gene_list + source) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Spearman", title = "Study–Study Reproducibility")


# 1. ast_pos all genes reproducibility: ast; exc  
# 2. exc_pos all genes reproducibility: ast; exc  
# 3. Thresholded high in xCell: ast; exc  
# 4. Thresholded high in xSubject and xCellxSubject: ast; exc   
# 5. Thresholded high in all: ast; exc  


