



test1 = data.table::fread("/cosmos/data/project-data/NW-rewiring/coExpr/ling-PEER4/PEER4-top1-overlaps/SNAP_vs_top1_overlapping_edges_ALL_MODES.csv")
test2 = data.table::fread("/cosmos/data/project-data/NW-rewiring/coExpr/ling-PEER4/PEER4-top1-overlaps/SNAP_vs_top1_overlap_summary_ALL_MODES.csv")



ling = test2 %>% filter(study == "ling-CTL", 
                        cell_type %in% c("Ast", "Exc"))


ling = test1 %>% filter(study == "ling-CTL", 
                        cell_type %in% c("Ast", "Exc"))

ling[, .(n_unique_genes = data.table::uniqueN(c(gene1, gene2))), 
     by = .(mode, cell_type)]


# Peak-driving edges vs Top1%  =============================================================
### get the data from compare-studies.R script 

ling_data <- combined1$`ling-CTL` %>% filter(
            source == "xSubject", 
            cell_type %in% c("Ast", "Exc"), 
            gene_list == "top_abs" , 
            corrected_cor_z >= 1 
)

top1 <- data.table::fread("/cosmos/data/project-data/NW-rewiring/coExpr/ling-PEER4/ling-CTL_PEER4_top1_gene_correlations_with_backgroundxSubject.csv")
top1 <- top1 %>% filter(cell_type %in% c("Ast", "Exc"), gene_list == "top1") 




# Function to create edge identifiers (handles gene order)
create_edge_id <- function(gene1, gene2) {
  pmin(gene1, gene2) %>% paste(pmax(gene1, gene2), sep = "_")
}

# 1. Add edge IDs to both datasets
ling_data <- ling_data %>%
  mutate(edge_id = create_edge_id(gene1, gene2))

top1 <- top1 %>%
  mutate(edge_id = create_edge_id(gene1, gene2))

# 2. Calculate overlap for each cell type
results <- lapply(c("Ast", "Exc"), function(ct) {
  
  # Filter peak-driving edges for this cell type
  peak_edges <- ling_data %>%
    filter(cell_type == ct, source == "xSubject") %>%
    pull(edge_id) %>%
    unique()
  
  # Filter top1% edges for this cell type
  top1_edges <- top1 %>%
    filter(cell_type == ct) %>%
    pull(edge_id) %>%
    unique()
  
  # Calculate overlap
  intersection <- intersect(peak_edges, top1_edges)
  union <- union(peak_edges, top1_edges)
  
  # Calculate metrics
  n_peak <- length(peak_edges)
  n_top1 <- length(top1_edges)
  n_intersection <- length(intersection)
  n_union <- length(union)
  
  jaccard <- n_intersection / n_union
  pct_peak_in_top1 <- (n_intersection / n_peak) * 100
  
  data.frame(
    cell_type = ct,
    n_peak_driving = n_peak,
    n_top1 = n_top1,
    n_intersection = n_intersection,
    n_union = n_union,
    jaccard_index = jaccard,
    pct_peak_in_top1 = pct_peak_in_top1
  )
})

# Combine results
overlap_results <- bind_rows(results)
print(overlap_results)





# (preserved) Peak-driving edges vs Top1%  =============================================================
top1 # from above 

gene_lists_summaryF = gene_lists_summary %>% filter(
  (cell_type == "Ast" & gene_list == "threshold_ast_pos_preserved_all") |
    (cell_type == "Exc" & gene_list == "threshold_exc_pos_preserved_all")
) %>% select(-c(source, correlation, corrected_cor_z, corrected_cor_mean)) %>% unique()
  # from compare-studies.R script 

# 1. Add edge IDs to both datasets
gene_lists_summaryF <- gene_lists_summaryF %>%
  mutate(edge_id = create_edge_id(gene1, gene2))

top1 <- top1 %>%
  mutate(edge_id = create_edge_id(gene1, gene2))

# 2. Calculate overlap and Jaccard for Ast
ast_threshold <- gene_lists_summaryF %>%
  filter(cell_type == "Ast", gene_list == "threshold_ast_pos_preserved_all") %>%
  pull(edge_id)

ast_top1 <- top1 %>%
  filter(cell_type == "Ast") %>%
  pull(edge_id)

ast_top <- (intersect(ast_threshold, ast_top1))
ast_overlap <- length(intersect(ast_threshold, ast_top1))
ast_union <- length(union(ast_threshold, ast_top1))
ast_jaccard <- ast_overlap / ast_union

# 3. Calculate overlap and Jaccard for Exc
exc_threshold <- gene_lists_summaryF %>%
  filter(cell_type == "Exc", gene_list == "threshold_exc_pos_preserved_all") %>%
  pull(edge_id)

exc_top1 <- top1 %>%
  filter(cell_type == "Exc") %>%
  pull(edge_id)

exc_top <- (intersect(exc_threshold, exc_top1))
exc_overlap <- length(intersect(exc_threshold, exc_top1))
exc_union <- length(union(exc_threshold, exc_top1))
exc_jaccard <- exc_overlap / exc_union

# Summary results
results <- data.frame(
  cell_type = c("Ast", "Exc"),
  threshold_edges = c(length(ast_threshold), length(exc_threshold)),
  top1_edges = c(length(ast_top1), length(exc_top1)),
  overlap = c(ast_overlap, exc_overlap),
  jaccard_index = c(ast_jaccard, exc_jaccard)
)

print(results)



# (preserved) Peak-driving edges cross study comparison summary  =============================================================
#from comapre-studies.R script 
d1_fixedF <- d1_fixed %>% filter(
  (cell_type == "Ast" & gene_list == "threshold_ast_pos_preserved_all") |
    (cell_type == "Exc" & gene_list == "threshold_exc_pos_preserved_all"))
  

# Calculate summary statistics across sources
ast_summary <- d1_fixedF %>% 
  filter(cell_type == "Ast") %>%
  filter(gene_list == "threshold_ast_pos_preserved_all") %>% 
  group_by(source) %>%
  summarise(
    n_comparisons = n(),
    mean_spearman = mean(spearman_cor),
    sd_spearman = sd(spearman_cor),
    median_spearman = median(spearman_cor),
    range_spearman = paste0("[", round(min(spearman_cor), 2), "-", round(max(spearman_cor), 2), "]"),
    mean_sign_preservation = mean(sign_preservation),
    sd_sign_preservation = sd(sign_preservation)
  )

exc_summary <- d1_fixedF %>% 
  filter(cell_type == "Exc") %>%
  filter(gene_list == "threshold_exc_pos_preserved_all") %>% 
  group_by(source) %>%
  summarise(
    n_comparisons = n(),
    mean_spearman = mean(spearman_cor),
    sd_spearman = sd(spearman_cor),
    median_spearman = median(spearman_cor),
    range_spearman = paste0("[", round(min(spearman_cor), 2), "-", round(max(spearman_cor), 2), "]"),
    mean_sign_preservation = mean(sign_preservation),
    sd_sign_preservation = sd(sign_preservation)
  )

# (preserved) Peak-driving edges overlap with top1 xstudy reproducability summary  =============================================================
df <- df %>%
  mutate(edge_id = paste(pmin(gene1, gene2), pmax(gene1, gene2), sep = "_"))

# Filter the data
df_filtered <- df %>%
  filter(
    (gene_list == "threshold_ast_pos_preserved_all" & edge_id %in% ast_top) |
      (gene_list == "threshold_exc_pos_preserved_all" & edge_id %in% exc_top)
  )

results_studies2 <- run_all_comparisons(df_filtered)

d11_fixed <- results_studies2 %>%
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

ast_summary1 <- d11_fixed %>% 
  filter(cell_type == "Ast") %>%
  filter(gene_list == "threshold_ast_pos_preserved_all") %>% 
  group_by(source) %>%
  summarise(
    n_comparisons = n(),
    mean_spearman = mean(spearman_cor),
    sd_spearman = sd(spearman_cor),
    median_spearman = median(spearman_cor),
    range_spearman = paste0("[", round(min(spearman_cor), 2), "-", round(max(spearman_cor), 2), "]"),
    mean_sign_preservation = mean(sign_preservation),
    sd_sign_preservation = sd(sign_preservation)
  )

exc_summary1 <- d11_fixed %>% 
  filter(cell_type == "Exc") %>%
  filter(gene_list == "threshold_exc_pos_preserved_all") %>% 
  group_by(source) %>%
  summarise(
    n_comparisons = n(),
    mean_spearman = mean(spearman_cor),
    sd_spearman = sd(spearman_cor),
    median_spearman = median(spearman_cor),
    range_spearman = paste0("[", round(min(spearman_cor), 2), "-", round(max(spearman_cor), 2), "]"),
    mean_sign_preservation = mean(sign_preservation),
    sd_sign_preservation = sd(sign_preservation)
  )






