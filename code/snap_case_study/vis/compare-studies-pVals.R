
### Option A: Use nlme package (already in base R) ===================================================
library(nlme)
library(emmeans)
library(tidyverse)

# 1. Descriptive statistics
d = d1_fixed %>% filter(gene_list != "ribo")

# 2. Fit mixed-effects model
model <- lme(spearman_cor ~ source, 
             random = ~ 1 | cell_type/gene_list,
             data = d,
             method = "REML")

# 3. Check results
summary(model)
anova(model)

# 4. Pairwise comparisons
emm <- emmeans(model, pairwise ~ source, adjust = "tukey")

print("\nEstimated Marginal Means:")
print(summary(emm$emmeans))

print("\nPairwise Comparisons:")
print(summary(emm$contrasts))

# 5. Get confidence intervals
confint(emm$contrasts)












### Option B: Use aggregated data with simple linear model ===================================================
# Aggregate to avoid pseudo-replication
d1_aggregated <- d1_fixed %>%
  group_by(cell_type, gene_list, source) %>%
  summarise(
    mean_pearson = mean(pearson_cor, na.rm = TRUE),
    n_comparisons = n(),
    .groups = 'drop'
  )

# Fit linear model
model <- lm(mean_pearson ~ source + cell_type + gene_list, 
            data = d1_aggregated)

summary(model)

# Pairwise comparisons
library(emmeans)
emm <- emmeans(model, pairwise ~ source, adjust = "tukey")
summary(emm$contrasts)

### Option C: Non-parametric approach (no special packages needed) ===================================================
# Aggregate first
d1_aggregated <- d1_fixed %>%
  group_by(cell_type, gene_list, source) %>%
  summarise(
    mean_pearson = mean(spearman_cor, na.rm = TRUE),
    .groups = 'drop'
  )

# 1. Kruskal-Wallis test (non-parametric alternative to ANOVA)
kruskal.test(mean_pearson ~ source, data = d1_aggregated)

# 2. Post-hoc pairwise Wilcoxon tests
pairwise.wilcox.test(d1_aggregated$mean_pearson, 
                     d1_aggregated$source, 
                     p.adjust.method = "bonferroni")






# Step 1: Kruskal-Wallis test on FULL data
kruskal_result <- kruskal.test(spearman_cor ~ source, data = d1_fixed)

# Step 2: ONLY if Kruskal-Wallis is significant, proceed to post-hoc
if(kruskal_result$p.value < 0.05) {
  
  # Option A: Pairwise Wilcoxon on SAME full data
  wilcox_result <- pairwise.wilcox.test(d1_fixed$spearman_cor, 
                                        d1_fixed$source, 
                                        p.adjust.method = "bonferroni")
  print(wilcox_result)
}






