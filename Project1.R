library(tidyverse)
library(janitor)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# 1. LOAD DATA
raw_data <- read_tsv(file.choose())

# 2. CLEANING
cleaned_data <- raw_data %>%
  filter(!grepl("contam", `Protein ID`, ignore.case = TRUE)) %>%
  filter(!grepl("rev_", `Protein ID`, ignore.case = TRUE))

# 3. SELECT QUANTITATION
quant_data <- cleaned_data %>%
  dplyr::select(`Protein ID`, Gene, contains("MaxLFQ")) %>%
  clean_names()

# 4. LOG TRANSFORMATION
log_data <- quant_data %>%
  mutate(across(where(is.numeric), ~ log2(. + 1)))

# 5. FILTERING
log_filtered <- log_data %>%
  filter(
    (rowSums(!is.na(dplyr::select(., starts_with("input_proteome_")))) >= 2) |
      (rowSums(!is.na(dplyr::select(., starts_with("streptavidin_pulldown_")))) >= 2)
  )

# 6. IMPUTATION
numeric_cols <- 3:8
min_val <- min(as.matrix(log_filtered[, numeric_cols]), na.rm = TRUE)

data_imputed <- log_filtered
data_imputed[, numeric_cols][is.na(data_imputed[, numeric_cols])] <- min_val * 0.9

# 7. EXPERIMENTAL DESIGN
groups <- c("Input", "Input", "Input", "Pulldown", "Pulldown", "Pulldown")
group_factor <- factor(groups, levels = c("Input", "Pulldown"))
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)

# 8. LINEAR MODEL (LIMMA)
contrast_matrix <- makeContrasts(Pulldown - Input, levels = design)
fit <- lmFit(data_imputed[, numeric_cols], design)
fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrast)

# 9. RESULTS ANNOTATION
results <- topTable(fit_ebayes, number = Inf, sort.by = "P")

results_annotated <- results %>%
  mutate(row_idx = as.numeric(rownames(.))) %>%
  left_join(
    data_imputed %>% 
      mutate(row_idx = row_number()) %>% 
      dplyr::select(row_idx, protein_id, gene), 
    by = "row_idx"
  ) %>%
  dplyr::select(protein_id, gene, everything(), -row_idx) %>%
  mutate(neg_log10_p = -log10(adj.P.Val)) %>%
  mutate(color_category = case_when(
    adj.P.Val < 0.05 & logFC > 1 ~ "Enriched (Pulldown)",
    adj.P.Val < 0.05 & logFC < -1 ~ "Depleted (Input)",
    TRUE ~ "Not Significant"
  ))

# 10. PLOTTING VOLCANO
ggplot(results_annotated, aes(x = logFC, y = neg_log10_p, color = color_category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Pulldown vs Input", subtitle = "Red = Hits")

ggsave("02_volcano_plot_results.png", width = 8, height = 6)

# 11. GO ENRICHMENT
enriched_genes <- results_annotated %>%
  filter(adj.P.Val < 0.05, logFC > 1) %>%
  pull(gene) %>%
  unique()

ego <- enrichGO(gene          = enriched_genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP", 
                pAdjustMethod = "BH")

dotplot(ego, showCategory=15)
ggsave("05_biological_pathways_dotplot.png", width = 10, height = 8)
