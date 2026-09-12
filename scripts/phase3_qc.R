# Re-run this only if expr_data/group aren't already loaded
source("load_data.R")
# ===== Phase 3: Quality Control =====

# Remove the 3 extra comparison columns, keep only Old and AD samples
count_data <- expr_data[, !(colnames(expr_data) %in% c("Young_vs_Old", "Old_vs_AD", "Young_vs_AD"))]
log_expr <- log2(as.matrix(count_data) + 1)

keep_samples <- group %in% c("Old", "AD")
log_expr_sub <- log_expr[, keep_samples]
group_sub <- factor(group[keep_samples], levels = c("Old", "AD"))
sample_names_sub <- colnames(log_expr_sub)
# Filter out genes with near-zero expression across all samples
# (keeps genes with mean log2 expression above a minimal threshold)
gene_means <- rowMeans(log_expr_sub)
keep_genes <- gene_means > 1  # keeps genes with reasonable expression

log_expr_filtered <- log_expr_sub[keep_genes, ]
dim(log_expr_filtered)  # check how many genes remain
pca <- prcomp(t(log_expr_filtered), scale. = TRUE)
summary(pca)
# Build a data frame for plotting
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = group_sub,
  Sample = sample_names_sub
)

library(ggplot2)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(title = "PCA of AD vs Old Samples",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
  theme_minimal()

ggsave("pca_plot.png", plot = p, width = 8, height = 6, dpi = 300)
print(p)
# Quick sensitivity check: how many of your significant DEGs still 
# show up if we exclude the outlier sample?
library(limma)
deg_limma <- read.csv("limma_significant_DEGs.csv", row.names = 1)
sum(rownames(deg_check) %in% rownames(deg_limma))
# Sample correlation heatmap
library(pheatmap)

sample_cor <- cor(log_expr_filtered, method = "pearson")

annotation_col <- data.frame(Group = group_sub)
rownames(annotation_col) <- colnames(log_expr_filtered)

png("sample_correlation_heatmap.png", width = 8, height = 8, units = "in", res = 300)
pheatmap(sample_cor, 
         annotation_col = annotation_col,
         main = "Sample Correlation Heatmap",
         fontsize_row = 8, fontsize_col = 8)
dev.off()
getwd()
file.exists("sample_correlation_heatmap.png")
list.files(pattern = "correlation")
exists("log_expr_filtered")
file.info("sample_correlation_heatmap.png")$size
file.remove("sample_correlation_heatmap.png")

png("sample_correlation_heatmap.png", width = 8, height = 8, units = "in", res = 300)
pheatmap(sample_cor, 
         annotation_col = annotation_col,
         main = "Sample Correlation Heatmap",
         fontsize_row = 8, fontsize_col = 8)
dev.off()
file.info("sample_correlation_heatmap.png")$size
sessionInfo()
