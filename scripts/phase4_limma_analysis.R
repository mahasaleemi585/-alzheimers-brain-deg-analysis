# ===== Phase 4: Differential Expression Analysis (limma) =====
# Uses normalized/processed expression data (not raw counts),
# so limma is used instead of DESeq2.
# NOTE: run load_data.R first - this script needs expr_data and group.

library(limma)

# Remove the 3 extra pre-computed comparison columns
count_data <- expr_data[, !(colnames(expr_data) %in% c("Young_vs_Old", "Old_vs_AD", "Young_vs_AD"))]

# Log2-transform (limma expects roughly normal, log-scale data)
log_expr <- log2(as.matrix(count_data) + 1)

# Keep only Old and AD samples
keep_samples <- group %in% c("Old", "AD")
log_expr_sub <- log_expr[, keep_samples]
group_sub <- factor(group[keep_samples], levels = c("Old", "AD"))

table(group_sub)  # should show Old: 10, AD: 12

# Build design matrix and run limma
design <- model.matrix(~ group_sub)
colnames(design) <- c("Intercept", "AD_vs_Old")

fit <- lmFit(log_expr_sub, design)
fit <- eBayes(fit)

res_limma <- topTable(fit, coef = "AD_vs_Old", number = Inf, sort.by = "P")

# Extract significant genes (adj. p < 0.05, |log2FC| > 1)
deg_limma <- res_limma[res_limma$adj.P.Val < 0.05 & abs(res_limma$logFC) > 1, ]
up <- deg_limma[deg_limma$logFC > 0, ]
down <- deg_limma[deg_limma$logFC < 0, ]

nrow(up)    # should be 73
nrow(down)  # should be 106

write.csv(deg_limma, "limma_significant_DEGs.csv")

# ===== Volcano plot =====
library(EnhancedVolcano)
v <- EnhancedVolcano(res_limma,
                     lab = rownames(res_limma),
                     x = "logFC",
                     y = "adj.P.Val",
                     title = "AD vs Old (limma)",
                     pCutoff = 0.05,
                     FCcutoff = 1)
ggsave("volcano_plot_limma.png", plot = v, width = 10, height = 8, dpi = 300)

# ===== Heatmap of top 30 DEGs =====
library(pheatmap)
top_genes <- rownames(res_limma)[1:30]

png("heatmap_limma.png", width = 8, height = 14, units = "in", res = 300)
pheatmap(log_expr_sub[top_genes, ], scale = "row", fontsize_row = 10, cellheight = 15)
dev.off()