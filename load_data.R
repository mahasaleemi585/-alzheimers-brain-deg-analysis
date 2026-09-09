expr_data <- read.delim("GSE104704_RNA-Seq_Table.txt.gz", header = TRUE, row.names = 1)
dim(expr_data)
head(expr_data[, 1:5])
sample_cols <- expr_data[, 1:30]
sample_names <- colnames(sample_cols)
group <- c(
  "Young","Young","Young","Young","Young","Young","Young","Young",
  "Old","Old","Old","Old","Old","Old","Old","Old","Old","Old",
  "AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD"
)
data.frame(sample_names, group)
group_lookup <- c(
  "2-12A"="Young", "3-17T"="Young", "4-13A"="Young", "5-18T"="Young",
  "6-14A"="Young", "7-19T"="Young", "8-15A"="Young", "9-16A"="Young",
  "10-8A"="Old", "11-10T"="Old", "12-6A"="Old", "13-11T"="Old",
  "14-7A"="Old", "15-13T"="Old", "16-14T"="Old", "17-9A"="Old",
  "18-10A"="Old", "19-11A"="Old",
  "20-1T"="AD", "21-1A"="AD", "22-2T"="AD", "23-2A"="AD",
  "24-3T"="AD", "25-5T"="AD", "26-3A"="AD", "27-5A"="AD",
  "28-8T"="AD", "29-6T"="AD", "30-9T"="AD", "31-7T"="AD"
)

sample_ids <- sub("^X", "", sample_names)
sample_ids <- sub("^([0-9]+\\.[0-9]+[A-Z])\\..*$", "\\1", sample_ids)
sample_ids <- sub("\\.", "-", sample_ids)

group <- group_lookup[sample_ids]

data.frame(sample_names, sample_ids, group)
exists("sample_names")
group_lookup <- c(
  "2-12A"="Young", "3-17T"="Young", "4-13A"="Young", "5-18T"="Young",
  "6-14A"="Young", "7-19T"="Young", "8-15A"="Young", "9-16A"="Young",
  "10-8A"="Old", "11-10T"="Old", "12-6A"="Old", "13-11T"="Old",
  "14-7A"="Old", "15-13T"="Old", "16-14T"="Old", "17-9A"="Old",
  "18-10A"="Old", "19-11A"="Old",
  "20-1T"="AD", "21-1A"="AD", "22-2T"="AD", "23-2A"="AD",
  "24-3T"="AD", "25-5T"="AD", "26-3A"="AD", "27-5A"="AD",
  "28-8T"="AD", "29-6T"="AD", "30-9T"="AD", "31-7T"="AD"
)

sample_ids <- sub("^X", "", sample_names)
sample_ids <- sub("^([0-9]+\\.[0-9]+[A-Z])\\..*$", "\\1", sample_ids)
sample_ids <- sub("\\.", "-", sample_ids)

group <- group_lookup[sample_ids]

data.frame(sample_names, sample_ids, group)
comparison_cols <- c("Young_vs_Old", "Old_vs_AD", "Young_vs_AD")
count_data <- expr_data[, !(colnames(expr_data) %in% comparison_cols)]

ncol(count_data)
identical(colnames(count_data), sample_names)
comparison_cols <- c("Young_vs_Old", "Old_vs_AD", "Young_vs_AD")
count_data <- expr_data[, !(colnames(expr_data) %in% comparison_cols)]

ncol(count_data)
identical(colnames(count_data), sample_names)
col_data <- data.frame(row.names = sample_names, group = group)
head(col_data)
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
library(DESeq2)
head(count_data[, 1:3])
library(DESeq2)
count_matrix <- round(as.matrix(count_data))
storage.mode(count_matrix) <- "integer"

head(count_matrix[, 1:3])
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ group)
dds
dds$group <- relevel(dds$group, ref = "Old")
levels(dds$group)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

nrow(dds)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "AD", "Old"))
summary(res)
res_ordered <- res[order(res$padj), ]
head(res_ordered)
write.csv(as.data.frame(res_ordered), file = "DESeq2_AD_vs_Old_results.csv")
getwd()
sum(res$padj < 0.1 & abs(res$log2FoldChange) > 1, na.rm = TRUE)
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "group")
deg <- res_ordered[which(res_ordered$padj < 0.05 & abs(res_ordered$log2FoldChange) > 1), ]
up <- deg[deg$log2FoldChange > 0, ]
down <- deg[deg$log2FoldChange < 0, ]

nrow(up)
nrow(down)
write.csv(as.data.frame(deg), "significant_DEGs.csv")
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
pca_data[order(pca_data$PC2), ]
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "padj",
                title = "AD vs Old",
                pCutoff = 0.05,
                FCcutoff = 1)
v <- EnhancedVolcano(res,
                     lab = rownames(res),
                     x = "log2FoldChange",
                     y = "padj",
                     title = "AD vs Old",
                     pCutoff = 0.05,
                     FCcutoff = 1)

ggsave("volcano_plot.png", plot = v, width = 10, height = 8, dpi = 300)
library(pheatmap)
top_genes <- rownames(res_ordered)[1:30]

p <- pheatmap(assay(vsd)[top_genes, ], scale = "row")

ggsave("heatmap.png", plot = p, width = 8, height = 10, dpi = 300)
install.packages("pheatmap")
png("heatmap.png", width = 8, height = 10, units = "in", res = 300)
pheatmap(assay(vsd)[top_genes, ], scale = "row")
dev.off()
png("heatmap.png", width = 8, height = 10, units = "in", res = 300)
pheatmap(assay(vsd)[top_genes, ], scale = "row")
dev.off()
cat(rownames(deg), sep = "\n")
cat(rownames(res_ordered)[1:30], sep = "\n")
