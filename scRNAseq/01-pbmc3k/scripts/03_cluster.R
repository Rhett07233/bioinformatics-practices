#!/usr/bin/env Rscript
# 步骤 3：归一化 → 高变基因 → 标准化 → PCA → 聚类 → UMAP
#
# 这是单细胞流程的主体。每一步为什么这么做，见 README 第 3 节。
set.seed(1)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# 结果输出根目录：R 的结果放 results/R，Python（Scanpy）的结果放 results/python。
RES_DIR <- Sys.getenv("RES_DIR", file.path("results", "R"))
dir.create(file.path(RES_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

pbmc <- readRDS(file.path(RES_DIR, "pbmc3k_qc.rds"))

# ---- 3.1 归一化 ----
# 每个细胞的测序深度不同，直接比表达量会把"测得深"当成"表达高"。
# LogNormalize：除以该细胞总 counts × 10000，再取 log。
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ---- 3.2 高变基因 ----
# 只用细胞间变异大的基因做下游，可以压低技术噪声、突出生物学信号。
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 高变基因：", paste(top10, collapse = ", "), "\n\n")

ggsave(file.path(RES_DIR, "figures/03_variable_features.png"),
       LabelPoints(VariableFeaturePlot(pbmc), points = top10, repel = TRUE),
       width = 8, height = 6, dpi = 150)

# ---- 3.3 标准化 ----
# 把每个基因的表达转成均值 0、方差 1，否则高表达基因会在 PCA 里占过大权重。
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# ---- 3.4 PCA ----
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

stdev <- Stdev(pbmc, reduction = "pca")
pct <- stdev / sum(stdev) * 100
write.csv(data.frame(PC = seq_along(pct),
                     pct_variance = round(pct, 3),
                     cum_pct = round(cumsum(pct), 3)),
          file.path(RES_DIR, "pca_variance.csv"), row.names = FALSE)

cat("前 10 个 PC 解释的方差（%）：", paste(round(pct[1:10], 2), collapse = ", "), "\n")
cat(sprintf("前 10 个 PC 累计解释 %.1f%%；前 20 个累计 %.1f%%\n\n",
            sum(pct[1:10]), sum(pct[1:20])))

ggsave(file.path(RES_DIR, "figures/04_elbow.png"), ElbowPlot(pbmc, ndims = 30),
       width = 7, height = 5, dpi = 150)
ggsave(file.path(RES_DIR, "figures/05_pca.png"),
       DimPlot(pbmc, reduction = "pca") + NoLegend(),
       width = 7, height = 6, dpi = 150)
ggsave(file.path(RES_DIR, "figures/06_pca_loadings.png"),
       VizDimLoadings(pbmc, dims = 1:2, reduction = "pca"),
       width = 8, height = 7, dpi = 150)

# ---- 3.5 选 PC 数和分辨率 ----
# 这两个是整条流程里最影响结果的参数：
#   NPCS —— 用前多少个主成分做聚类。肘部图（04_elbow.png）会告诉你拐点在哪。
#           官方教程在这套数据上取 10（拐点在 PC9-10）。
#           宁可取大不要取小：只取 5 个 PC 会明显破坏结果。
#   RES  —— 聚类分辨率。越大 cluster 越多。3K 细胞常用 0.4~1.2。
# 换成别的值重跑：NPCS=15 RES=0.8 Rscript scripts/03_cluster.R
NPCS <- as.integer(Sys.getenv("NPCS", "10"))
RES  <- as.numeric(Sys.getenv("RES",  "0.5"))
if (NPCS > length(pct)) {
  stop(sprintf("NPCS=%d 超过了实际算出的 PC 数（%d）", NPCS, length(pct)))
}
cat(sprintf("== 使用 dims = 1:%d，resolution = %.2f ==\n", NPCS, RES))

# ---- 3.6 KNN 图 + Louvain 聚类 ----
# FindNeighbors：在 PCA 空间里按欧氏距离建 KNN 图，再用 Jaccard 相似度重新加权边。
# FindClusters：在这张图上做模块度优化（默认 Louvain），把图切成社区。
pbmc <- FindNeighbors(pbmc, dims = 1:NPCS)
pbmc <- FindClusters(pbmc, resolution = RES)

# ---- 3.7 UMAP ----
# 注意：UMAP 只用于可视化，不要只凭 UMAP 上"看起来分得开"就下生物学结论。
# 它保留的是局部距离，全局距离没有意义——两个 cluster 之间的远近不能当结论用。
pbmc <- RunUMAP(pbmc, dims = 1:NPCS)

ggsave(file.path(RES_DIR, "figures/07_umap_clusters.png"),
       DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 5) + NoLegend(),
       width = 7.5, height = 6, dpi = 150)

tb <- table(Idents(pbmc))
cat(sprintf("\n共 %d 个 cluster：\n", length(tb)))
print(tb)

write.csv(data.frame(cluster = names(tb), n_cells = as.integer(tb)),
          file.path(RES_DIR, "cluster_sizes.csv"), row.names = FALSE)

saveRDS(pbmc, file.path(RES_DIR, "pbmc3k_clustered.rds"))
cat("\n已保存：", file.path(RES_DIR, "pbmc3k_clustered.rds"), "\n")
