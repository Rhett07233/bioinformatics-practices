#!/usr/bin/env Rscript
# 步骤 4：找每个 cluster 的 marker 基因，并用经典 marker 注释细胞类型
#
# 思路：无监督聚类只给出编号（0、1、2…），编号本身没有生物学意义。
# 要知道"这一群是什么细胞"，得看它们高表达哪些基因——
# 这就是 marker 基因（也叫 cluster biomarker）。
set.seed(1)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# 结果输出根目录：R 的结果放 results/R，Python（Scanpy）的结果放 results/python。
RES_DIR <- Sys.getenv("RES_DIR", file.path("results", "R"))
dir.create(file.path(RES_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

pbmc <- readRDS(file.path(RES_DIR, "pbmc3k_clustered.rds"))

# ---- 4.1 每个 cluster vs 其余所有细胞的差异表达 ----
# only.pos      只保留在该 cluster 里上调的基因
# min.pct       至少要在 25% 的细胞里检测到，否则可能是噪声
# logfc.threshold 平均 log2 倍数变化的下限
markers <- FindAllMarkers(pbmc, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(RES_DIR, "cluster_markers.csv"), row.names = FALSE)

cat(sprintf("共 %d 条 marker 记录，覆盖 %d 个 cluster\n\n",
            nrow(markers), length(unique(markers$cluster))))

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()
write.csv(top10, file.path(RES_DIR, "cluster_top10_markers.csv"), row.names = FALSE)

cat("== 每个 cluster 的 top5 marker ==\n")
print(as.data.frame(
  top10 %>% group_by(cluster) %>% slice_head(n = 5) %>% ungroup() %>%
    select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)
))

# ---- 4.2 经典 marker 的表达模式 ----
# 这些是免疫细胞领域公认的标记基因。看它们在哪些 cluster 高表达，
# 就能把无监督的编号对上真实的细胞类型。
canonical <- c("IL7R", "CCR7", "CD14", "LYZ", "S100A4", "MS4A1", "CD8A",
               "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP")

ggsave(file.path(RES_DIR, "figures/09_dotplot_canonical.png"),
       DotPlot(pbmc, features = canonical) + RotatedAxis(),
       width = 10, height = 5.5, dpi = 150)

avg <- AverageExpression(pbmc, features = canonical, assays = "RNA")$RNA
write.csv(round(as.matrix(avg), 3),
          file.path(RES_DIR, "canonical_marker_avg_expression.csv"))

cat("\n== 经典 marker 在各 cluster 的平均表达（用于人工判断细胞类型）==\n")
print(round(as.matrix(avg), 2))

ggsave(file.path(RES_DIR, "figures/08_heatmap_top10.png"),
       DoHeatmap(pbmc, features = top10$gene) + NoLegend(),
       width = 12, height = 10, dpi = 150)

# ---- 4.3 把编号换成细胞类型 ----
# 官方教程在这套数据（resolution=0.5，dims=1:10）上得到的对应关系：
#   0 Naive CD4 T / 1 CD14+ Mono / 2 Memory CD4 T / 3 B / 4 CD8 T /
#   5 FCGR3A+ Mono / 6 NK / 7 DC / 8 Platelet
#
# 但 cluster 编号不是稳定契约：换了分辨率、换了 Seurat 版本都可能变。
# 所以这里先检查 cluster 数是否符合预期，不符就只给提示、不硬套。
new.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
             "FCGR3A+ Mono", "NK", "DC", "Platelet")

if (nlevels(pbmc) == length(new.ids)) {
  names(new.ids) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, new.ids)

  ggsave(file.path(RES_DIR, "figures/10_umap_annotated.png"),
         DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5) +
           NoLegend() + ggtitle("PBMC 3k — cell type annotation"),
         width = 8, height = 6, dpi = 150)

  cat("\n== 注释后的细胞类型分布 ==\n")
  print(table(Idents(pbmc)))

  saveRDS(pbmc, file.path(RES_DIR, "pbmc3k_final.rds"))
  cat("\n已保存：", file.path(RES_DIR, "pbmc3k_final.rds"), "\n")
} else {
  warning(sprintf("cluster 数是 %d，不是教程里的 9 个，跳过自动注释。\n", nlevels(pbmc)),
          "请对照 results/R/canonical_marker_avg_expression.csv 人工判断每个 cluster 的细胞类型。",
          call. = FALSE)
}
