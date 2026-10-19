#!/usr/bin/env Rscript
# 步骤 2：读入 count matrix，计算 QC 指标并过滤低质量细胞
#
# 对应官方教程的 "QC and selecting cells for further analysis"。
#
# 三个常用的 QC 指标：
#   nFeature_RNA —— 每个细胞检测到的基因数（太少 = 空液滴/死细胞；异常多 = 双细胞）
#   nCount_RNA   —— 每个细胞的 UMI 总数（和 nFeature 强相关）
#   percent.mt   —— 线粒体基因的 counts 占比（高 = 细胞膜破损、濒死）
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)   # 47-48 行用 + 拼两个 FeatureScatter，靠 patchwork 提供 +.gg
})

set.seed(1)

# 结果输出根目录：R 的结果放 results/R，Python（Scanpy）的结果放 results/python。
# 两套实现互不覆盖，05_verify 会拿两个目录下的同名 CSV 逐项对比。
RES_DIR <- Sys.getenv("RES_DIR", file.path("results", "R"))
dir.create(file.path(RES_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

DATA_DIR <- Sys.getenv("PBMC3K_DIR", "data/filtered_gene_bc_matrices/hg19")
if (!dir.exists(DATA_DIR)) {
  stop("找不到数据目录：", DATA_DIR, "\n请先运行：bash scripts/01_download.sh")
}

# ---- 2.1 读入并建对象 ----
# min.cells / min.features 是基因和细胞的初筛门槛：
# 一个基因至少在 3 个细胞里出现；一个细胞至少检测到 200 个基因。
counts <- Read10X(data.dir = DATA_DIR)
pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k",
                           min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

cat(sprintf("读入完成：%d 个基因 × %d 个细胞\n", nrow(pbmc), ncol(pbmc)))

# ---- 2.2 看图，再决定阈值 ----
# 不要照抄教程的阈值。先看自己的图，再定。
# 官方教程在这一步用的是：nFeature 200~2500，percent.mt < 5
ggsave(file.path(RES_DIR, "figures/01_qc_violin_before.png"),
       VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3, pt.size = 0),
       width = 12, height = 4.5, dpi = 150)

ggsave(file.path(RES_DIR, "figures/02_qc_scatter.png"),
       FeatureScatter(pbmc, "nCount_RNA", "percent.mt") +
         FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA"),
       width = 12, height = 5, dpi = 150)

# ---- 2.3 过滤 ----
med_before <- sapply(pbmc@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")], median)
n_before <- ncol(pbmc)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

med_after <- sapply(pbmc@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")], median)

cat(sprintf("过滤后：%d 个基因 × %d 个细胞（去掉 %d 个，占 %.1f%%）\n",
            nrow(pbmc), ncol(pbmc), n_before - ncol(pbmc),
            (n_before - ncol(pbmc)) / n_before * 100))

ggsave(file.path(RES_DIR, "figures/02b_qc_violin_after.png"),
       VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3, pt.size = 0),
       width = 12, height = 4.5, dpi = 150)

qc <- data.frame(
  step             = c("过滤前", "过滤后"),
  cells            = c(n_before, ncol(pbmc)),
  genes            = c(nrow(pbmc), nrow(pbmc)),
  median_nCount    = round(c(med_before["nCount_RNA"],   med_after["nCount_RNA"])),
  median_nFeature  = round(c(med_before["nFeature_RNA"], med_after["nFeature_RNA"])),
  median_percent_mt = round(c(med_before["percent.mt"],  med_after["percent.mt"]), 2),
  row.names = NULL
)
write.csv(qc, file.path(RES_DIR, "qc_summary.csv"), row.names = FALSE)
cat("\n== QC 汇总 ==\n"); print(qc)

saveRDS(pbmc, file.path(RES_DIR, "pbmc3k_qc.rds"))
cat("\n已保存：", file.path(RES_DIR, "pbmc3k_qc.rds"), "\n")
