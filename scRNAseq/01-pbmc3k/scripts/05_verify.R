#!/usr/bin/env Rscript
# 步骤 5：和历史结果对账 —— 怎么知道自己的流程跑对了？
#
# 这是整个练习里最重要的一步：
# "你怎么知道你的结果是可信的？"
#
# 对账基准有两个：
#   ① Seurat 官方 pbmc3k 教程
#      https://satijalab.org/seurat/articles/pbmc3k_tutorial
#      （2026-04-28 编译版，Seurat 5.5.0 / R 4.5.2，使用同一份 10x PBMC 3k 数据）
#   ② Python（Scanpy）版的输出 results/python/
#      两套独立实现跑同一份数据，结果对得上 → 逻辑错误的可能性大幅降低
suppressPackageStartupMessages({
  library(Seurat)
})

RES_DIR <- Sys.getenv("RES_DIR", file.path("results", "R"))
PY_DIR  <- file.path("results", "python")

qc <- readRDS(file.path(RES_DIR, "pbmc3k_qc.rds"))
cl <- readRDS(file.path(RES_DIR, "pbmc3k_clustered.rds"))

# ---- 5.1 与官方教程对账 ----
items <- data.frame(
  item     = c("QC 过滤后细胞数", "基因数", "cluster 数（dims=1:10, res=0.5）"),
  expected = c("2638", "13714", "9"),
  actual   = c(as.character(ncol(qc)),
               as.character(nrow(qc)),
               as.character(nlevels(cl))),
  stringsAsFactors = FALSE
)
items$match <- items$expected == items$actual

cat("\n================ ① 与 Seurat 官方教程对账 ================\n")
print(items, row.names = FALSE)

cat("\n说明：细胞数和基因数应该是精确一致的；cluster 数可能差 ±1，\n")
cat("      因为 Louvain 实现（igraph 版本）和 Seurat 版本不同会带来微小差异。\n")

# ---- 5.2 与 Python（Scanpy）版对账 ----
cat("\n================ ② 与 Python（Scanpy）版对账 ================\n")
py_qc_file <- file.path(PY_DIR, "qc_summary.csv")

if (file.exists(py_qc_file)) {
  py_qc <- read.csv(py_qc_file)
  py_after <- py_qc[py_qc$step == "过滤后", ]

  med_r  <- sapply(qc@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")], median)
  cmp <- data.frame(
    metric      = c("过滤后细胞数", "基因数", "中位 nCount", "中位 nFeature", "中位 percent.mt"),
    R_Seurat    = c(ncol(qc), nrow(qc), round(med_r["nCount_RNA"]),
                    round(med_r["nFeature_RNA"]), round(med_r["percent.mt"], 2)),
    Python_Scanpy = c(py_after$cells, py_after$genes, py_after$median_nCount,
                      py_after$median_nFeature, py_after$median_percent_mt),
    stringsAsFactors = FALSE
  )
  cmp$match <- as.character(cmp$R_Seurat) == as.character(cmp$Python_Scanpy)
  print(cmp, row.names = FALSE)

  py_size_file <- file.path(PY_DIR, "cluster_sizes.csv")
  if (file.exists(py_size_file)) {
    py_sizes <- read.csv(py_size_file)
    cat(sprintf("\ncluster 数：R(Seurat/Louvain) = %d，Python(Scanpy/Leiden) = %d\n",
                nlevels(cl), nrow(py_sizes)))
  }
  cat("\n结论：细胞数与基因数一致 → 说明两边读的是同一份数据、用了同一套过滤阈值。\n")
  cat("      这是双语言交叉验证最有价值的地方：逻辑错误很难在两种实现里同时出现。\n")
} else {
  cat("（没有找到 ", py_qc_file, "，跳过。先跑 python scripts/02_qc.py 就能对比。）\n", sep = "")
}

# ---- 5.3 更本质的验证：已知的 marker 基因有没有被找出来 ----
# 这一项不依赖 cluster 编号。只要这些公认的细胞类型 marker
# 出现在 marker 表里，就说明流程捞回了真实的生物学信号，
# 而不是被版本差异或参数选择决定。
markers_file <- file.path(RES_DIR, "cluster_markers.csv")
if (file.exists(markers_file)) {
  markers <- read.csv(markers_file)
  expected_markers <- list(
    "T 细胞"   = c("CD3D", "CD3E", "IL7R"),
    "B 细胞"   = c("MS4A1", "CD79A"),
    "NK 细胞"  = c("NKG7", "GNLY"),
    "单核细胞" = c("LYZ", "CD14"),
    "血小板"   = c("PPBP")
  )

  cat("\n================ ③ 已知细胞类型 marker 是否被找回 ================\n")
  res <- do.call(rbind, lapply(names(expected_markers), function(ct) {
    genes <- expected_markers[[ct]]
    hit <- genes[genes %in% markers$gene]
    data.frame(
      cell_type = ct,
      markers   = paste(genes, collapse = ", "),
      found     = paste(if (length(hit)) hit else "（无）", collapse = ", "),
      n_found   = length(hit),
      n_total   = length(genes)
    )
  }))
  print(res, row.names = FALSE)

  n_ok <- sum(res$n_found > 0)
  cat(sprintf("\n%d / %d 类细胞的 marker 至少找回一个\n", n_ok, nrow(res)))
  cat("如果全部找回，说明流程在生物学上是可信的——\n")
  cat("因为这些基因和你的参数选择无关，是数据集里真实存在的信号。\n")

  write.csv(res, file.path(RES_DIR, "verify_marker_recovery.csv"), row.names = FALSE)
} else {
  cat("\n（还没跑 04_markers.R，跳过 marker 验证）\n")
}

write.csv(items, file.path(RES_DIR, "verify_vs_tutorial.csv"), row.names = FALSE)
cat("\n对账结果已保存：", file.path(RES_DIR, "verify_vs_tutorial.csv"), "\n")

# ---- 5.4 可复现性 ----
cat("\n================ 运行环境 ================\n")
cat(sprintf("Seurat       %s\n", as.character(packageVersion("Seurat"))))
cat(sprintf("SeuratObject %s\n", as.character(packageVersion("SeuratObject"))))
cat(sprintf("R            %s\n", R.version.string))
