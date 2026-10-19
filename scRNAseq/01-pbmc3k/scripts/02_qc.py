#!/usr/bin/env python
"""步骤 2（Python / Scanpy）：读入 count matrix，计算 QC 指标并过滤低质量细胞

和 scripts/02_qc.R 是同一件事的两种实现，参数完全对齐：
    nFeature_RNA > 200 且 < 2500，percent.mt < 5

结果写到 results/python/，R 的结果在 results/R/，两边互不覆盖。
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")          # 无图形界面环境（服务器/WSL）必须
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

RES = Path(os.environ.get("RES_DIR", "results/python"))
FIG = RES / "figures"
FIG.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 1

DATA_DIR = Path(os.environ.get("PBMC3K_DIR", "data/filtered_gene_bc_matrices/hg19"))
if not (DATA_DIR / "matrix.mtx").exists():
    sys.exit(f"找不到数据目录：{DATA_DIR}\n请先运行：bash scripts/01_download.sh")

# ---- 2.1 读入 ----
# 和 Seurat 的 Read10X() 读的是同一份文件。
adata = sc.read_10x_mtx(DATA_DIR, var_names="gene_symbols", cache=False)
adata.var_names_make_unique()

# 对齐 Seurat 的 CreateSeuratObject(min.cells = 3, min.features = 200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)

n_genes_raw = adata.n_vars
n_cells_raw = adata.n_obs
print(f"读入完成：{n_genes_raw} 个基因 × {n_cells_raw} 个细胞")

# ---- 2.2 QC 指标 ----
# 对应 Seurat 的 PercentageFeatureSet(pattern = "^MT-")
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                           log1p=False, inplace=True)
# 生成的 obs 列：n_genes_by_counts（= nFeature_RNA）、total_counts（= nCount_RNA）、
#                 pct_counts_mt（= percent.mt）

def violin(adata, cols, path, title):
    sc.pl.violin(adata, cols, jitter=0, multi_panel=True, show=False)
    plt.gcf().suptitle(title, y=1.02)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")

violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
       FIG / "01_qc_violin_before.png", "Before QC filtering")

sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False)
plt.savefig(FIG / "02a_qc_scatter_mt.png", dpi=150, bbox_inches="tight")
plt.close("all")

sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=False)
plt.savefig(FIG / "02b_qc_scatter_genes.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ---- 2.3 过滤 ----
med_before = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].median()

keep = (
    (adata.obs["n_genes_by_counts"] > 200)
    & (adata.obs["n_genes_by_counts"] < 2500)
    & (adata.obs["pct_counts_mt"] < 5)
)
adata = adata[keep].copy()

med_after = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].median()

removed = n_cells_raw - adata.n_obs
print(f"过滤后：{adata.n_vars} 个基因 × {adata.n_obs} 个细胞"
      f"（去掉 {removed} 个，占 {removed / n_cells_raw * 100:.1f}%）")

violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
       FIG / "02c_qc_violin_after.png", "After QC filtering")

# ---- 2.4 汇总 ----
# 列名和 R 版保持一致，这样两个脚本的输出可以直接逐行对比。
qc = pd.DataFrame({
    "step":              ["过滤前", "过滤后"],
    "cells":             [n_cells_raw, adata.n_obs],
    "genes":             [n_genes_raw, adata.n_vars],
    "median_nCount":     [int(med_before["total_counts"]),      int(med_after["total_counts"])],
    "median_nFeature":   [int(med_before["n_genes_by_counts"]), int(med_after["n_genes_by_counts"])],
    "median_percent_mt": [round(med_before["pct_counts_mt"], 2), round(med_after["pct_counts_mt"], 2)],
})
qc.to_csv(RES / "qc_summary.csv", index=False)
print("\n== QC 汇总 ==")
print(qc.to_string(index=False))

adata.write_h5ad(RES / "pbmc3k_qc.h5ad")
print(f"\n已保存：{RES / 'pbmc3k_qc.h5ad'}")
