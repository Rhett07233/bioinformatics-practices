#!/usr/bin/env python
"""步骤 3（Python / Scanpy）：归一化 → 高变基因 → PCA → 聚类 → UMAP

和 scripts/03_cluster.R 做同一件事。

两处刻意的差异（不是 bug，是工具惯例不同）：
  1. **邻接图邻居数**：Seurat 的 FindNeighbors() 默认 k.param = 20，
     Scanpy 的 sc.pp.neighbors() 默认 n_neighbors = 15。这里显式设成 20 向 Seurat 对齐。
  2. **聚类算法**：Seurat 默认 Louvain；Scanpy 现在默认 Leiden（Louvain 的改进版，
     保证社区内部连通）。所以两边的 cluster 数不保证完全一致，差 ±1 属于正常。
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

RES = Path(os.environ.get("RES_DIR", "results/python"))
FIG = RES / "figures"
FIG.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(RES / "pbmc3k_qc.h5ad")

# ---- 3.1 归一化 ----
# 对应 Seurat 的 NormalizeData(scale.factor = 10000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ---- 3.2 高变基因 ----
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
print(f"高变基因：{int(adata.var.highly_variable.sum())} 个")
print("Top 10 高变基因：",
      ", ".join(adata.var_names[adata.var.highly_variable][:10]))

sc.pl.highly_variable_genes(adata, show=False)
plt.savefig(FIG / "03_variable_features.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ---- 3.3 只保留高变基因，标准化后做 PCA ----
# 关键一步：先把「全基因的 log-normalized 数据」存进 .raw，
# 后面找 marker、画图都从 .raw 取，不会被 scale 后的 z-score 污染。
adata.raw = adata
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)

NPCS = int(os.environ.get("NPCS", "10"))
RESOLUTION = float(os.environ.get("RES", "0.5"))
N_COMP_PCA = max(50, NPCS + 10)

sc.tl.pca(adata, n_comps=N_COMP_PCA, svd_solver="arpack")

vr = adata.uns["pca"]["variance_ratio"] * 100
pd.DataFrame({
    "PC": np.arange(1, len(vr) + 1),
    "pct_variance": np.round(vr, 3),
    "cum_pct": np.round(np.cumsum(vr), 3),
}).to_csv(RES / "pca_variance.csv", index=False)

print("前 10 个 PC 解释的方差（%）：", ", ".join(f"{v:.2f}" for v in vr[:10]))
print(f"前 10 个 PC 累计解释 {vr[:10].sum():.1f}%；前 20 个累计 {vr[:20].sum():.1f}%")

# 肘部图（对应 Seurat 的 ElbowPlot）
sc.pl.pca_variance_ratio(adata, n_pcs=min(30, len(vr)), log=True, show=False)
plt.savefig(FIG / "04_elbow.png", dpi=150, bbox_inches="tight")
plt.close("all")

sc.pl.pca(adata, show=False)
plt.savefig(FIG / "05_pca.png", dpi=150, bbox_inches="tight")
plt.close("all")

sc.pl.pca_loadings(adata, components="1,2", show=False)
plt.savefig(FIG / "06_pca_loadings.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ---- 3.4 邻居图 + Leiden 聚类 + UMAP ----
print(f"\n== 使用 n_pcs = {NPCS}，resolution = {RESOLUTION}，n_neighbors = 20 ==")
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=NPCS)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=RESOLUTION, key_added="leiden")

sc.pl.umap(adata, color="leiden", legend_loc="on data", frameon=False,
           title="PBMC 3k — Leiden clusters", show=False)
plt.savefig(FIG / "07_umap_clusters.png", dpi=150, bbox_inches="tight")
plt.close("all")

counts = adata.obs["leiden"].value_counts().sort_index()
print(f"\n共 {len(counts)} 个 cluster：")
print(counts.to_string())

pd.DataFrame({"cluster": counts.index, "n_cells": counts.values.astype(int)}) \
  .to_csv(RES / "cluster_sizes.csv", index=False)

adata.write_h5ad(RES / "pbmc3k_clustered.h5ad")
print(f"\n已保存：{RES / 'pbmc3k_clustered.h5ad'}")
