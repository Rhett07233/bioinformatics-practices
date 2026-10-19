#!/usr/bin/env python
"""步骤 4（Python / Scanpy）：marker 基因 + 细胞类型注释

和 scripts/04_markers.R 做同一件事。

列名对照（两个脚本的输出可以这样对上）：
    Scanpy                    Seurat / FindAllMarkers
    names                  →  gene
    logfoldchanges         →  avg_log2FC
    pct_nz_group           →  pct.1
    pct_nz_reference       →  pct.2
    pvals_adj              →  p_val_adj
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

adata = sc.read_h5ad(RES / "pbmc3k_clustered.h5ad")
GROUP = "leiden"

# ---- 4.1 每个 cluster vs 其余所有细胞的差异表达 ----
# use_raw=True：从 .raw 取全基因的 log-normalized 数据，
# 而不是被 scale 成 z-score 的 .X。这一步很容易搞错，结果会全歪。
# pts=True 必须显式传：新版本 scanpy 已把 pts 的默认值改为 False（本环境 1.12.4 实测）。
# 不传的话 adata.uns["rank_genes_groups"] 里不会有 pts / pts_rest，
# sc.get.rank_genes_groups_df() 就不会生成 pct_nz_group / pct_nz_reference 两列，
# 下游一引用就是 KeyError: 'pct_nz_group'。
sc.tl.rank_genes_groups(adata, GROUP, method="wilcoxon", use_raw=True, pts=True)

markers = sc.get.rank_genes_groups_df(adata, group=None)
markers = markers.dropna(subset=["logfoldchanges"])

# 防御：列名随 scanpy 版本变过，缺了就早点抛出可读信息，而不是一句 KeyError
for _col in ("pct_nz_group", "pct_nz_reference"):
    if _col not in markers.columns:
        raise RuntimeError(
            f"找不到列 {_col!r}；实际列为 {list(markers.columns)}。"
            "多半是 rank_genes_groups() 没传 pts=True。"
        )

markers.to_csv(RES / "cluster_markers_all.csv", index=False)

# 只保留正向上的、且在该 cluster 里有一定检出比例的基因
sig = markers[(markers["logfoldchanges"] > 0.25) & (markers["pct_nz_group"] > 0.25)].copy()
sig.to_csv(RES / "cluster_markers.csv", index=False)
print(f"共 {len(sig)} 条 marker 记录，覆盖 {sig['group'].nunique()} 个 cluster")

top10 = (sig.sort_values("scores", ascending=False)
            .groupby("group", observed=True).head(10))
top10.to_csv(RES / "cluster_top10_markers.csv", index=False)

print("\n== 每个 cluster 的 top5 marker ==")
show = (top10.groupby("group", observed=True).head(5)
        [["group", "names", "logfoldchanges", "pct_nz_group", "pct_nz_reference", "pvals_adj"]])
show.columns = ["cluster", "gene", "log2FC", "pct.1", "pct.2", "p_val_adj"]
print(show.to_string(index=False))

# ---- 4.2 经典 marker 的表达模式 ----
canonical = ["IL7R", "CCR7", "CD14", "LYZ", "S100A4", "MS4A1", "CD8A",
             "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP"]

sc.pl.dotplot(adata, canonical, groupby=GROUP, use_raw=True, show=False)
plt.savefig(FIG / "09_dotplot_canonical.png", dpi=150, bbox_inches="tight")
plt.close("all")

# 每个 cluster 里这些基因的平均表达（人工判断细胞类型用）
raw = adata.raw.to_adata()
present = [g for g in canonical if g in raw.var_names]
sub = raw[:, present].to_df()
avg = sub.groupby(adata.obs[GROUP].values).mean()
avg.to_csv(RES / "canonical_marker_avg_expression.csv")
print("\n== 经典 marker 在各 cluster 的平均表达 ==")
print(avg.round(2).to_string())

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, use_raw=True,
                                groupby=GROUP, show=False)
plt.savefig(FIG / "08_heatmap_top10.png", dpi=150, bbox_inches="tight")
plt.close("all")

# ---- 4.3 注释细胞类型 ----
# 官方教程（resolution=0.5）给出的对应关系：
#   0 Naive CD4 T / 1 CD14+ Mono / 2 Memory CD4 T / 3 B / 4 CD8 T /
#   5 FCGR3A+ Mono / 6 NK / 7 DC / 8 Platelet
#
# 注意：Leiden 的编号和 Louvain 不一定一致，cluster 数也可能差 1 个。
# 所以这里不硬套，先检查数量对不对；不对就提示人工判断。
new_ids = ["Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
           "FCGR3A+ Mono", "NK", "DC", "Platelet"]
clusters = sorted(adata.obs[GROUP].unique(), key=lambda x: int(x))

if len(clusters) == len(new_ids):
    mapping = dict(zip(clusters, new_ids))
    adata.obs["cell_type"] = adata.obs[GROUP].map(mapping)

    sc.pl.umap(adata, color="cell_type", legend_loc="on data", frameon=False,
               title="PBMC 3k — cell type annotation", show=False)
    plt.savefig(FIG / "10_umap_annotated.png", dpi=150, bbox_inches="tight")
    plt.close("all")

    print("\n== 注释后的细胞类型分布 ==")
    print(adata.obs["cell_type"].value_counts().to_string())

    adata.write_h5ad(RES / "pbmc3k_final.h5ad")
    print(f"\n已保存：{RES / 'pbmc3k_final.h5ad'}")
else:
    print(f"\n注意：cluster 数是 {len(clusters)}，不是教程里的 9 个，跳过自动注释。")
    print("请对照 canonical_marker_avg_expression.csv 人工判断每个 cluster 的类型。")
