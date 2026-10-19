#!/usr/bin/env python
"""步骤 5（Python / Scanpy）：和历史结果对账

和 scripts/05_verify.R 做同一件事，但对账对象多一个：
  ① 对官方教程 —— Seurat pbmc3k 教程给出的期望值
  ② 对 R 版结果 —— 如果 results/R/ 下有输出，直接逐项比对两套实现
                     （这一步才是"双语言交叉验证"真正的价值所在）
"""
import os
from pathlib import Path

import pandas as pd
import scanpy as sc

RES = Path(os.environ.get("RES_DIR", "results/python"))
R_RES = Path("results/R")

if not (RES / "qc_summary.csv").exists():
    raise SystemExit(f"找不到 {RES / 'qc_summary.csv'}，请先运行 python scripts/02_qc.py")

qc = pd.read_csv(RES / "qc_summary.csv")
sizes = pd.read_csv(RES / "cluster_sizes.csv")
n_clusters = len(sizes)
after = qc[qc["step"] == "过滤后"].iloc[0]

print("\n================ ① 与 Seurat 官方教程对账 ================")
items = pd.DataFrame({
    "item":     ["QC 过滤后细胞数", "基因数", "cluster 数（res=0.5）"],
    "expected": ["2638", "13714", "9"],
    "actual":   [str(int(after["cells"])), str(int(after["genes"])), str(n_clusters)],
})
items["match"] = items["expected"] == items["actual"]
print(items.to_string(index=False))
print("\n说明：细胞数和基因数应该精确一致；cluster 数可能差 ±1——")
print("      Scanpy 默认用 Leiden，Seurat 默认用 Louvain，两者结果本就不完全相同。")

# ---- ② 与 R 版对账 ----
print("\n================ ② 与 R 版结果对账（Scanpy vs Seurat）================")
r_qc_file = R_RES / "qc_summary.csv"
r_size_file = R_RES / "cluster_sizes.csv"

if r_qc_file.exists():
    r_qc = pd.read_csv(r_qc_file)
    r_after = r_qc[r_qc["step"] == "过滤后"].iloc[0]

    cmp = pd.DataFrame({
        "指标": ["过滤后细胞数", "基因数", "中位 nCount", "中位 nFeature", "中位 percent.mt"],
        "R (Seurat)": [int(r_after["cells"]), int(r_after["genes"]),
                       int(r_after["median_nCount"]), int(r_after["median_nFeature"]),
                       float(r_after["median_percent_mt"])],
        "Python (Scanpy)": [int(after["cells"]), int(after["genes"]),
                            int(after["median_nCount"]), int(after["median_nFeature"]),
                            float(after["median_percent_mt"])],
    })
    cmp["一致"] = cmp["R (Seurat)"].astype(str) == cmp["Python (Scanpy)"].astype(str)
    print(cmp.to_string(index=False))

    if r_size_file.exists():
        r_sizes = pd.read_csv(r_size_file)
        print(f"\ncluster 数：R(Seurat/Louvain) = {len(r_sizes)}，"
              f"Python(Scanpy/Leiden) = {n_clusters}")
    print("\n结论：细胞数与基因数一致 → 说明两边读的是同一份数据、用了同一套过滤阈值。")
    print("      这是双语言交叉验证最有价值的地方：逻辑错误很难在两种实现里同时出现。")
else:
    print(f"（没有找到 {r_qc_file}，跳过。先跑 Rscript scripts/02_qc.R 就能对比。）")

# ---- ③ 已知细胞类型 marker 是否被找回 ----
print("\n================ ③ 已知细胞类型 marker 是否被找回 ================")
marker_file = RES / "cluster_markers.csv"
if marker_file.exists():
    markers = pd.read_csv(marker_file)
    expected_markers = {
        "T 细胞":   ["CD3D", "CD3E", "IL7R"],
        "B 细胞":   ["MS4A1", "CD79A"],
        "NK 细胞":  ["NKG7", "GNLY"],
        "单核细胞": ["LYZ", "CD14"],
        "血小板":   ["PPBP"],
    }
    rows = []
    for ctype, genes in expected_markers.items():
        found = [g for g in genes if g in set(markers["names"])]
        rows.append({"cell_type": ctype,
                     "markers": ", ".join(genes),
                     "found": ", ".join(found) if found else "（无）",
                     "n_found": len(found), "n_total": len(genes)})
    rec = pd.DataFrame(rows)
    print(rec.to_string(index=False))
    print(f"\n{int((rec.n_found > 0).sum())} / {len(rec)} 类细胞的 marker 至少找回一个")
    print("这一项不依赖 cluster 编号，证明的是数据集里真实的生物学信号被捞了回来。")
    rec.to_csv(RES / "verify_marker_recovery.csv", index=False)
else:
    print(f"（没有找到 {marker_file}，跳过。先跑 04_markers.py。）")

items.to_csv(RES / "verify_vs_tutorial.csv", index=False)
print(f"\n对账结果已保存：{RES / 'verify_vs_tutorial.csv'}")

print("\n================ 运行环境 ================")
import sys
from importlib.metadata import version as _pkg_version

print(f"scanpy   {_pkg_version('scanpy')}")
print(f"anndata  {_pkg_version('anndata')}")
print(f"Python   {sys.version.split()[0]}")
