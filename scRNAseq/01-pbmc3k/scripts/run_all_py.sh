#!/usr/bin/env bash
# 一键跑完整个单细胞分析流程（Python / Scanpy 版）
#
#   conda activate scanpy
#   bash scripts/run_all_py.sh
#
# 和 run_all.sh（R / Seurat 版）跑的是同一份数据、同一套阈值。
# 两边都跑一遍，05_verify 会把结果逐项对比——这才是"双语言"的意义。
set -euo pipefail
cd "$(dirname "$0")/.."

echo "########## 1/5 下载数据 ##########";       bash scripts/01_download.sh
echo "########## 2/5 质控与过滤 ##########";     python scripts/02_qc.py
echo "########## 3/5 降维与聚类 ##########";     python scripts/03_cluster.py
echo "########## 4/5 marker 与注释 ##########"; python scripts/04_markers.py
echo "########## 5/5 与官方教程对账 ##########"; python scripts/05_verify.py

echo
echo "== 全部完成 =="
echo "结果表：results/python/*.csv"
echo "结果图："
ls -1 results/python/figures/ 2>/dev/null || true
