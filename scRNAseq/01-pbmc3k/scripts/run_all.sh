#!/usr/bin/env bash
# 一键跑完整个单细胞分析流程
#
#   conda activate scrna
#   bash scripts/run_all.sh
#
# 全程在笔记本上几分钟到十几分钟，不需要服务器。
# 如果想换聚类分辨率或 PC 数（03_cluster.R 的两个关键参数）：
#   NPCS=15 RES=0.8 bash scripts/run_all.sh
set -euo pipefail
cd "$(dirname "$0")/.."

echo "########## 1/5 下载数据 ##########";     bash scripts/01_download.sh
echo "########## 2/5 质控与过滤 ##########";   Rscript scripts/02_qc.R
echo "########## 3/5 降维与聚类 ##########";   Rscript scripts/03_cluster.R
echo "########## 4/5 marker 与注释 ##########"; Rscript scripts/04_markers.R
echo "########## 5/5 与官方教程对账 ##########"; Rscript scripts/05_verify.R

echo
echo "== 全部完成 =="
echo "结果表：results/R/*.csv"
echo "结果图："
ls -1 results/R/figures/ 2>/dev/null || true
