#!/usr/bin/env bash
# 一键跑完全流程
#
#   bash scripts/run_all.sh
#
# 默认使用 config/ecoli.env。换数据集时用 CONFIG 指定：
#   CONFIG=config/其他.env bash scripts/run_all.sh
#
# 注意：CONFIG 必须 export，子脚本才会读到同一个配置文件。
set -euo pipefail
cd "$(dirname "$0")/.."

export CONFIG="${CONFIG:-config/ecoli.env}"
echo "########## 使用配置：$CONFIG ##########"

echo "########## 1/5 下载 ##########";  bash scripts/01_download.sh
echo "########## 2/5 质控 ##########";  bash scripts/02_qc.sh
echo "########## 3/5 比对 ##########";  bash scripts/03_align.sh
echo "########## 4/5 变异检测 ##########"; bash scripts/04_call.sh
echo "########## 5/5 注释 ##########";  bash scripts/05_annotate.sh

echo
echo "== 全部完成 =="
ls -1 vcf/ 2>/dev/null || true
