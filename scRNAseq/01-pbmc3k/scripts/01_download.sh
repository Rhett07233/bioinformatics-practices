#!/usr/bin/env bash
# 步骤 1：下载 10x Genomics 的 PBMC 3k 数据并解压
#
# 这是 Seurat 官方入门教程使用的数据集：
#   https://satijalab.org/seurat/articles/pbmc3k_tutorial
# 来源：10x Genomics 公开数据，2700 个 PBMC 单细胞，Illumina NextSeq 500。
#
# 注意：单细胞分析**不需要**从 fastq 开始做比对——
# CellRanger 已经把 reads 比对并计数好了，直接给出「基因 × 细胞」的 count matrix。
# 所以这一步只是下载一个很小的矩阵文件，不是几十 GB 的原始数据。
set -euo pipefail
cd "$(dirname "$0")/.."

URL="${PBMC3K_URL:-https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz}"
TARBALL="data/pbmc3k_filtered_gene_bc_matrices.tar.gz"
DEST="data/filtered_gene_bc_matrices/hg19"

mkdir -p data

if [ -s "$DEST/matrix.mtx" ]; then
  echo "[skip] 数据已存在：$DEST"
else
  if [ ! -s "$TARBALL" ]; then
    echo "== 下载 PBMC 3k（压缩包只有几 MB）=="
    wget -c -O "$TARBALL" "$URL"
  fi
  echo "== 解压 =="
  tar -xzf "$TARBALL" -C data
fi

echo
echo "== 数据就绪，目录内容 =="
ls -1 "$DEST"
echo
echo "下一步：Rscript scripts/02_qc.R"
