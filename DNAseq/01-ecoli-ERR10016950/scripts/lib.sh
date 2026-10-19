#!/usr/bin/env bash
# 公共加载：切到项目根目录、读取配置文件、建好目录
# 用法（在每个脚本开头）： source "$(dirname "$0")/lib.sh"
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

CONFIG="${CONFIG:-config/ecoli.env}"
if [ ! -f "$CONFIG" ]; then
  echo "找不到配置文件：$CONFIG" >&2
  echo "用法：CONFIG=config/ecoli.env bash scripts/03_align.sh" >&2
  exit 1
fi
# shellcheck disable=SC1090
source "$CONFIG"

THREADS="${THREADS:-8}"
REF="${REF:-ref/${REF_NAME}.fa}"
SAMPLES="${SAMPLES:?配置文件里必须设置 SAMPLES}"
# 倍性：二倍体物种为 2（默认）；细菌等单倍体在对应 config 里设 PLOIDY=1。
# 用二倍体模型去调用单倍体基因组，会产生 0/1 杂合基因型（生物学上不可能）
# 和虚高的 SNP 数，同时显著压低 indel 检出。
PLOIDY="${PLOIDY:-2}"

export THREADS CONFIG SAMPLES REF PLOIDY
mkdir -p fastq clean bam vcf ref logs qc/fastp qc/fastqc

echo ">>> 配置：$CONFIG  样本表：$SAMPLES  参考：$REF  线程：$THREADS"
