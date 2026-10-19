#!/usr/bin/env bash
# 步骤 1：从 ENA 下载原始 fastq
source "$(dirname "$0")/lib.sh"

if command -v aria2c >/dev/null 2>&1; then
  DL() { aria2c -x 8 -s 8 -c -d fastq --auto-file-renaming=false "$1"; }
else
  DL() { wget -c -P fastq "$1"; }
fi

tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r run name fq1 fq2; do
  [ -n "${run:-}" ] || continue
  for f in "$fq1" "$fq2"; do
    b="$(basename "$f")"
    if [ -s "fastq/$b" ]; then echo "[skip] $b"; continue; fi
    echo "[get ] $b"
    DL "$f"
  done
done

echo "== 下载完成 =="
du -sh fastq
