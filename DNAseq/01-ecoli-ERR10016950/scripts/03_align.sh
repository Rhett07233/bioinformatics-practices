#!/usr/bin/env bash
# 步骤 3：建索引 + BWA-MEM 比对 + 排序 + 去重
#
# 这里是 DNA 流程和 RNA-seq 流程最不一样的地方：
#   - 比对用 BWA-MEM（DNA 比对器，不使用剪接）
#   - 去重是标准步骤（PCR 重复会造成假阳性变异），用 samtools markdup
#   - markdup 要求先做 fixmate，所以顺序必须是：名字排序 → fixmate → 坐标排序 → markdup
source "$(dirname "$0")/lib.sh"

# ---------- 3.1 参考基因组 ----------
if [ ! -s "$REF" ]; then
  if [ -z "${REF_URL:-}" ]; then
    echo "配置文件里没有设置 REF_URL" >&2
    exit 1
  fi
  echo "== 下载参考基因组 =="
  wget -c -O "${REF}.gz" "$REF_URL"
  gunzip -c "${REF}.gz" > "$REF"
fi
samtools faidx "$REF"
[ -s "${REF}.bwt" ] || { echo "== bwa index（一次性） =="; bwa index "$REF" 2> logs/bwa-index.log; }

# ---------- 3.2 比对与去重 ----------
tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r run name fq1 fq2; do
  [ -n "${run:-}" ] || continue
  [ -s "clean/${run}_R1.clean.fq.gz" ] || continue
  echo "== bwa mem: $run =="

  bwa mem -t "$THREADS" \
      -R "@RG\tID:${run}\tSM:${name}\tLB:${name}\tPL:ILLUMINA\tPU:${run}" \
      "$REF" "clean/${run}_R1.clean.fq.gz" "clean/${run}_R2.clean.fq.gz" \
      2> "logs/${run}.bwa.log" \
    | samtools sort -n -@ "$THREADS" -m 1G -o "bam/${run}.namesort.bam" -

  samtools fixmate -m -@ "$THREADS" "bam/${run}.namesort.bam" "bam/${run}.fixmate.bam"
  samtools sort -@ "$THREADS" -m 1G -o "bam/${run}.sorted.bam" "bam/${run}.fixmate.bam"
  samtools markdup -@ "$THREADS" -s \
      "bam/${run}.sorted.bam" "bam/${run}.markdup.bam" 2> "logs/${run}.markdup.log"
  samtools index -@ "$THREADS" "bam/${run}.markdup.bam"
  samtools flagstat -@ "$THREADS" "bam/${run}.markdup.bam" > "logs/${run}.flagstat.txt"

  rm -f "bam/${run}.namesort.bam" "bam/${run}.fixmate.bam" "bam/${run}.sorted.bam"
done

echo "== 比对完成，比对率与去重情况： =="
grep -H "mapped (" logs/*.flagstat.txt || true
grep -H "DUPLICATE" logs/*.markdup.log 2>/dev/null | head -5 || true
