#!/usr/bin/env bash
# 步骤 2：质控与接头过滤（fastp）+ 汇总报告（MultiQC）
source "$(dirname "$0")/lib.sh"

tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r run name fq1 fq2; do
  [ -n "${run:-}" ] || continue
  [ -s "fastq/$(basename "$fq1")" ] || continue
  echo "== fastp: $run =="
  fastp \
    -i "fastq/$(basename "$fq1")" -I "fastq/$(basename "$fq2")" \
    -o "clean/${run}_R1.clean.fq.gz" -O "clean/${run}_R2.clean.fq.gz" \
    --detect_adapter_for_pe \
    -q 15 -u 40 -l 36 -n 5 \
    -w "$THREADS" \
    -j "qc/fastp/${run}.json" -h "qc/fastp/${run}.html" \
    2> "logs/${run}.fastp.log"
done

# 注意：FastQC 0.12+ 的输出目录参数只接受 --outdir，旧版本的 -o 会报
# "unrecognized arguments: -o"，这里用长参数以兼容新旧版本。
fastqc -t "$THREADS" --outdir qc/fastqc clean/*.fq.gz 2> logs/fastqc.log
multiqc -f -o qc/multiqc qc/fastp qc/fastqc > logs/multiqc.log 2>&1

echo "== 质控完成，报告见 qc/multiqc/multiqc_report.html =="
