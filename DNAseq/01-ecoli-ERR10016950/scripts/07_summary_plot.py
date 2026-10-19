#!/usr/bin/env python
"""07_summary_plot.py —— 从已有结果生成一张 QC / 变异检测汇总图

只读取已有文件，不改动任何中间产物：
    qc/fastp/<run>.json          reads 过滤前后统计（fastp）
    logs/<run>.flagstat.txt      比对统计（samtools flagstat）
    vcf/<run>.stats.txt          变异数与 Ti/Tv（bcftools stats）

输出：
    qc/summary/<run>.summary.png

用法（在项目根目录）：
    conda activate dna
    python scripts/07_summary_plot.py

依赖：matplotlib（已包含在 envs/dna.yaml 中；老环境执行
      `conda env update -f envs/dna.yaml` 补装）。

    # 指定 run（默认取 qc/fastp/ 下第一个 json 对应的名字）
    RUN=ERR10016950 python scripts/07_summary_plot.py

注意：本脚本只需要 matplotlib，不需要 bwa/samtools/bcftools，
      所以放在哪个 Python 环境里跑都可以，不一定要用 dna 环境。
"""
import json
import os
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # 无图形界面（服务器 / WSL）必须
import matplotlib.pyplot as plt

RES = Path(os.environ.get("RES_DIR", "."))
QC_DIR = RES / "qc"
LOG_DIR = RES / "logs"
VCF_DIR = RES / "vcf"
OUT_DIR = QC_DIR / "summary"

MUTED = "#4C72B0"
WARN = "#C44E52"
LIGHT = "#8C8C8C"


def die(msg):
    sys.exit(msg)


# ---------- 定位 run ----------
run = os.environ.get("RUN", "").strip()
if not run:
    cands = sorted((QC_DIR / "fastp").glob("*.json")) if (QC_DIR / "fastp").is_dir() else []
    if not cands:
        die(f"找不到 {QC_DIR / 'fastp'}/*.json —— 请先跑完 02_qc.sh")
    run = cands[0].stem
print(f">>> run = {run}")


def read_json(p):
    with open(p) as fh:
        return json.load(fh)


def one_int(pattern, text, default=None):
    m = re.search(pattern, text, re.M)
    return int(m.group(1)) if m else default


def one_float(pattern, text, default=None):
    m = re.search(pattern, text, re.M)
    return float(m.group(1)) if m else default


# ---------- 读取三个来源 ----------
fastp_p = QC_DIR / "fastp" / f"{run}.json"
flag_p = LOG_DIR / f"{run}.flagstat.txt"
stats_p = VCF_DIR / f"{run}.stats.txt"

for p in (fastp_p, flag_p, stats_p):
    if not p.exists():
        print(f"⚠️  缺少 {p}（对应图格将留空）", file=sys.stderr)

fastp = read_json(fastp_p) if fastp_p.exists() else {}
flag = flag_p.read_text(errors="replace") if flag_p.exists() else ""
stats = stats_p.read_text(errors="replace") if stats_p.exists() else ""

bf = fastp.get("summary", {}).get("before_filtering", {})
af = fastp.get("summary", {}).get("after_filtering", {})

# flagstat 每行形如 "3532755 + 0 in total (QC-passed reads + QC-failed reads)"
total_flag = one_int(r"^(\d+) \+ \d+ in total", flag, 0)
mapped = one_int(r"^(\d+) \+ \d+ mapped \(", flag, 0)
mapped_pct = one_float(r"^\d+ \+ \d+ mapped \(([\d.]+)%", flag)
proper_pct = one_float(r"properly paired \(([\d.]+)%", flag)
dups = one_int(r"^(\d+) \+ \d+ duplicates", flag, 0)

# bcftools stats：SN/TSTV 行以制表符分隔
n_snp = one_int(r"number of SNPs:\s+(\d+)", stats)
n_indel = one_int(r"number of indels:\s+(\d+)", stats)
n_rec = one_int(r"number of records:\s+(\d+)", stats)
titv = one_float(r"^TSTV\t\d+\t\d+\t\d+\t([\d.]+)", stats)

print(f"    reads 过滤前/后 : {bf.get('total_reads', 'NA')} / {af.get('total_reads', 'NA')}")
print(f"    Q30 过滤前/后   : {bf.get('q30_rate', 'NA')} / {af.get('q30_rate', 'NA')}")
print(f"    mapped          : {mapped_pct}%  properly paired: {proper_pct}%  重复: {dups}")
print(f"    SNP / indel / records : {n_snp} / {n_indel} / {n_rec}   Ti/Tv: {titv}")

# ---------- 画图 ----------
fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.6), dpi=150)
fig.suptitle(f"{run} — QC & variant calling summary", fontsize=14, fontweight="bold")

# ① reads 数量：过滤前 vs 过滤后
ax = axes[0][0]
if bf and af:
    vals = [bf["total_reads"] / 1e6, af["total_reads"] / 1e6]
    bars = ax.bar(["Before\nfiltering", "After\nfiltering"], vals, color=[LIGHT, MUTED], width=0.55)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.2f}M", ha="center", va="bottom", fontsize=10)
    ax.set_ylim(0, max(vals) * 1.18)
ax.set_ylabel("Reads (millions)")
ax.set_title("Reads passing fastp", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)

# ② 比对率：堆叠条形
ax = axes[0][1]
if mapped_pct is not None:
    ax.barh(["Aligned"], [mapped_pct], color=MUTED, height=0.45, label="mapped")
    ax.barh(["Aligned"], [100 - mapped_pct], left=[mapped_pct], color="#DDDDDD", height=0.45,
            label="unmapped")
    ax.text(mapped_pct / 2, 0, f"{mapped_pct:.2f}%", ha="center", va="center",
            color="white", fontsize=11, fontweight="bold")
    ax.text(mapped_pct + (100 - mapped_pct) / 2, 0, f"{100 - mapped_pct:.2f}%",
            ha="center", va="center", fontsize=10)
    if proper_pct is not None:
        ax.text(0.5, -0.42, f"properly paired: {proper_pct:.2f}%", ha="center",
                fontsize=9, color="#555555")
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.7, 0.7)
    ax.set_xlabel("Percent of reads")
ax.set_title("Alignment (BWA-MEM → K-12 MG1655)", fontsize=11)
ax.spines[["top", "right", "left"]].set_visible(False)
ax.set_yticks([])

# ③ 变异数量（对数轴）
ax = axes[1][0]
if n_snp is not None:
    vals = [n_snp, n_indel or 0]
    bars = ax.bar(["SNPs", "Indels"], vals, color=[MUTED, WARN], width=0.5)
    ax.set_yscale("log")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:,}", ha="center", va="bottom", fontsize=10)
    ax.set_ylim(1, max(vals) * 4)
    if titv is not None:
        ax.text(0.98, 0.95, f"Ti/Tv = {titv:.2f}", transform=ax.transAxes, ha="right", va="top",
                fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="#F2F2F2", ec="#CCCCCC"))
ax.set_ylabel("Count (log scale)")
ax.set_title("Variants after filtering (ploidy = 1)", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)

# ④ Q30 比例
ax = axes[1][1]
if bf and af:
    vals = [bf.get("q30_rate", 0) * 100, af.get("q30_rate", 0) * 100]
    bars = ax.bar(["Before\nfiltering", "After\nfiltering"], vals, color=[LIGHT, MUTED], width=0.55)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f}%", ha="center", va="bottom", fontsize=10)
    ax.set_ylim(min(vals) * 0.9, 100)
ax.set_ylabel("Q30 bases (%)")
ax.set_title("Base quality (Q30)", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)

fig.text(
    0.5, 0.012,
    "Note: the sample's strain metadata is empty and no matching closed genome exists; "
    "the ~21% unmapped reads are most likely strain-specific accessory sequence (see README §8.5).",
    ha="center", fontsize=8, color="#666666",
)

fig.tight_layout(rect=(0, 0.035, 1, 0.94))
OUT_DIR.mkdir(parents=True, exist_ok=True)
out = OUT_DIR / f"{run}.summary.png"
fig.savefig(out)
print(f"\n已保存：{out}")
