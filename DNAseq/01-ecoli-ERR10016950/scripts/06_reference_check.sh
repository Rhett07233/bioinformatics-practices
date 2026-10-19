#!/usr/bin/env bash
# 步骤 6（诊断）：参考基因组筛选
#
# 背景：本数据集（ERR10016950，真实 E. coli，约 123×）比对到 K-12 MG1655 只有
# 78.57%，而 E. coli 通常应有 >95%。已知 E. coli 泛基因组差异极大
# （O157:H7 相对 K-12 有约 1.4 Mb 菌株特有序列），所以首要怀疑是**参考株不匹配**，
# 而不是数据质量或流程实现有问题。
#
# 本脚本用**同一份 reads** 比对 6 个代表性 E. coli 完整基因组，比较比对率：
#   - 若某个参考明显更高        → 参考株不匹配，换参考即可
#   - 若「6 个合并」参考能到 ~99% → 确认为菌株特有序列，不是技术缺陷
#   - 若全都停在 78–85% 且互相接近 → 未比对部分并非 E. coli 参考差异，
#                                 提示样本混合/污染（本项目后来实测证实混入 K. pneumoniae）
#
# 为控制耗时，先随机抽 20 万对 reads（seqtk 同种子抽样，R1/R2 严格配对）。
#
# 用法：
#   CONFIG=config/ecoli.env bash scripts/06_reference_check.sh
#   SUB_READS=500000 CONFIG=config/ecoli.env bash scripts/06_reference_check.sh
source "$(dirname "$0")/lib.sh"

OUT="qc/refcheck"
REFDIR="ref/screen"
SUB_READS="${SUB_READS:-200000}"
SEED="${SEED:-42}"
mkdir -p "$OUT" "$REFDIR" logs

command -v seqtk >/dev/null 2>&1 || { echo "缺少 seqtk：conda install -c bioconda seqtk" >&2; exit 1; }

# ---------- 6.1 取第一个样本的 run 名 ----------
run="$(awk 'NR==2{print $1; exit}' "$SAMPLES")"
[ -n "$run" ] || { echo "样本表为空：$SAMPLES" >&2; exit 1; }
[ -s "clean/${run}_R1.clean.fq.gz" ] || {
  echo "找不到 clean/${run}_R1.clean.fq.gz —— 请先跑完 02_qc.sh" >&2; exit 1; }
echo ">>> 样本 run = $run，抽样 $SUB_READS 对（seed=$SEED）"

# ---------- 6.2 抽样子集（R1/R2 同种子，保持配对） ----------
if [ ! -s "$OUT/sub_R1.fq.gz" ] || [ ! -s "$OUT/sub_R2.fq.gz" ]; then
  echo "== 抽样 =="
  seqtk sample -s "$SEED" "clean/${run}_R1.clean.fq.gz" "$SUB_READS" | gzip > "$OUT/sub_R1.fq.gz"
  seqtk sample -s "$SEED" "clean/${run}_R2.clean.fq.gz" "$SUB_READS" | gzip > "$OUT/sub_R2.fq.gz"
fi
n_sub=$(zcat "$OUT/sub_R1.fq.gz" | wc -l)
echo ">>> 抽样后 R1 行数 = $n_sub（应为 $((${SUB_READS} * 4))）"

# ---------- 6.3 候选参考（均为 NCBI RefSeq 完整基因组） ----------
# 格式：短名|GCF号|组装名
REFS="
K-12_MG1655|GCF_000005845.2|ASM584v2
O157H7_Sakai|GCF_000008865.2|ASM886v2
CFT073|GCF_000007445.1|ASM744v1
UTI89|GCF_000013265.1|ASM1326v1
Nissle1917|GCF_000714595.1|ASM71459v1
ETEC_H10407|GCF_000210475.1|ASM21047v1
"

echo "== 下载参考基因组（每个约 1.4 MB 压缩） =="
echo "$REFS" | while IFS='|' read -r short acc asmname; do
  [ -n "${short:-}" ] || continue
  fa="$REFDIR/${short}.fa"
  [ -s "$fa" ] && { echo "  [skip] $short"; continue; }

  d="${acc#GCF_}"; d="${d%%.*}"
  url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${d:0:3}/${d:3:3}/${d:6:3}/${acc}_${asmname}/${acc}_${asmname}_genomic.fna.gz"

  echo "  [get ] $short  <- $acc"
  if wget -q -c -O "${fa}.gz" "$url"; then
    gunzip -c "${fa}.gz" > "$fa"
    samtools faidx "$fa"
  else
    echo "    下载失败，跳过：$short" >&2
    rm -f "${fa}.gz"
  fi
done

# ---------- 6.4 合并参考（模拟「泛基因组」） ----------
combined="$REFDIR/all_combined.fa"
if [ ! -s "$combined" ]; then
  echo "== 构建合并参考 =="
  for fa in "$REFDIR"/*.fa; do
    [ -s "$fa" ] || continue
    b="$(basename "$fa" .fa)"
    awk -v p="$b" '/^>/{sub(/^>/, ">" p "|")}1' "$fa"
  done > "$combined"
  samtools faidx "$combined"
fi

# ---------- 6.5 逐个比对，统计比对率 ----------
: > "$OUT/reference_screening.tsv"
printf 'reference\tmapped_pct\tprimary_mapped_pct\n' >> "$OUT/reference_screening.tsv"

map_one() {
  local tag="$1" fa="$2"
  [ -s "$fa" ] || return 0
  echo "== bwa mem: $tag =="
  [ -s "${fa}.bwt" ] || bwa index "$fa" 2> "logs/refcheck_${tag}.index.log"

  bwa mem -t "$THREADS" "$fa" "$OUT/sub_R1.fq.gz" "$OUT/sub_R2.fq.gz" \
      2> "logs/refcheck_${tag}.bwa.log" \
    | samtools view -@ "$THREADS" -b - \
    | samtools flagstat -@ "$THREADS" - > "$OUT/${tag}.flagstat.txt"

  local mp pm
  mp=$(sed -n 's/.*[^a-z]mapped (\([0-9.]*\)%.*/\1/p' "$OUT/${tag}.flagstat.txt" | head -1)
  pm=$(sed -n 's/.*primary mapped (\([0-9.]*\)%.*/\1/p' "$OUT/${tag}.flagstat.txt" | head -1)
  printf '%s\t%s\t%s\n' "$tag" "${mp:-NA}" "${pm:-NA}" >> "$OUT/reference_screening.tsv"
}

# 已有的 K-12 参考（带索引，直接复用）
map_one "K-12_MG1655(existing)" "ref/${REF_NAME}.fa"

# 下载到的各菌株
echo "$REFS" | while IFS='|' read -r short acc asmname; do
  [ -n "${short:-}" ] || continue
  [ "$short" = "K-12_MG1655" ] && continue      # 上面已单独比对
  map_one "$short" "$REFDIR/${short}.fa"
done

# 合并参考
map_one "ALL_COMBINED" "$combined"

# ---------- 6.6 报告 ----------
echo
echo "================ 参考基因组筛选结果 ================"
column -t -s $'\t' "$OUT/reference_screening.tsv" 2>/dev/null || cat "$OUT/reference_screening.tsv"
echo
echo "判读："
echo "  · 某菌株明显高于 78.57%   → 参考株不匹配，换该参考重跑 03_align.sh"
echo "  · ALL_COMBINED 接近 99%   → 未比对 reads 属于菌株特有序列，属生物学差异"
echo "  · ALL_COMBINED 仍停在 ~84% → 未比对部分不来自 E. coli，提示样本含其他物种"
echo "    （本样本已证实混入约 16–20% 的 K. pneumoniae，因此上限约 84%，到不了 99%）"
echo "  · 全部停在 78% 左右且彼此无差异 → 数据本身混合/受污染，需换数据集"
echo
echo "明细日志：logs/refcheck_*.bwa.log   结果表：$OUT/reference_screening.tsv"
