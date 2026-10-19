#!/usr/bin/env bash
# 步骤 8（附加分析）：基因内容注释 —— 毒力 / 耐药 / 质粒
#
# 为什么不能靠步骤 4 的 SNP 结果回答"这个菌为什么致病"：
#   本样本是尼日利亚来源的 EAEC（肠聚集性大肠杆菌），致病性来自
#   accessory genome —— pAA 质粒、AAF 黏附素、EAST-1 毒素等。
#   这些序列在 K-12 MG1655 参考里根本不存在，比对阶段那 21% 未比对的 reads
#   正是它们。所以 SNP 集里看不到致病信息，必须走
#   de novo 组装 → 基因注释 这条路。
#
# 产出：
#   assembly/spades/contigs.fasta        组装结果
#   results/gene_content/{vfdb,card,resfinder,plasmidfinder}.tsv
#   results/gene_content/summary.tsv     过滤去重后的基因清单
#
# 用法：
#   CONFIG=config/ecoli.env bash scripts/08_gene_content.sh
#   SUB_PAIRS=1000000 SPADES_MEM=8 CONFIG=config/ecoli.env bash scripts/08_gene_content.sh
source "$(dirname "$0")/lib.sh"

ASM_DIR="assembly"
OUT="results/gene_content"
SUB_PAIRS="${SUB_PAIRS:-800000}"     # 组装 50x 左右足够，且快很多
SPADES_MEM="${SPADES_MEM:-8}"        # GB
mkdir -p "$ASM_DIR" "$OUT" logs

# ---------- 8.0 工具检查 ----------
missing=""
for t in spades.py abricate seqtk; do
  command -v "$t" >/dev/null 2>&1 || missing="$missing $t"
done
if [ -n "$missing" ]; then
  echo "缺少工具:$missing" >&2
  echo "安装：conda env update -f envs/dna.yaml" >&2
  exit 1
fi

run="$(awk 'NR==2{print $1; exit}' "$SAMPLES")"
[ -s "clean/${run}_R1.clean.fq.gz" ] || {
  echo "找不到 clean/${run}_R1.clean.fq.gz —— 请先跑完 02_qc.sh" >&2; exit 1; }

# ---------- 8.1 抽样子集 ----------
# 组装用不到 245x；降到约 50x 可把 SPAdes 耗时压到几分钟。
sub1="logs/${run}.sub_R1.fq.gz"
sub2="logs/${run}.sub_R2.fq.gz"
if [ ! -s "$sub1" ] || [ ! -s "$sub2" ]; then
  echo "== 抽样 $SUB_PAIRS 对 reads（约 50x）=="
  seqtk sample -s 42 "clean/${run}_R1.clean.fq.gz" "$SUB_PAIRS" | gzip > "$sub1"
  seqtk sample -s 42 "clean/${run}_R2.clean.fq.gz" "$SUB_PAIRS" | gzip > "$sub2"
fi

# ---------- 8.2 de novo 组装 ----------
contigs="$ASM_DIR/spades/contigs.fasta"
if [ ! -s "$contigs" ]; then
  echo "== SPAdes 组装（--isolate，$THREADS 线程，${SPADES_MEM} GB）=="
  rm -rf "$ASM_DIR/spades"
  spades.py --isolate \
      -1 "$sub1" -2 "$sub2" \
      -o "$ASM_DIR/spades" \
      -t "$THREADS" -m "$SPADES_MEM" \
      > "logs/spades.log" 2>&1 || {
        echo "SPAdes 失败，见 logs/spades.log" >&2; exit 1; }
fi

n_contig=$(grep -c '^>' "$contigs")
asm_len=$(awk '/^>/{next}{n+=length($0)}END{print n}' "$contigs")
echo ">>> 组装完成：$n_contig 条 contig，总长 $((asm_len/1000000)) Mb"

# ---------- 8.3 按覆盖度去除污染物种 ----------
# 本数据集的组装呈双峰：
#   覆盖度 2–4×、GC 57.3%、约 5 Mb  —— 另一个物种（疑似 Klebsiella 类）
#   覆盖度 8–40×、GC 50–52%、约 5.8 Mb —— E. coli 主基因组 + accessory
# 两峰在 4× 处干净分开，所以按 MIN_COV 过滤，只保留 E. coli 部分再注释。
MIN_COV="${MIN_COV:-4}"
filtered="$ASM_DIR/spades/contigs.filtered.fasta"
if [ ! -s "$filtered" ]; then
  awk -v minc="$MIN_COV" '
    /^>/ {
      name = $1; sub(/^>/, "", name)
      c = 0
      if (match($0, /cov_[0-9.]+/)) c = substr($0, RSTART+4, RLENGTH-4)
      if (c + 0 >= minc) print name
    }' "$contigs" > "$ASM_DIR/spades/keep.list"
  seqtk subseq "$contigs" "$ASM_DIR/spades/keep.list" > "$filtered"
fi
kept=$(grep -c '^>' "$filtered")
kept_len=$(awk '!/^>/{n+=length($0)}END{print n}' "$filtered")
echo ">>> 覆盖度过滤（>=${MIN_COV}x）：保留 $kept / $n_contig 条，$((kept_len/1000000)) Mb"
echo "    （被过滤掉的低覆盖度 contig 属另一物种，不是 E. coli 的 accessory 区）"

# ---------- 8.4 abricate 数据库 ----------
if ! abricate --list 2>/dev/null | tail -n +2 | grep -q .; then
  echo "== 初始化 abricate 数据库（一次性）=="
  abricate --setupdb > "logs/abricate_setupdb.log" 2>&1 || {
    echo "数据库初始化失败，见 logs/abricate_setupdb.log" >&2; exit 1; }
fi

# ---------- 8.5 逐库扫描 ----------
# 注意：只扫描过滤后的 E. coli 部分，避免把另一物种的基因算进来。
# vfdb        毒力因子
# card        耐药基因（综合库）
# resfinder   耐药基因（水平转移相关，适合临床）
# plasmidfinder 质粒复制子（EAEC 的 pAA 质粒常在此出现）
for db in vfdb card resfinder plasmidfinder; do
  echo "== abricate --db $db =="
  abricate --db "$db" "$filtered" > "$OUT/${db}.tsv" 2> "logs/abricate_${db}.log" || true
done

# ---------- 8.5 汇总 ----------
# abricate 列：#FILE SEQUENCE START END STRAND GENE COVERAGE COVERAGE_MAP GAPS %COVERAGE %IDENTITY DATABASE ACCESSION PRODUCT RESISTANCE
# 过滤标准：覆盖度 ≥80% 且一致性 ≥90%，同一基因只保留最佳一条
summary="$OUT/summary.tsv"
printf 'database\tgene\tidentity\tcoverage\tproduct\n' > "$summary"
for db in vfdb card resfinder plasmidfinder; do
  f="$OUT/${db}.tsv"
  [ -s "$f" ] || continue
  # 列位置随 abricate 版本变过（1.4 的输出只有 12 列，没有 STRAND / PRODUCT），
  # 所以按表头名字定位列，而不是写死列号。
  awk -F'\t' -v db="$db" '
    NR==1 {
      gi=0; ci=0; ii=0; pi=0
      for (i=1; i<=NF; i++) {
        if ($i == "GENE")      gi = i
        if ($i == "%COVERAGE") ci = i
        if ($i == "%IDENTITY") ii = i
        if ($i == "PRODUCT")   pi = i
      }
      if (!gi) gi = 5
      if (!ci) ci = 9
      if (!ii) ii = 10
      next
    }
    {
      cov = $(ci) + 0
      idt = $(ii) + 0
      if (cov >= 80 && idt >= 90) {
        k = db SUBSEP $(gi)
        if (!(k in best) || idt > iden[k]) {
          best[k] = 1; iden[k] = idt
          prod = (pi ? $(pi) : "")
          line[k] = sprintf("%s\t%s\t%.1f\t%.1f\t%s", db, $(gi), idt, cov, prod)
        }
      }
    }
    END { for (k in line) print line[k] }
  ' "$f" >> "$summary"
done

# 排序（表头保留）
{ head -1 "$summary"; tail -n +2 "$summary" | sort -t$'\t' -k1,1 -k2,2; } > "$summary.tmp" && mv "$summary.tmp" "$summary"

echo
echo "================ 基因内容注释结果 ================"
for db in vfdb card resfinder plasmidfinder; do
  n=$(awk -F'\t' -v d="$db" 'NR>1 && $1==d' "$summary" | wc -l)
  printf '  %-15s %s 个基因\n' "$db" "$n"
done
echo
echo "== 毒力基因（vfdb）=="
awk -F'\t' 'NR>1 && $1=="vfdb" {printf "  %-14s id=%.0f%% cov=%.0f%%  %s\n", $2,$3,$4,$5}' "$summary" | head -25
echo
echo "== 耐药基因（card / resfinder）=="
awk -F'\t' 'NR>1 && ($1=="card"||$1=="resfinder") {printf "  [%s] %-16s id=%.0f%%  %s\n", $1,$2,$3,$5}' "$summary" | head -20
echo
echo "== 质粒复制子（plasmidfinder）=="
awk -F'\t' 'NR>1 && $1=="plasmidfinder" {printf "  %-14s %s\n", $2,$5}' "$summary" | head -10
echo
echo "完整表：$summary"
echo "原始命中：$OUT/{vfdb,card,resfinder,plasmidfinder}.tsv"
