#!/usr/bin/env bash
# 步骤 4：变异检测与过滤（bcftools）
#
#   mpileup  → 逐位点 pileup，输出基因型似然
#   call -mv → 多项模型调用变异位点
#   norm     → 左对齐 + 拆分多等位基因（不同工具对 indel 的表示不同，规范化后才能比较）
#   filter   → 按 QUAL 和深度过滤
source "$(dirname "$0")/lib.sh"

# 外显子/目标区域限制（可选）。不设置就在所有有 reads 的位置调用。
R_OPT=""
if [ -n "${BED:-}" ]; then
  [ -s "$BED" ] || { echo "BED 文件不存在：$BED" >&2; exit 1; }
  R_OPT="-R $BED"
  echo ">>> 限定目标区域：$BED"
else
  echo ">>> 未设置 BED，将在全基因组范围调用（targeted 数据会慢一些）"
fi

tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r run name fq1 fq2; do
  [ -n "${run:-}" ] || continue
  [ -s "bam/${run}.markdup.bam" ] || continue
  echo "== bcftools: $run（ploidy=$PLOIDY）=="

  # shellcheck disable=SC2086
  bcftools mpileup -f "$REF" $R_OPT \
      -q 20 -Q 20 -d 250 -a AD,DP -Ou \
      "bam/${run}.markdup.bam" 2> "logs/${run}.mpileup.log" \
    | bcftools call -mv --ploidy "$PLOIDY" -Oz -o "vcf/${run}.raw.vcf.gz"
  bcftools index -t "vcf/${run}.raw.vcf.gz"

  bcftools norm -f "$REF" -m -both -Oz \
      -o "vcf/${run}.norm.vcf.gz" "vcf/${run}.raw.vcf.gz" 2> "logs/${run}.norm.log"
  bcftools index -t "vcf/${run}.norm.vcf.gz"

  bcftools filter -i 'QUAL>=20 && INFO/DP>=10' -Oz \
      -o "vcf/${run}.filtered.vcf.gz" "vcf/${run}.norm.vcf.gz"
  bcftools index -t "vcf/${run}.filtered.vcf.gz"

  bcftools stats "vcf/${run}.filtered.vcf.gz" > "vcf/${run}.stats.txt"
done

echo
echo "== 变异数量（SN=SNP 数，indels）=="
grep -H "^SN" vcf/*.stats.txt | grep -E "number of (SNPs|indels)" || true
echo
echo "== 转换/颠换比 Ti/Tv（仅在与参考高度一致的数据中有判读意义；本项目见 README 8.5）=="
grep -H "^TSTV" vcf/*.stats.txt || true
