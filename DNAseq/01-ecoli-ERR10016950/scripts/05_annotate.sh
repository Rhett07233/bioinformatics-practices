#!/usr/bin/env bash
# 步骤 5（可选）：变异注释（SnpEff）
#
# 这一步需要下载注释数据库（物种相关，几百 MB 到 1 GB 不等），是可选步骤。
# 先查看可用库名：
#   snpEff databases | grep -i "Escherichia"    # 大肠杆菌
# 然后把库名填进 config 的 SNPEFF_DB。
source "$(dirname "$0")/lib.sh"

if [ -z "${SNPEFF_DB:-}" ]; then
  cat <<'EOF'
配置里没有设置 SNPEFF_DB，跳过注释步骤。

启用方法：
  1) snpEff databases | grep -i Escherichia   # 找一个可用库名
  2) 在该 config 文件里加一行： SNPEFF_DB=<库名>
  3) 重跑： CONFIG=config/ecoli.env bash scripts/05_annotate.sh
EOF
  exit 0
fi

echo "== 检查/下载 SnpEff 数据库：$SNPEFF_DB =="
snpEff download "$SNPEFF_DB" 2> logs/snpeff-download.log || \
  echo "(数据库已存在或下载失败，继续尝试使用本地库)"

for vcf in vcf/*.filtered.vcf.gz; do
  [ -s "$vcf" ] || continue
  run="$(basename "$vcf" .filtered.vcf.gz)"
  echo "== snpEff: $run =="
  snpEff -nodownload -noStats "$SNPEFF_DB" "$vcf" > "vcf/${run}.ann.vcf" 2> "logs/${run}.snpeff.log"

  echo "--- $run 变异后果类型分布（前 15）---"
  grep -v '^#' "vcf/${run}.ann.vcf" \
    | grep -o 'ANN=[^;]*' \
    | cut -d'|' -f2 | tr '&' '\n' \
    | sort | uniq -c | sort -rn | head -15
done

echo
echo "== 注释完成，结果见 vcf/*.ann.vcf =="
