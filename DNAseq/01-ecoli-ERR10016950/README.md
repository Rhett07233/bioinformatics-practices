# DNA 重测序变异检测流程

我从 ENA 下载原始测序数据，独立走完
**质控 → 比对 → 排序去重 → 变异检测 → 过滤 → 注释** 全流程。

和另一个 RNA-seq 项目（暂未纳入本仓库）互补：
那条线练的是表达量分析，这条线练的是 **BWA + samtools + bcftools** 这套 DNA 变异检测工具链。

---

## 1. 数据（已核实）

### 数据：大肠杆菌 WGS

| 项目 | 值 |
|---|---|
| ENA run | **ERR10016950** |
| BioSample | `SAMEA9333667`（二级号 `ERS7066877`） |
| 所属 study | **PRJEB8667** — *Pathogenic lineages of enteric bacteria in Nigeria*（Sanger Institute） |
| 物种 / 类型 | *Escherichia coli* / WGS，双端 |
| 测序平台 | Illumina HiSeq X Ten，151 bp |
| 数据量 | ENA 元数据记 3.77M read pairs；**实测文件为 298 万 pairs**（见第 8 节），约 270 MB |
| 覆盖度 | 约 **245×**（按 4.64 Mb 基因组推算） |
| 参考基因组 | E. coli K-12 MG1655，**GCF_000005845.2**（1.3 MB） |
| 端到端耗时 | **20–30 分钟** |

> ⚠️ **数据集局限（重要）**
>
> 该 BioSample 的 `strain`、`serotype`、`sample_description` **全部为空**，
> 也没有关联的已组装基因组。结合 study 描述可判断它属于
> **肠聚集性大肠杆菌（EAEC）**，但**无法确定到具体菌株**。
>
> 因此本项目的定位是「**验证流程能否跑通**」，**不是**「变异检测基准」：
> 它没有 truth set，无法判断检出的 21,137 个 SNP 里有多少是正确的。
> 若需要可量化的验证，应换用「元数据完整 + 有配套完整基因组」的数据集。
>
> 详细核查过程见 [第 8.6 节](#86-数据溯源核查)。

相关页面：
[ENA run](https://www.ebi.ac.uk/ena/browser/view/ERR10016950) ·
[study](https://www.ebi.ac.uk/ena/browser/view/PRJEB8667) ·
[sample](https://www.ebi.ac.uk/ena/browser/view/SAMEA9333667)

```bash
CONFIG=config/ecoli.env bash scripts/run_all.sh
```

下载链接来自 ENA 官方 API 返回的 `fastq_ftp` 字段，已核对为双端配对的两个文件。

---

## 2. 参考基因组（已查证）

| 文件 | 大小 |
|---|---|
| `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz` | 1.3 MB |

参考基因组取自 NCBI RefSeq（E. coli K-12 MG1655，GCF_000005845.2）。
路径写在 `config/ecoli.env` 里，脚本检测到文件不存在会自动下载。

---

## 3. 环境

```bash
wsl -d Ubuntu
cd /path/to/bioinformatics-practices/DNAseq/01-ecoli-ERR10016950     # 换成自己的实际路径
conda env create -f envs/dna.yaml
conda activate dna
```

**资源预算**（大肠杆菌 WGS）

| 项目 | 需求 |
|---|---|
| 磁盘 | 约 1 GB（fastq 0.25 GB + BAM 0.27 GB + 参考与索引 < 0.1 GB） |
| 内存 | 2 GB 足够（bwa index 与比对） |
| 时间 | 全流程 20–30 分钟（下载数分钟、比对约 10 分钟、变异检测数分钟） |

---

## 4. 流程

| 步骤 | 工具 | 产出 |
|---|---|---|
| 1 下载 | aria2 / wget | `fastq/*.fastq.gz` |
| 2 质控 | fastp + FastQC + MultiQC | `qc/multiqc/multiqc_report.html` |
| 3 比对去重 | **bwa mem** → samtools | `bam/<run>.markdup.bam` |
| 4 变异检测 | bcftools mpileup/call/norm/filter | `vcf/<run>.filtered.vcf.gz` |
| 5 注释（可选） | SnpEff | `vcf/<run>.ann.vcf` |
| 6 参考筛选（诊断，可选） | bwa mem + samtools flagstat | `qc/refcheck/reference_screening.tsv` |
| 7 汇总图（可选） | matplotlib | `qc/summary/<run>.summary.png` |

### 为什么顺序必须是「名字排序 → fixmate → 坐标排序 → markdup」

`samtools markdup` 需要 mate 的坐标信息才能判断一对 reads 是否重复，
所以流程是：

```
bwa mem → samtools sort -n（按名字排序，让 mate 相邻）
        → samtools fixmate -m（补齐 mate 标签、MC/ms 标记）
        → samtools sort（按坐标排序）
        → samtools markdup（标记重复）
```

漏掉 `fixmate` 直接 markdup，samtools 会直接报错。

### 三个关键概念

**① 为什么 DNA 要去重，而 RNA-seq 不用？**
PCR 扩增产生的重复片段会被误判成「支持变异的证据」，造成假阳性，所以 DNA 必须去重。
RNA-seq 的重复主要来自高表达基因，是真实信号，去重反而有害。

**② 为什么 `bcftools norm -f ref -m -both`？**
不同工具对 indel 的表示方式不同（同一位点可能写成 `AT>A` 或 `A>-`），
`norm` 做**左对齐（left-align）**并拆分多等位基因，保证结果可与其他工具/数据库比较。

**③ Ti/Tv 比值是什么？**
转换（A↔G、C↔T）与颠换的比值。在样本与参考高度一致的数据里，真实变异通常在 **2.0 以上**；
只有 1.2 左右往往意味着假阳性偏多，因此常被当作质控参考。`04_call.sh` 会自动打印。

> ⚠️ **适用边界**：Ti/Tv 作为 QC 指标只在「样本与参考足够接近」时成立。
> 本项目是"致病分离株 vs 实验室 K-12"，SNP 集本质是谱系分歧，Ti/Tv 反映突变谱而非错误率，
> **不能据此判断流程好坏**。详见第 8.5 节。

---

## 5. 跑完怎么自检

| 检查项 | 参考值 / 判读 |
|---|---|
| 比对率（flagstat `mapped`） | 同种参考通常 >95%；本数据集实测 **78.57%**，原因见第 8.5 节 |
| 重复率（markdup 报告） | 一般 5%–20%；本数据集实测 **5.7%** |
| 变异数量 | 与 MG1655 参考不同的分离株会有几千到几万个变异，属正常；本数据集实测 **21,137 SNP + 251 indel** |
| Ti/Tv | 细菌无公认基准，**不建议作为判断依据**（见第 8.5 节） |

`04_call.sh` 跑完会把 SNP 数、indel 数、Ti/Tv 直接打在屏幕上；本项目把这些数字连同运行环境与参数一并记录，便于复现和横向对比。

---

## 6. 输出与目录

```
01-ecoli-ERR10016950/
├── README.md
├── config/          ecoli.env   ← 换数据集只改这里
├── samples/         样本表（含 ENA 直链）
├── envs/dna.yaml
├── scripts/         lib.sh + 01~07 + run_all.sh
├── fastq/ clean/ bam/ ref/     中间文件（.gitignore 已排除）
├── qc/              MultiQC 报告 + summary/（汇总图）+ refcheck/（参考筛选结果）
└── vcf/             raw / norm / filtered VCF + stats
```

> **生成汇总图**（只读已有结果，不重跑流程）：
>
> ```bash
> conda activate scanpy          # 该环境自带 matplotlib
> python scripts/07_summary_plot.py
> ```
>
> 输出 `qc/summary/<run>.summary.png`：reads 过滤前后、比对率、变异数、Q30 四联图。

### 结果概览

![QC 与变异检测汇总](qc/summary/ERR10016950.summary.png)

四格依次为：reads 过滤前后（3.77M → 3.52M）、比对率（78.57%，其中 properly paired 76.25%）、
变异数（21,137 SNP / 251 indel，对数轴）、Q30 碱基比例（85.6% → 88.2%）。
图中"21.43% 未比对"的成因见 [第 8.5 节](#85-两个异常的归因2026-10-修订其一已定位)。

入库内容：`scripts/`、`config/`、`samples/`、`envs/`、`README.md`、
`vcf/*.vcf.gz`（约 2 MB）、`vcf/*.stats.txt`、`qc/` 下的报告与 TSV。
BAM、fastq、参考基因组及其索引、`.rds`/`.h5ad` 都在 `.gitignore` 里。

---

## 7. 这套流程的边界与已知取舍

- **单样本**：真实项目里多样本会做 joint calling（GATK GenomicsDBImport 或 `bcftools merge`），
  而不是各自调用再合并。
- **没有 BED 时是全基因组调用**：targeted（外显子组 / panel）数据因为只有目标区有 reads，
  结果基本等价，但会慢不少。生产流程一定用试剂盒对应的 target BED。
- **没有做 BQSR / 变异质控建模**：GATK 最佳实践里有，bcftools 流程通常省略。
- **没做 CNV / SV**：那是另一套工具（CNVkit、Manta、Delly）。

---

## 8. 已知限制与排查记录（2026-10 大肠杆菌自测）

用 ERR10016950 完整跑通全流程后的**实测结果**，以及排查过程中发现的两个 bug。
数字全部来自实际运行，限制如实写出。

### 8.1 实测结果

| 指标 | 实测值 |
|---|---|
| 原始 reads | 5,969,474（约 298 万 read pairs） |
| fastp 后通过过滤 | 3,524,108 条（**41% 因低质量被丢弃**，Q40 bases = 57.3%） |
| 比对率（flagstat `mapped`） | **78.57%** |
| 重复率（flagstat `duplicates`） | 202,648 / 3,524,108 ≈ **5.7%** |
| 最终 SNP | **21,137** |
| 最终 indel | **251** |
| Ti/Tv | **2.51** |

> ENA 元数据里写的 "3.77M read pairs" 与实际文件（298 万 pairs）不符。
> **文档中一律引用实测值**，不要直接抄元数据。

### 8.2 已修复的 bug ①：`bcftools call` 缺少 `--ploidy 1`

初版命令是 `bcftools call -mv`，即使用**默认的二倍体模型**。大肠杆菌是单倍体，后果很大：

| | 修复前（ploidy=2） | 修复后（ploidy=1） |
|---|---|---|
| SNP | 84,425 | **21,137** |
| indel | 338 | 251 |
| 基因型 | 出现 `0/1` 杂合 | 全部 `1/1` |

**84,425 个 SNP 里有约 6.3 万个是二倍体模型造出来的假阳性**——单倍体基因组不可能有杂合位点。
`bcftools call --help` 里 `--ploidy` 的默认值明确写着 `[2]`，所以这个参数必须显式指定。

现已把倍性做成配置项，避免二倍体物种与细菌（单倍体）互相污染：

- `scripts/lib.sh`：`PLOIDY="${PLOIDY:-2}"`（默认二倍体），并 export
- `config/ecoli.env`：显式设 `PLOIDY=1`
- `scripts/04_call.sh`：`bcftools call -mv --ploidy "$PLOIDY"`

### 8.3 已修复的 bug ②：`fastqc -o` 在新版本上失效

FastQC 0.12+ 的输出目录参数只接受 `--outdir`，`-o` 会直接报
`fastqc: error: unrecognized arguments: -o`。因为脚本开头有 `set -e`，
这个报错会让整个流程停在第 2 步末尾（MultiQC 也一并没跑）。

已把两个项目的 `02_qc.sh` 统一改成 `--outdir`。

> 这类"老参数在新版本静默失效"的问题不容易第一时间发现：`envs/dna.yaml` 没有固定版本，
> 每次 `conda env create` 都可能换上更新的工具。

### 8.4 一个已排除的"疑似问题"：indel 是不是丢了？

21,137 个 SNP 只配 251 个 indel（约 84:1），乍看像 indel 检出失败。
排查后得到一个**有边界的结论**：**调用阶段没有丢 indel，但"样本本身 indel 就少"仍未定论**（见下）。

**直接证据：用 CIGAR 统计（独立于变异调用器）**

```bash
samtools view bam/ERR10016950.markdup.bam | head -200000 | cut -f6 \
  | awk '/[ID]/{n++} END{printf "%d / %d = %.1f%%\n", n, NR, n/NR*100}'
# 输出：2091 / 200000 = 1.0%
```

推算：277 万条 mapped reads × 1.0% ≈ 2.78 万条跨 indel 的比对；按约 90× 覆盖度，
每个 indel 位点被约 90 条 reads 覆盖 → **约 300 个 indel 位点**，与实测 251 个吻合。

**这条证据能证明什么、不能证明什么（重要，别说过头）：**

- ✅ **能证明**：比对器已经表示成 `I`/`D` 的 indel，调用器一个没丢——
  CIGAR 是 bcftools 之前就定下来的，不受任何调用参数影响。
- ❌ **不能证明**：样本本身的 indel 就少。如果 bwa 把 indel 用**错配**或 **soft clip**
  糊过去（短 indel 落在重复序列里时很常见），CIGAR 里就不会出现 `I`/`D`，
  这套统计对此是**盲的**。

所以准确的结论是：**"调用阶段没丢"成立；"比对阶段没藏、样本本身 indel 就少"仍未定论。**

**若要继续区分**，可看 soft clip 比例（若异常高，说明比对器把序列切掉而不是用 gap 表示）：

```bash
samtools view bam/ERR10016950.markdup.bam | head -200000 | cut -f6 \
  | awk '/S/{s++} /[ID]/{n++} END{printf "含S: %d (%.1f%%)  含I/D: %d (%.1f%%)\n", s, s/NR*100, n, n/NR*100}'
```

**同时排除了另外三个假设**：

| 假设 | 验证方法 | 结果 |
|---|---|---|
| 被 `bcftools filter` 砍掉 | 对比过滤前后 | 过滤前 345 → 过滤后 338，**不是 filter 造成的** |
| 重复序列区误比对 | 按 MQ 阈值过滤 | MQ≥60 时 21,137 个 SNP 仍剩 16,105 个，**不是假象** |
| 先验突变率过低 | `-P` 从 1.1e-3 提到 5e-2 | indel 始终 20–21 个，**纹丝不动** |

> 复盘（两处判断失误）：
> 1. 最初判断"84:1 异常"，依据是"细菌 SNP:indel 约 5–15:1"这个**未经核实的经验比例**——
>    这个判断本身不严谨。观测比例 = **真实比例 × 检出效率差**，SNP 和 indel 的检出难度
>    根本不同，光看比例判断不了对错。
> 2. 后来用 CIGAR 做直接测量，方向是对的，但**结论一度说强了**：它排除了"调用阶段丢失"，
>    没有排除"比对阶段隐藏"。

### 8.5 两个异常的归因（2026-10 修订：其一已定位）

**① 比对率 78.57% —— 已定位为「参考株不匹配」（强推测，可实测验证）**

最初归因于"数据质量/纯度问题"，**这个归因是错的**。核查 ENA 元数据后：

- 该 run 所属 study 为 **PRJEB8667 — "Pathogenic lineages of enteric bacteria in Nigeria"**，
  研究描述明确写的是 *non-typhoidal Salmonella and **enteroaggregative E. coli*** ——
  样本属于 **EAEC（肠聚集性大肠杆菌）**，是临床/致病分离株，不是实验室驯化的 K-12。
- *E. coli* 泛基因组差异极大（不同菌株间可有 1–2 Mb 特有序列）。EAEC 常携带
  pAA 质粒、AAF 黏附素、EAST-1 毒素等序列，**这些在 K-12 参考里根本不存在，当然比对不上**。
  这属于生物学差异，不是数据缺陷。

验证脚本：`scripts/06_reference_check.sh` —— 用**同一份 reads** 比对 6 个代表性
*E. coli* 完整基因组，并额外比对一个"6 株合并参考"。

判据：若合并参考的比对率接近 99%，即可确认未比对的 reads 属于菌株特有序列。

> 本仓库**未附该验证的结果表**（脚本需从 NCBI 下载约 9 MB 参考基因组，
> 国内网络条件下可能较慢）。上面的结论是**基于元数据的推断**，脚本保留在仓库里，可自行实测。

**② Ti/Tv = 2.51 —— 这个指标在本场景下本来就不该当 QC 用**

Ti/Tv 的常用预期区间（约 2.0–2.1）是针对**"与参考高度一致的样本"**建立的。
本项目的 SNP 集是"致病分离株 vs 实验室 K-12"，本质是**谱系分歧**，反映突变谱而非测序错误率。
所以"细菌应约 1.0"这个预期本来就没有依据（当时也未找到出处）。

**结论修订**：78.57% 不再是"未解之谜"，而是**可解释的参考选择问题**。
通用原则：**遇到数字异常，先查元数据，再怀疑数据质量。**

### 8.6 数据溯源核查

用 ENA 官方 API 拉取样本 XML 后的完整结果（2026-10 核查）：

| 字段 | 值 |
|---|---|
| BioSample | `SAMEA9333667`（二级号 `ERS7066877`） |
| 提交者内部 ID | `3534STDY10261424`（Sanger 测序批次号，**不是菌株名**） |
| 物种 | *Escherichia coli*（taxid 562） |
| `strain` | **空** |
| `serotype` / `isolate` / `sample_description` | **空** |
| 关联的已组装基因组 | **无**（仅提供 FASTQ 链接） |

**含义：该样本无法从公开记录确定到菌株级别。** 若要判断谱系，只能从 reads 计算
（7 基因 Achtman MLST，或 8 基因 Pasteur 方案定 phylogroup）。

**因此本数据集的适用范围是「流程连通性验证」，不是「变异检测基准」**：
没有 truth set，就无法判断 21,137 个 SNP 中有多少是真的。

**结论**：这条流程的逻辑与实现是通的——完整跑过 6 个步骤、产出 VCF、参数化可切换数据集。
数据侧的保留意见只有一条：**没有基准真值**。若要升级为可量化验证，
应换用「元数据完整 + 有配套完整基因组」的样本，这样比对率与变异集才有明确的判读标准。

> **这轮的结论**：选数据前先确认「元数据完整 + 有配套完整基因组」。
> 元数据的 `strain` 字段是否为空，一条 ENA API 请求两分钟就能查到——这次是事后才查的。

---

## 9. 设计取舍与常见问题

### 9.1 工具与参数为什么这么选

**为什么用 BWA-MEM 而不是 BWA-backtrack？**

`backtrack` 是早期为 35–100 bp 短读长设计的算法，不支持 split alignment
（把一条 read 拆成两段分别比对到基因组不同位置，跨结构变异或剪接时必需），
读长一长效率就明显下降。`mem` 现在是事实标准，也是本项目用的那一个——没有理由再用 backtrack。

**`bcftools mpileup` 的 `-q` 和 `-Q` 有什么区别？**

| 参数 | 含义 | 回答的问题 |
|---|---|---|
| `-q` | 最小**比对质量** MAPQ | 这条 read 放对地方了吗 |
| `-Q` | 最小**碱基质量** BASEQ | 这个碱基读准了吗 |

两者作用在不同层面：一条 read 可以比对得很确定（高 `-q`），但其中某个碱基测序质量很低（低 `-Q`）。
本项目都是 `20`，即低于 20 的直接不参与调用。

**为什么要显式写 `-d 250`？**

`-d/--max-depth` 是每个位点参与调用计算的 reads 上限，bcftools 默认值就是 250。
极高覆盖度的位点（重复区、rRNA 区）如果不设上限，内存和耗时会失控。
本项目把它显式写出来，是为了让参数**可见、可调**，而不是藏在默认值里。

### 9.2 结果怎么判断可信

**比对率多少算正常？偏低时按什么顺序排查？**

同种参考通常 >95%（去重前）；本数据集实测 **78.57%**，原因见第 8.5 节。明显偏低时按证据强度依次排查：

1. **接头污染**——看 fastp 报告的接头检测结果
2. **样本污染 / 物种不符**——比对上的 reads 里是否有大量其他物种
3. **参考基因组版本错**——染色体命名（`chr1` vs `1`）、组装版本不匹配会导致大面积比对失败
4. **参考里缺序列**——用了 `primary_assembly` 就不含 alt/scaffold，那是预期行为不是 bug

**怎么判断一个变异是真的？**

证据强度从低到高：

| 层级 | 依据 | 性质 |
|---|---|---|
| 1 | QUAL（位点整体置信度）、DP（覆盖深度，过低不可信） | 生信内部检查 |
| 2 | AD（等位基因深度比，接近 0.5 或 1.0 才像真变异） | 生信内部检查 |
| 3 | 链偏好性——真变异正负链支持应均衡 | 生信内部检查 |
| 4 | 是否在 dbSNP / gnomAD / 1000G 中有记录 | 外部数据库 |
| 5 | 独立方法验证：Sanger 测序 | **实验验证** |

前四层都只是**生信内部的一致性检查**，能提高可信度但都不能证明变异真实存在；
只有第五层是直接证据。生产流程会用正交数据（同一材料的另一种测序）替代 Sanger。

### 9.3 已在本文档其他位置说明的问题

| 问题 | 位置 |
|---|---|
| 为什么顺序必须是「名字排序 → fixmate → 坐标排序 → markdup」 | 第 4 节 |
| 为什么 DNA 要去重，RNA-seq 不用 | 第 4 节 ① |
| `bcftools norm` 在做什么、为什么要左对齐 | 第 4 节 ② |
| Ti/Tv 比值是什么、什么情况下才能用它判断 | 第 4 节 ③、第 8.5 节 |
| 单样本与多样本的差异（joint calling） | 第 7 节 |
| 流程怎么复用、中间文件怎么管理 | 第 6 节 |
| 参数化与倍性配置踩过的坑 | 第 8.2 / 8.3 节 |
