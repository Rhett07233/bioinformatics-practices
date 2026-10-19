# 生物信息学练习集 · Bioinformatics Practices

*End-to-end bioinformatics pipelines built from public raw data — fully scripted, reproducible, and reconciled against reference results.*

这是我在学习生物信息学分析过程中完成的两个端到端练习。每个练习都从**公开原始数据**出发，
我自己写脚本走完完整分析链，并对结果做了逐项核对——不只是"跑出图了"。
全部脚本都能在笔记本上重跑，不依赖服务器。

**技术栈**：BWA · samtools · bcftools · SnpEff · Seurat (R) · Scanpy (Python) · conda · bash / WSL

## English Summary

Two self-contained exercises, each a numbered script chain plus a one-command `run_all.sh`:

- **scRNA-seq** (`scRNAseq/01-pbmc3k/`) — 10x PBMC 3k implemented **twice and independently**
  in Seurat (R) and Scanpy (Python), including a cross-implementation reconciliation step:
  2,700 → 2,638 cells, 13,714 genes, 9 clusters, identical median QC metrics on both sides.
- **DNA-seq variant calling** (`DNAseq/01-ecoli-ERR10016950/`) —
  FASTQ → fastp/FastQC/MultiQC → BWA-MEM → samtools fixmate/markdup → bcftools → filtered VCF.
  Real *E. coli* WGS, ~245×, 20–30 min on a laptop. Ships a documented debugging log
  (a ploidy-parameter bug and a tool-version breakage, both fixed).

Design: configuration-driven, scripts locate themselves (`cd "$(dirname "$0")/.."`),
raw data and intermediate files excluded via `.gitignore` and re-downloaded on demand.

## 结果速览

**单细胞（S01）** — 细胞类型注释与 QC 指标分布：

| 细胞类型注释 | QC 三大核心指标 |
|---|---|
| ![细胞类型注释 UMAP](scRNAseq/01-pbmc3k/results/R/figures/10_umap_annotated.png) | ![QC 小提琴图](scRNAseq/01-pbmc3k/results/R/figures/01_qc_violin_before.png) |

**DNA 重测序（D01）** — reads 过滤、比对、变异检出、碱基质量四联图：

![DNA-seq QC 与变异检测汇总](DNAseq/01-ecoli-ERR10016950/qc/summary/ERR10016950.summary.png)

## 练习索引

| 编号 | 练习 | 方向 | 工具链 | 数据 | 耗时 |
|---|---|---|---|---|---|
| S01 | [pbmc3k](scRNAseq/01-pbmc3k/README.md) | 单细胞转录组 | Seurat (R) · Scanpy (Python) | 10x PBMC 3k，2700 细胞 | 5–15 分钟 |
| D01 | [ecoli-ERR10016950](DNAseq/01-ecoli-ERR10016950/README.md) | 重测序变异检测 | BWA · samtools · bcftools · SnpEff | EAEC 临床分离株 WGS，298 万 pairs，245× | 20–30 分钟 |

> 另一个 bulk RNA-seq 项目（表达定量）暂未纳入本仓库。

### 各练习的定位与验证强度

| 编号 | 定位 | 能说明什么 | 边界 |
|---|---|---|---|
| **S01** | **可量化验证** | 与 Seurat 官方教程逐项对账（2700 → 2638 细胞、13714 基因、9 个 cluster），并有 R/Python 双实现交叉验证 | 仅覆盖 PBMC 3k 这一套数据，未涉及批次校正与双细胞检测 |
| **D01** | **流程连通性验证 + 基因内容分析** | 从原始 reads 完整跑通 6 步，定位并修复两个真实 bug；组装后做了毒力/耐药基因注释，并查出样本中混入 *K. pneumoniae* | 该样本 `strain` 元数据为空、无配套完整基因组、**无基准真值**，不宣称变异检出准确率。归因见其 README 第 8–9 节 |

## 目录约定

```
<方向>/<两位编号>-<数据集ID>/
```

- **方向** —— `DNAseq` / `scRNAseq`，以后可扩展 `RNAseq` 等。根目录只保留这几个分类，练到几十份也不会乱。
- **编号** —— 两位补零（`01` 而不是 `1`），按加入顺序递增；排序稳定，引用时可简写为 D01 / S01。
- **数据集 ID** —— 优先用官方 accession（`ERR10016950`、`pbmc3k`、`HG002`），比描述性短名更不容易混淆。
- 同一数据集换方法做多次时，在末尾补方法标签：`01-ecoli-ERR10016950-minimap2`。

命名约束：**不用空格、中文、括号、`#`**（括号是 bash 元字符，`cd DNAseq(1)` 会直接语法错误）。

每个练习目录独立自洽、互不依赖，可以单独拷走使用：

```
<方向>/<编号>-<数据集>/
├── README.md          # 数据来源、步骤、期望结果、能力边界
├── envs/              # conda 环境定义
├── scripts/           # 编号脚本 + run_all.sh
└── config/ 或 data/   # (DNAseq) 运行配置 / (scRNAseq) 下载的矩阵
```

所有脚本都用 `cd "$(dirname "$0")/.."` 自定位，**移动或改目录名不影响运行**。

原始数据、参考基因组和中间文件（`.bam` / `fastq` / `.h5ad` / `.rds`）**不入库**，
规则见各练习目录下的 `.gitignore`。克隆后按各自 README 重新下载即可复现。

## 共同的设计原则

1. **不跳步**：DNA 侧从测序 reads 开始（比对、去重、call 变异），单细胞侧从 count matrix 开始（CellRanger 的比对环节由 10x 官方完成，这也是真实工作流的起点）。
2. **脚本可分段执行**：编号脚本 `01_…` → `06_…`，既能单步跑、也能 `bash scripts/run_all.sh` 一键跑完。
3. **结果要和可核对的东西比对**，而不是"看起来对"：
   - `S01 pbmc3k` 与 Seurat 官方教程逐项比对（2700 → 2638 细胞、13714 基因、9 个 cluster），并与 Scanpy 版交叉验证。
   - `D01 ecoli` 做到"实测值与文档记录一致 + 异常有归因"；因无基准真值，不宣称变异检出准确率。见上方定位表。
4. **双实现交叉验证**（`S01`）：同一份数据用 R/Seurat 和 Python/Scanpy 各跑一遍，两份独立实现互相核对——一个逻辑错误很难在两种实现里同时出现。
5. **边界写清楚**：每个练习的 README 都有「已知限制」一节，说明这套流程没做什么、哪些结论不能下。

## 快速开始

工具链在 Linux 下运行（WSL 的 Ubuntu 即可）：

```bash
cd scRNAseq/01-pbmc3k
conda env create -f envs/scanpy.yaml               # 一次性
bash scripts/run_all_py.sh                         # 5–15 分钟

cd ../../DNAseq/01-ecoli-ERR10016950
CONFIG=config/ecoli.env bash scripts/run_all.sh    # 20–30 分钟
```

各练习的详细依赖、参数说明和结果解读，见其自身 README。

## License

MIT — 见 [LICENSE](LICENSE)。
