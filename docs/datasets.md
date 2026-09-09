# 数据集说明

本发布包 **不附带** 原始测序或大型 `.h5ad`。以下记录本项目实际用过的主要数据，便于复现与引用。

工程内处理后的 ATAC h5ad（含 `varm['indices']`）通常位于完整仓库：

`12.atac_impute/1.source_from_zyx/1.datasets/output/`

预处理流水线：原始/公开矩阵 → 过滤 → peak 映射 ChromBERT bin → `output/*.h5ad`。

---

## 主基准（细胞状态 / Fig.2）

| 名称 | 物种/build | 约规模 (cells × peaks) | 来源备注（工程内） | 输出文件名 |
|------|------------|------------------------|--------------------|------------|
| 10XPBMC | human / hg38 | 10032 × 106286 | 10x Multiome PBMC；工程标注 from scOpen 流程 | `10XPBMC.h5ad` |
| Hemato | human / hg38 | 2210 × 106900 | from scOpen | `Hemato.h5ad` |
| Tcells | human / hg38 | ~765 × 47569 | from scOpen | `Tcells.h5ad` |
| cell_line | human / hg38 | 1224 × 121541 | from scOpen | `cell_line.h5ad` |
| BMMC | human / hg38 | 大（多 site） | GEO **GSE194122** | `GSE194122_BMMC*.h5ad` |

请在论文 Methods 中补全各自原始论文 / GEO / 10x 官方下载链接与许可证。

---

## 深度下采样（Fig.3）

- 基于 10x PBMC ATAC fragments（granulocyte sorted 10k 一类文库）
- 细胞集合固定约 **9631**（与正式 `10XPBMC.h5ad` 10032 细胞集合不同，勿混用）
- 深度：100 / 75 / 50 / 25 / 10%（另有更低深度探索）
- 每深度独立 Binomial 下采样 fragments multiplicity，并 **重新 call peaks**

工程路径：`12.atac_impute/2.plot/4.read_downsample/`

---

## 区域恢复 / 外推（Fig.4）

- PBMC 低深度矩阵 + ENCODE 风格 cCRE（工程内 `ccre_process.bed`）
- Strategy3：与低深度 peak 无 overlap 的随机 cCRE（约 50k）
- Candidate ±5 kb：邻域 cCRE 集合（约 35 万）+ ChromBERT extrapolate

---

## 跨模态（Fig.5）

- 小鼠胚胎 **E7.25** multiome（工程命名 `E7_25` / `2.e725`）
- build：**mm10**
- 配对 barcode 约 7914
- **不要**与 NMP 模块的 E8.x 数据混称为同一阶段

工程路径：`12.atac_impute/4.cross_modal/`

---

## NMP 应用（Fig.6）

- 小鼠胚胎对象名 `E7_5_to_E8_75`，分析过滤为 **E8.0 / E8.5 / E8.75**
- build：mm10
- NMP 子集约 2458 cells

工程路径：`12.atac_impute/5.driver_tf/`

---

## 上传 GitHub 时建议放什么

| 放入仓库 | 不要放入 |
|----------|----------|
| 本包代码与 docs | `*.h5ad` / fragments / BAM / FASTQ |
| 小型 BED 定义（若 < 数 MB 且许可允许） | 完整 ChromBERT zarr |
| `configs/` 与示例命令 | 私有绝对路径日志 |

可在 `data/README.md` 写下载脚本指针，用 Git LFS 仅当必要时托管中等大小的示例矩阵。
