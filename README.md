# ChromBERT-scAI

稀疏 scATAC-seq 可及性矩阵插补：用固定的 ChromBERT 区域表征 \(E\) 约束低秩重建

\[
Y \approx E A B
\]
<img width="1600" height="900" alt="image" src="https://github.com/user-attachments/assets/e92fd101-7e8c-4ccb-b727-94deac5b2c2d" />

本目录是从完整分析工程中整理出的 **GitHub 发布包**，只包含模型代码、关键脚本与数据集说明，不含大型矩阵 / fragments / checkpoint。


---

## 目录结构

```
ChromBERT-scAI_github/
├── README.md
├── requirements.txt
├── .gitignore
├── chrombert_scai/          # 核心模型
│   ├── impute_als.py        # 生产版 ALS（推荐）
│   ├── impute_als_ablation.py  # real / random / shuffled E 消融
│   └── impute_ridge.py      # 早期同 E 的 ridge 对照
├── scripts/
│   ├── run_impute.sh        # 命令行示例
│   ├── prepare_peak_indices.py  # peak → ChromBERT bin index
│   └── evaluate_clustering.py   # ARI/NMI 等聚类评估
├── configs/
│   └── default_impute.yaml
├── docs/
│   ├── model.md             # 模型公式与实现细节
│   ├── datasets.md          # 使用的数据集与获取方式
│   └── embedding.md         # E (zarr) 说明
├── annotations/             # ChromBERT bin BED（hg38 / mm10）
│   ├── hg38_6k_1kb_region.bed
│   └── mm10_5k_1kb_region.bed
└── examples/
    └── minimal_usage.md
```

---

## 快速开始

### 1. 环境

```bash
conda create -n chrombert_scai python=3.10 -y
conda activate chrombert_scai
pip install -r requirements.txt
# 另需 bedtools（peak 索引映射）
```

### 2. 输入要求

输入 `.h5ad`（cell × peak）需包含：

- `adata.X`：稀疏/稠密 peak 计数（或二值）
- `adata.varm['indices']`：每个 peak 对应 ChromBERT 基因组 bin 的整数下标  
  （用 `scripts/prepare_peak_indices.py` 生成）

Bin BED（已随仓库提供）：

- hg38: `annotations/hg38_6k_1kb_region.bed`
- mm10: `annotations/mm10_5k_1kb_region.bed`

另需 ChromBERT 区域 embedding 的 zarr 目录（含 `full` 或 `PCA/<dim>`）。
### 3. 运行插补

```bash
python -u chrombert_scai/impute_als.py \
  --input  data/10XPBMC.h5ad \
  --output output/10XPBMC_imputed_ALS.h5ad \
  --emb_dir /path/to/atac_embs_dim.zarr \
  --emb_dim 768 \
  --topics 128 \
  --lambdaA 1 \
  --lambdaB 0.1 \
  --als_iters 20 \
  --cg_iters 500 \
  --binary \
  --seed 42
```

或编辑后运行：

```bash
bash scripts/run_impute.sh
```

### 4. 输出字段

| 字段 | 含义 |
|------|------|
| `adata.X` | 插补后的相对可及性分数（float16；**可为负**） |
| `adata.obsm['latent']` | 细胞表征 \(B^\top\) |
| `adata.varm['topic']` | 区域表征 \(EA\) |
| `adata.uns['ALS_A']` | 因子 \(A\) |

对未见峰外推：`E_test @ ALS_A @ latent.T`。

---

## 引用与先验来源

- 本方法：ChromBERT-scAI（固定 \(E\) + ALS 低秩重建）
- 区域表征来自 ChromBERT 预训练产物；**请勿在未核对 checkpoint 前称为纯 DNA sequence language model**（已发表 ChromBERT 主要基于调控因子 ChIP-seq 共关联预训练）

---

## 许可

代码按你仓库选定的开源协议发布。数据集请遵循各自原始数据许可与引用要求（见 `docs/datasets.md`）。
