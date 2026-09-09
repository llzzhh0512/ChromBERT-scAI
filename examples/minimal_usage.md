# 最小使用示例

## 前置

1. 安装 `requirements.txt` 与 `bedtools`
2. 准备 ChromBERT embedding zarr（见 `docs/embedding.md`）
3. 准备基因组 bin BED，放入 `annotations/`

## 步骤

```bash
# 1) 给 peak 矩阵写 varm['indices']
# hg38:
BIN_BED=annotations/hg38_6k_1kb_region.bed
# mm10:
# BIN_BED=annotations/mm10_5k_1kb_region.bed

python scripts/prepare_peak_indices.py \
  --input data/raw_peaks.h5ad \
  --bins "$BIN_BED" \
  --output data/peaks_indexed.h5ad \
  --drop-unmapped

# 2) ALS 插补
export EMB_DIR=/path/to/atac_embs_dim.zarr
export INPUT=data/peaks_indexed.h5ad
export OUTPUT=output/peaks_imputed_ALS.h5ad
bash scripts/run_impute.sh

# 3) （可选）聚类指标；需 obs['cell_type']
# 按 scripts/evaluate_clustering.py 的 CLI 参数调用
```

## 读取结果

```python
import scanpy as sc
ad = sc.read_h5ad("output/peaks_imputed_ALS.h5ad")
X = ad.X            # imputed scores
B = ad.obsm["latent"]
Z = ad.varm["topic"]
A = ad.uns["ALS_A"]
```

## 消融（real / random / shuffled E）

```bash
python chrombert_scai/impute_als_ablation.py \
  --input data/peaks_indexed.h5ad \
  --output output/ablation_shuffled.h5ad \
  --emb_dir $EMB_DIR \
  --emb_dim 768 \
  --emb_mode shuffled \
  --emb_ablation_seed 42 \
  --topics 64 --lambdaA 1 --lambdaB 1 \
  --als_iters 20 --cg_iters 500 --binary --seed 42
```
