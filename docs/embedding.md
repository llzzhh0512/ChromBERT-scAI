# ChromBERT 区域 embedding（\(E\)）

## 本仓库提供什么

仅提供 zarr **元数据**（`docs/embedding_zarr_meta/`），用于说明数组形状，**不含** embedding 数值。

| 数组 | 形状（hg38） | dtype |
|------|--------------|-------|
| `full` | 2,137,894 × 768 | float16 |
| `PCA/32` 等 | 2,137,894 × {32,50,100} | float64 |

小鼠：`atac_embs_dim_mm10.zarr`，`full` 约为 1,531,858 × 768。

## 使用方式

`impute_als.py`：

```text
emb_dim == 768 → {emb_dir}/full
else           → {emb_dir}/PCA/{emb_dim}
E = E_src[adata.varm['indices']]
```

生产实验默认 **`--emb_dim 768`（full，不做 PCA）**。

## 获取方式

embedding 由 ChromBERT 预训练管线生成，工程内常引用路径形如：

- hg38: `atac_embs_dim.zarr`
- mm10: `atac_embs_dim_mm10.zarr`

发布前请补充：**下载链接 / DOI / 生成脚本 / license**。  
在说明中应写清预训练数据类型（调控因子 ChIP-seq 共关联等），避免未经核验写成 “sequence-only LM”。

## Peak 与 bin 对齐

用 `scripts/prepare_peak_indices.py` 把 peak 中心映射到 `bin_index`，写入 `adata.varm['indices']`。

标准 bin BED **已随本仓库提供**：

| Genome | Path | #bins |
|--------|------|------:|
| hg38 | `annotations/hg38_6k_1kb_region.bed` | 2,137,894 |
| mm10 | `annotations/mm10_5k_1kb_region.bed` | 1,530,871 |

格式：`chrom  start  end  bin_index`（无表头）。hg38 行数与 zarr `full` 第一维一致。

原始拷贝来源见 `annotations/README.md`。
