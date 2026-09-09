# 本目录与完整工程的关系

| 完整工程路径 | 本发布包对应 |
|--------------|--------------|
| `1.source_from_zyx/2.benchmark/ChromBERT_topic/source/impute_ALS.py` | `chrombert_scai/impute_als.py` |
| `2.plot/1.model_test/source/impute_ALS.py` | `chrombert_scai/impute_als_ablation.py` |
| `1.source_from_zyx/2.benchmark/ChromBERT/source/impute.py` | `chrombert_scai/impute_ridge.py` |
| `1.source_from_zyx/3.utils/1.cluster_metrics/1.cluster_metric.py` | `scripts/evaluate_clustering.py` |
| `1.source_from_zyx/1.datasets/` | 仅文档：`docs/datasets.md` |
| `chrombert_source/data/config/hg38_6k_1kb_region.bed` | `annotations/hg38_6k_1kb_region.bed` |
| `chrombert_source/mouse_data/config/mm10_5k_1kb_region.bed` | `annotations/mm10_5k_1kb_region.bed` |
| `paper_audit_20260908_server/` | 不纳入 GitHub 包（审计笔记） |

上传建议：将 `ChromBERT-scAI_github/` 作为独立 git 仓库根目录，不要把整个 `12.atac_impute` 推上去。
