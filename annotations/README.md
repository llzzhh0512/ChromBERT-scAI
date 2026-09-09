# ChromBERT genomic bin BED files (tracked in this repo)

| Genome | File in repo | #bins | Format |
|--------|--------------|------:|--------|
| hg38 | `annotations/hg38_6k_1kb_region.bed` | 2,137,894 | `chrom  start  end  bin_index` |
| mm10 | `annotations/mm10_5k_1kb_region.bed` | 1,530,871 | same |

Copied from ChromBERT source configs:

- `/mnt/Storage2/home/lizhanhao/chrombert_source/data/config/hg38_6k_1kb_region.bed`
- `/mnt/Storage2/home/lizhanhao/chrombert_source/mouse_data/config/mm10_5k_1kb_region.bed`

## Usage

```bash
# human
python scripts/prepare_peak_indices.py \
  --input data/peaks.h5ad \
  --bins annotations/hg38_6k_1kb_region.bed \
  --output data/peaks_indexed.h5ad \
  --drop-unmapped

# mouse
python scripts/prepare_peak_indices.py \
  --input data/peaks_mm10.h5ad \
  --bins annotations/mm10_5k_1kb_region.bed \
  --output data/peaks_mm10_indexed.h5ad \
  --drop-unmapped
```
