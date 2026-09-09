# 模型结构与实现

## 目标

对 peak × cell 目标矩阵 \(Y\)（通常为二值化后的 TF-IDF）：

\[
Y \approx E A B,\quad
Z = E A\ (\text{peak latent}),\quad
B\ (\text{cell latent})
\]

- \(E \in \mathbb{R}^{P\times L}\)：固定 ChromBERT 区域表征（不训练）
- \(A \in \mathbb{R}^{L\times K}\)：可训练
- \(B \in \mathbb{R}^{K\times C}\)：可训练
- 正则：\(\lambda_A\|A\|_F^2 + \lambda_B\|B\|_F^2\)

## \(Y\) 构造（`impute_als.py`）

1. 读入 h5ad → 转置为 peak × cell
2. 可选 `--binary`：`(X > 0)`
3. `sklearn.TfidfTransformer(smooth_idf=False, norm='l2')`，并对 `idf_ -= 1`

## ALS 更新

预计算 \(G_E = E^\top E\)，\(G_Y = E^\top Y\)：

- **A 步**：解 \(G_E A (BB^\top) + \lambda_A A = G_Y B^\top\)（矩阵共轭梯度）
- **B 步**：解 \((A^\top G_E A + \lambda_B I) B = A^\top G_Y\)

默认固定迭代次数，无早停。初始化：\(A\sim\mathcal{N}(0,1)\times 0.02\)，再闭式求 \(B\)。

## 输出注意

- 预测 **无非负约束**（clip 代码默认注释掉），分数可为负
- 这是相对 / TF-IDF-like accessibility，**不是** fragment counts，也未经概率校准
- 下游若 `log1p`，需先处理负值

## 外推未见峰

仅用训练峰拟合得到 `ALS_A` 与 `latent`（\(B^\top\)）后：

```python
imputed_test = (E_test @ A @ B.T).T   # cells × test_peaks
```

严格留出要求：测试峰行不得进入训练时的 \(Y\)、\(E^\top E\)、TF-IDF 与 loss。

## 对照实现

| 文件 | 作用 |
|------|------|
| `impute_als.py` | 生产 ALS |
| `impute_als_ablation.py` | `--emb_mode {real,random,shuffled}` |
| `impute_ridge.py` | 早期 \(Y\approx E W\) ridge（非主路径） |

## 常用超参（基准实验）

`emb_dim=768, topics=128, lambdaA=1, lambdaB=0.1, als_iters=20, cg_iters=500, binary`
