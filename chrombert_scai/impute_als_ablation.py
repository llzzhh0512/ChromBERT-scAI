#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, argparse, math, zarr
import numpy as np
import pandas as pd
import scanpy as sc

def binarize_matrix(X):
    return (X > 0).astype(np.int8)

def tfidf_transform(X_bin):
    from sklearn.feature_extraction.text import TfidfTransformer
    tfidf = TfidfTransformer(smooth_idf=False, norm='l2')
    tfidf = tfidf.fit(X_bin.T)
    tfidf.idf_ -= 1        
    X_tfidf = tfidf.transform(X_bin.T).T
    return X_tfidf.toarray().astype(np.float32)

def precompute_GE(E, blk=20000):
    P, L = E.shape
    GE = np.zeros((L, L), dtype=np.float32)
    for i in range(0, P, blk):
        j = min(P, i + blk)
        GE += E[i:j].T @ E[i:j]
    return GE

def block_ETY(E, Y, blk=20000):
    # GY = E^T Y, 分块计算
    P, L = E.shape
    _, C = Y.shape
    GY = np.zeros((L, C), dtype=np.float32)
    for i in range(0, P, blk):
        j = min(P, i + blk)
        GY += E[i:j].T @ Y[i:j]
    return GY

def mat_cg(Aop, B, X0=None, max_iter=200, tol=1e-5, verbose=False):
    """
    矩阵共轭梯度：解 Aop(X)=B（Frobenius 内积）
    适合 A 步的 Sylvester 形算子：X -> GE @ (X @ (B B^T)) + λ_A X
    """
    X = np.zeros_like(B) if X0 is None else X0.copy()
    R = B - Aop(X)
    P = R.copy()
    r2 = float((R*R).sum())
    if r2 < tol*tol: return X
    for it in range(1, max_iter+1):
        AP = Aop(P)
        denom = float((P*AP).sum()) + 1e-12
        alpha = r2 / denom
        X += alpha * P
        R -= alpha * AP
        r2_new = float((R*R).sum())
        if verbose:
            print(f"  CG it={it:03d}  rel_res={math.sqrt(r2_new)/(math.sqrt(r2)+1e-12):.3e}")
        if r2_new < tol*tol: break
        beta = r2_new / (r2 + 1e-12)
        P = R + beta * P
        r2 = r2_new
    return X

def sample_mae(E, AB, Y, idx_p, idx_c):
  """在固定 peak/cell 子集上计算 sample MAE。"""
  pred = E[idx_p] @ AB[:, idx_c]
  return float(np.abs(pred - Y[idx_p][:, idx_c]).mean())


def find_stable_iter(mae_history, rel_tol=0.01, window=3):
  """当连续 window 轮相对变化 < rel_tol 时，认为趋于稳定。"""
  if len(mae_history) < window + 1:
    return None
  vals = [x["sample_mae"] for x in mae_history]
  for i in range(window, len(vals)):
    recent = vals[i - window:i + 1]
    rel_change = abs(recent[-1] - recent[0]) / (abs(recent[0]) + 1e-12)
    if rel_change < rel_tol:
      return mae_history[i]["iter"]
  return mae_history[-1]["iter"]


def lowrank_impute_no_graph(
    E, Y, n_topics=64, lambda_A=1e-2, lambda_B=1e-2,
    als_iters=10, cg_iters=200, cg_tol=1e-5, seed=42,
    mae_sample_size=1000, log_every=1,
):
    rng = np.random.default_rng(seed)
    P, L = E.shape
    _, C = Y.shape
    K = n_topics
    mae_history = []

    idx_p = rng.integers(0, P, size=mae_sample_size)
    idx_c = rng.integers(0, C, size=mae_sample_size)

    GE = precompute_GE(E)
    GY = block_ETY(E, Y)

    A = (rng.standard_normal((L, K)).astype(np.float32)) * 0.02
    M0 = A.T @ GE @ A + lambda_B * np.eye(K, dtype=np.float32)
    RHS0 = A.T @ GY
    B  = np.linalg.solve(M0 + 1e-6*np.eye(K, dtype=np.float32), RHS0)

    AB = A @ B
    init_mae = sample_mae(E, AB, Y, idx_p, idx_c)
    mae_history.append({"iter": 0, "sample_mae": init_mae})
    print(f"[ALS 0/{als_iters}] sample MAE = {init_mae:.4e}")

    for t in range(1, als_iters+1):
        # --- A step: (GE) A (BB^T) + λ_A A = GY B^T
        BBT = B @ B.T                               # (K,K)
        RHS_A = GY @ B.T                            # (L,K)
        def Aop(X): return GE @ (X @ BBT) + lambda_A * X
        A = mat_cg(Aop, RHS_A, X0=A, max_iter=cg_iters, tol=cg_tol, verbose=False)

        # --- B step: (A^T GE A + λ_B I) B = A^T GY
        M = A.T @ GE @ A + lambda_B * np.eye(K, dtype=np.float32)   # (K,K)
        RHS = A.T @ GY                                              # (K,C)
        B  = np.linalg.solve(M + 1e-6*np.eye(K, dtype=np.float32), RHS)

        if t % log_every == 0 or t == als_iters:
            AB = A @ B
            err = sample_mae(E, AB, Y, idx_p, idx_c)
            mae_history.append({"iter": t, "sample_mae": err})
            print(f"[ALS {t}/{als_iters}] sample MAE = {err:.4e}")

    AB = A @ B                                     # (L,C)
    X_pred = np.zeros((P, C), dtype=np.float32)
    blk = max(20000, min(P, 50000))
    for i in range(0, P, blk):
        j = min(P, i + blk)
        X_pred[i:j] = E[i:j] @ AB
    return X_pred, A, B, mae_history
def prepare_embedding(E, mode, seed=42):
    """mode: real | random | shuffled"""
    rng = np.random.default_rng(seed)
    P, L = E.shape

    if mode == "real":
        return E.copy()

    if mode == "random":
        # 推荐：匹配真实 E 的 per-dimension 均值/方差，避免尺度差异干扰
        mu = E.mean(axis=0)
        std = E.std(axis=0) + 1e-8
        E_rand = rng.standard_normal((P, L)).astype(np.float32)
        return E_rand * std + mu

    if mode == "shuffled":
        perm = rng.permutation(P)
        return E[perm].copy()

    raise ValueError(f"Unknown emb_mode: {mode}")
# ---------- CLI ----------
def main():
    parser = argparse.ArgumentParser(description="Low-rank imputation without graph regularization")
    parser.add_argument('--input', required=True, help='.h5ad 或 .txt（行=peaks/coords，列=cells）')
    parser.add_argument('--output', required=True, help='输出 .h5ad')
    parser.add_argument('--emb_dir', required=True, help='zarr 根目录（含 full 或 PCA/<dim>）')
    parser.add_argument('--emb_dim', type=int, default=768, help='embedding 维度（768 或 PCA 维）')
    parser.add_argument('--topics', type=int, default=30)
    parser.add_argument('--lambdaA', type=float, default=1e-2)
    parser.add_argument('--lambdaB', type=float, default=1e-2)
    parser.add_argument('--als_iters', type=int, default=10)
    parser.add_argument('--cg_iters', type=int, default=200)
    parser.add_argument('--cg_tol',   type=float, default=1e-5)
    parser.add_argument('--binary', action='store_true', help='把输入二值化后再做 TF-IDF（默认用原始0/计数→TF-IDF）')
    parser.add_argument('--mae_log', type=str, default=None, help='sample MAE 曲线输出 csv；默认与 output 同目录')
    parser.add_argument('--mae_sample_size', type=int, default=1000, help='每轮 sample MAE 抽样数')
    parser.add_argument('--log_every', type=int, default=1, help='每隔几轮 ALS 记录一次 MAE')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--emb_mode', choices=['real', 'random', 'shuffled'], default='real')
    parser.add_argument('--emb_ablation_seed', type=int, default=42,
                    help='random/shuffled 消融用的随机种子')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    if args.input.endswith('.h5ad'):
        ad = sc.read_h5ad(args.input)
        X = ad.X.toarray().T        # (P,C)
        emb_idx = ad.varm['indices']
   
   
    if args.emb_dim == 768:
        E_src = zarr.open_array(f"{args.emb_dir}/full", mode="r")
    else:
        E_src = zarr.open_array(f"{args.emb_dir}/PCA/{args.emb_dim}", mode="r")
    E = np.asarray(E_src[emb_idx], dtype=np.float32)
    E = prepare_embedding(E, mode=args.emb_mode, seed=args.emb_ablation_seed)
    print(f"[ablation] emb_mode={args.emb_mode}, seed={args.emb_ablation_seed}, E shape={E.shape}")
    X_bin = binarize_matrix(X) if args.binary else X
    Y = tfidf_transform(X_bin)                        # (P,C) float32

    print("Running low-rank ALS ...")
    X_pred, A, B, mae_history = lowrank_impute_no_graph(
        E=E, Y=Y, n_topics=args.topics,
        lambda_A=args.lambdaA, lambda_B=args.lambdaB,
        als_iters=args.als_iters, cg_iters=args.cg_iters, cg_tol=args.cg_tol,
        seed=args.seed, mae_sample_size=args.mae_sample_size, log_every=args.log_every,
    )

    dataset_name = os.path.splitext(os.path.basename(args.input))[0]
    mae_log = args.mae_log or os.path.join(
        os.path.dirname(args.output) or ".", "als_mae_history.csv"
    )
    df_mae = pd.DataFrame(mae_history)
    df_mae["dataset"] = dataset_name
    stable_iter = find_stable_iter(mae_history)
    df_mae["stable_iter"] = stable_iter
    df_mae.to_csv(mae_log, index=False)

    init_mae = mae_history[0]["sample_mae"]
    final_mae = mae_history[-1]["sample_mae"]
    print(f"[summary] dataset={dataset_name}")
    print(f"[summary] sample MAE: {init_mae:.4e} -> {final_mae:.4e}")
    if stable_iter is not None:
        print(f"[summary] stable around iter {stable_iter}")
    print(f"[summary] mae log saved to {mae_log}")

    # 写出
    out = ad.copy()
    out.X = X_pred.T.astype('float16')
    out.varm['topic'] = E@A
    out.uns['ALS_A'] = A
    out.obsm['latent'] = B.T
    # mae_history 已写入 csv；list[dict] 无法直接存入 uns
    out.write_h5ad(args.output)

if __name__ == "__main__":
    main()
