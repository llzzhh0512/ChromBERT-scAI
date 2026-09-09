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

def lowrank_impute_no_graph(
    E, Y, n_topics=64, lambda_A=1e-2, lambda_B=1e-2,
    als_iters=10, cg_iters=200, cg_tol=1e-5, seed=42
):
    rng = np.random.default_rng(seed)
    P, L = E.shape
    _, C = Y.shape
    K = n_topics

    GE = precompute_GE(E)
    GY = block_ETY(E, Y)

    A = (rng.standard_normal((L, K)).astype(np.float32)) * 0.02
    M0 = A.T @ GE @ A + lambda_B * np.eye(K, dtype=np.float32)
    RHS0 = A.T @ GY
    B  = np.linalg.solve(M0 + 1e-6*np.eye(K, dtype=np.float32), RHS0)

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

        if t % 2 == 0:
            AB = A @ B
            idx_p = rng.integers(0, P, size=1000)
            idx_c = rng.integers(0, C, size=1000)
            err = np.abs((E[idx_p] @ AB[:, idx_c]) - Y[idx_p][:, idx_c]).mean()
            print(f"[ALS {t}/{als_iters}] sample MAE = {err:.4e}")

    AB = A @ B                                     # (L,C)
    X_pred = np.zeros((P, C), dtype=np.float32)
    blk = max(20000, min(P, 50000))
    for i in range(0, P, blk):
        j = min(P, i + blk)
        X_pred[i:j] = E[i:j] @ AB
    return X_pred, A, B

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
    parser.add_argument('--seed', type=int, default=42, help='ALS initialization seed')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    if args.input.endswith('.h5ad'):
        ad = sc.read_h5ad(args.input)
        X = ad.X.toarray().T        # (P,C)
        emb_idx = ad.varm['indices']
   
   
    if args.emb_dim == 768:
        E_src = zarr.open_array(f"{args.emb_dir}/full", mode="r")
    else:
        E_src = zarr.open_array(f"{args.emb_dir}/PCA/{args.emb_dim}", mode="r")
    E = E_src[emb_idx]                               # (P,L)
    E = E.astype(np.float32, copy=False)

    X_bin = binarize_matrix(X) if args.binary else X
    Y = tfidf_transform(X_bin)                        # (P,C) float32

    print("Running low-rank ALS ...")
    X_pred, A, B = lowrank_impute_no_graph(
        E=E, Y=Y, n_topics=args.topics,
        lambda_A=args.lambdaA, lambda_B=args.lambdaB,
        als_iters=args.als_iters, cg_iters=args.cg_iters, cg_tol=args.cg_tol,
        seed=args.seed
    )

    # X_pred = np.where(X_pred < 0, 0, X_pred)

    # 写出
    out = ad.copy()
    out.X = X_pred.T.astype('float16')
    out.varm['topic'] = E@A
    out.uns['ALS_A'] = A
    out.obsm['latent'] = B.T
    out.write_h5ad(args.output)

if __name__ == "__main__":
    main()
