#!/usr/bin/env python3

import os
import time
import argparse
import zarr
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from sklearn.decomposition import NMF
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import TruncatedSVD
import scanpy as sc
from sklearn.linear_model import Ridge
import h5py
import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
# from sklearn.linear_model import RidgeCV


def binarize_matrix(X):
    return (X > 0).astype(int)

def tfidf_transform(X_bin):
    model = TfidfTransformer(smooth_idf=False, norm='l2')
    model = model.fit(X_bin.T)
    model.idf_ -= 1
    X_tfidf = model.transform(X_bin.T).T
    return X_tfidf #.toarray()

def ridge_reconstruct_batch(E, X_tfidf, alpha=1.0, batch_size=None):
    n_targets = X_tfidf.shape[1]
    if batch_size is None:
        batch_size = n_targets
    W_list = []
    for i in range(0, n_targets, batch_size):
        start_time = time.time()
        j = min(i + batch_size, n_targets)  # 防止越界
        y_batch = X_tfidf[:, i:j]
        if not isinstance(y_batch, np.ndarray):
            y_batch = y_batch.toarray().astype(np.float16)

        model = Ridge(alpha=alpha, fit_intercept=False)
        model.fit(E, y_batch)
        # print("Best alpha:", model.alpha_)
        W_list.append(model.coef_.T)  # shape: (emb_dim, batch)
        end_time = time.time()
        print(f"Batch {i} reconstruction time: {end_time - start_time} seconds")

    W = np.hstack(W_list)
    print("reconstructed W shape: ", W.shape)            # shape: (emb_dim, n_targets)
    X_pred = E @ W                   # shape: (n_samples, n_targets)
    return W, X_pred


def main():
    parser = argparse.ArgumentParser(description='Impute ATAC-seq data using ChromBERT method')
    parser.add_argument('--input', required=True, help='Input file')
    parser.add_argument('--output', required=True, help='Output path for imputed data')
    parser.add_argument('--coord_index', required=False, default=None, help='Coordinate to index mapping file')
    parser.add_argument('--emb_dir', required=True, help='Directory containing embeddings')
    parser.add_argument('--emb_dim', type=int, default=50, help='Embedding dimension')
    parser.add_argument('--batch_size', type=int, default=None, help='Batch size for ridge reconstruction')
    
    args = parser.parse_args()

   
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    # assert os.path.exists(args.coord_index), f"Coordinate to index mapping file {args.coord_index} does not exist"
    assert os.path.exists(args.input), f"Input file {args.input} does not exist"
    assert os.path.exists(args.emb_dir), f"Embedding directory {args.emb_dir} does not exist"

    # Load data
    if args.coord_index is not None:
        df_coord_to_idx = pd.read_csv(args.coord_index)
        dict_coord_to_idx = dict(zip(df_coord_to_idx['coord'], df_coord_to_idx['bin_index']))
        if args.input.endswith('.txt'):
            df_raw = pd.read_csv(args.input, sep='\t', index_col=0).reset_index(names=['coord'])
            df_raw['bin_index'] = df_raw['coord'].map(dict_coord_to_idx)
            df_raw_valid = df_raw[~df_raw['bin_index'].isna()]
            X = df_raw_valid.drop(columns=['coord', 'bin_index']).values
        elif args.input.endswith('.h5ad'):
            data = sc.read_h5ad(args.input)
            df_raw = pd.DataFrame(data.var.index.tolist(), columns=['coord'])
            df_raw['bin_index'] = df_raw['coord'].map(dict_coord_to_idx)
            df_raw_valid = df_raw[~df_raw['bin_index'].isna()]
            try:
                X = data.X[:, df_raw_valid.index].toarray().T
            except:
                X = data.X[:, df_raw_valid.index].T
            emb_indices = df_raw_valid['bin_index'].values.astype(int)
    else:
        assert args.input.endswith('.h5ad'), f"Input file {args.input} is not a valid file"
        data = sc.read_h5ad(args.input)
        X = data.X.T
        emb_indices = data.varm['indices']
    
    # Load embeddings
    if args.emb_dim == 768:
        emb_source = zarr.open_array(f"{args.emb_dir}/full", mode="r")
    else:
        emb_source = zarr.open_array(f"{args.emb_dir}/PCA/{str(args.emb_dim)}", mode="r")


    print("Running Ridge...")
    
    X_binarized = binarize_matrix(X)
    X_tfidf = tfidf_transform(X_binarized)
    E = emb_source[emb_indices]
    W, X_chrombert = ridge_reconstruct_batch(E, X_tfidf, alpha=1.0, batch_size=args.batch_size)
    X_chrombert = np.where(X_chrombert < 0, 0, X_chrombert)
    print("writing to h5ad...")
    data.X = X_chrombert.T
    data.obsm['latent'] = W.T
    data.write_h5ad(args.output)

if __name__ == "__main__":
    main() 