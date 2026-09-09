import os
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import Counter, defaultdict
from sklearn.metrics import f1_score, precision_score, recall_score
from sklearn.metrics.cluster import adjusted_rand_score, normalized_mutual_info_score, adjusted_mutual_info_score, homogeneity_score, silhouette_score

import warnings
warnings.filterwarnings("ignore")

def binarize_matrix(X):
    return (X > 0).astype(np.int8)

def tfidf_transform(X_bin):
    from sklearn.feature_extraction.text import TfidfTransformer
    tfidf = TfidfTransformer(smooth_idf=False, norm='l2')
    tfidf = tfidf.fit(X_bin.T)
    tfidf.idf_ -= 1        
    X_tfidf = tfidf.transform(X_bin.T).T
    return X_tfidf.toarray().astype(np.float32)


def map_clusters_to_labels_by_majority(labels_true, labels_pred):
    cluster_to_true = defaultdict(list)
    for true_label, pred_label in zip(labels_true, labels_pred):
        cluster_to_true[pred_label].append(true_label)
    cluster_to_label = {cluster: Counter(true_list).most_common(1)[0][0]
                       for cluster, true_list in cluster_to_true.items()}
    aligned_pred = [cluster_to_label[cl] for cl in labels_pred]
    return aligned_pred, cluster_to_label

def get_metric(adata, odir, oname, n_neighbors=10, resolution=1.5, celltypes=None, use_latent=False):
    if not use_latent:
        sc.pp.scale(adata)
        adata.X = np.nan_to_num(adata.X)
        sc.tl.pca(adata, n_comps=100, svd_solver='arpack')
        sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=n_neighbors)
        # sc.tl.tsne(adata, random_state=0, use_rep='X_pca')
    else:
        assert "latent" in adata.obsm, "latent not found in adata.obsm"
        sc.pp.neighbors(adata, use_rep="latent", n_neighbors=n_neighbors, random_state=0)
        # sc.tl.tsne(adata, random_state=0, use_rep="latent")

    sc.tl.umap(adata, random_state=0)
    sc.tl.leiden(adata, resolution=resolution)

    true_labels = adata.obs['cell_type'].to_numpy()
    cluster_labels = adata.obs['leiden'].astype(int).to_numpy()
    aligned_labels, _ = map_clusters_to_labels_by_majority(true_labels, cluster_labels)
    adata.obs['leiden_aligned'] = aligned_labels

    with PdfPages(os.path.join(odir, f"{oname}.pdf")) as pdf:
        for method in ["leiden", "leiden_aligned"]:
            for func in [sc.pl.umap]:
                fig = func(adata, color=[method, "cell_type"], show=False, return_fig=True)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
    
    ari = adjusted_rand_score(true_labels, cluster_labels)
    nmi = normalized_mutual_info_score(true_labels, cluster_labels)
    ami = adjusted_mutual_info_score(true_labels, cluster_labels)
    homo = homogeneity_score(true_labels, cluster_labels)
    silh = silhouette_score(adata.obsm['X_umap'], cluster_labels)
    f1 = f1_score(true_labels, aligned_labels, average='micro')
    
    list_metrics = []
    if celltypes is not None:
        true_labels = pd.Series(true_labels)
        aligned_labels = pd.Series(aligned_labels)
        for celltype in celltypes:
            inner_true_labels = (true_labels == celltype).astype(int)
            inner_aligned_labels = (aligned_labels == celltype).astype(int)
            celltype_f1 = f1_score(inner_true_labels, inner_aligned_labels, zero_division=0)
            precision = precision_score(inner_true_labels, inner_aligned_labels, zero_division=0)
            recall = recall_score(inner_true_labels, inner_aligned_labels, zero_division=0)
            list_metrics.append({
                "ARI": ari, "NMI": nmi, "total_F1": f1, "celltype": celltype,
                "F1": celltype_f1, "Precision": precision, "Recall": recall
            })
    else:
        list_metrics.append({
                "ARI": ari, "NMI": nmi, "AMI": ami, "Homogeneity": homo, "Silhouette": silh, "total_F1": f1
            })
    return list_metrics

def main():
    parser = argparse.ArgumentParser(description="Compute clustering metrics for imputed ATAC data with optional latent representations.")
    parser.add_argument('--data', type=str, required=True, help='File prefix (dataset name)')
    parser.add_argument('--odir', type=str, required=True, help='Output directory for results')
    parser.add_argument('--rare_ratio', type=float, default=0.03, help='Threshold ratio for rare cell types')
    parser.add_argument('--count_cell', action="store_true", default=False, help='Count cell type')
    parser.add_argument('--methods', type=str, nargs='+', default=["ChromBERT", "ChromBERT_latent", "scOpen", "scOpen_latent", "scAGDE", "scAGDE_latent", "SCALE", "SCALE_latent", "cisTopic", "cisTopic_latent"], help='List of methods to evaluate')
    args = parser.parse_args()
    
    path_scopen = "/mnt/Storage2/home/lizhanhao/project/12.atac_impute/1.source_from_zyx/2.benchmark/scopen/output"
    path_chrombert = "/mnt/Storage2/home/lizhanhao/project/12.atac_impute/1.source_from_zyx/2.benchmark/ChromBERT/output"
    path_chrombert_topic = "/mnt/Storage/home/zhangyuxuan/projects/scATACseqImpute/2.benchmark/ChromBERT_topic/output"
    path_scagde = "/mnt/Storage2/home/lizhanhao/project/12.atac_impute/1.source_from_zyx/2.benchmark/scagde/output"
    path_scale = "/mnt/Storage2/home/lizhanhao/project/12.atac_impute/1.source_from_zyx/2.benchmark/scale/output"
    path_cistopic = "/mnt/Storage/home/zhangyuxuan/projects/scATACseqImpute/2.benchmark/cisTopic/output"
    path_baseline = "/mnt/Storage2/home/lizhanhao/project/12.atac_impute/1.source_from_zyx/1.datasets/output"
    file_baseline = os.path.join(path_baseline, f"{args.data}")
    file_scopen = os.path.join(path_scopen, args.data, "imputed.h5ad")
    file_chrombert = os.path.join(path_chrombert, args.data, "imputed.h5ad")
    file_scagde = os.path.join(path_scagde, args.data, "imputed.h5ad")
    file_scale = os.path.join(path_scale, args.data, "adata.h5ad")
    file_cistopic = os.path.join(path_cistopic, f"{args.data}",  "imputed.h5ad")
    file_chrombert_topic = os.path.join(path_chrombert_topic, args.data, "imputed_ALS.h5ad")
    #file_chrombert_topic = os.path.join(path_chrombert_topic, args.data, "imputed.h5ad")
    os.makedirs(args.odir, exist_ok=True)
    # adata = sc.read_h5ad(file_chrombert)
    celltypes = None
    if args.count_cell:
        try:
            total_counts = adata.obs['cell_type'].value_counts().sum()
            df = pd.DataFrame(adata.obs['cell_type'].value_counts()).reset_index()
            df.columns = ['cell_type', 'count']
            df['ratio'] = df['count'] / total_counts
            celltypes = df.query("ratio <= 0.03")['cell_type'].tolist()
            dict_ct_to_ratio = dict(zip(df['cell_type'], df['ratio']))
        except:
            assert False, "cell_type column not found in adata"

    list_metrics = []
    for n_neighbors in [15, 10, 5]:
        for resolution in [0.5, 1, 1.5]:
            for method in args.methods:
                oname = f"{args.data.split('/')[-1]}_{method}_n{n_neighbors}_r{resolution}"
                print(oname)
                if method == "ChromBERT_latent":
                    adata = sc.read_h5ad(file_chrombert)
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes, use_latent=True)
                elif method == "ChromBERT_topic_latent":
                    adata = sc.read_h5ad(file_chrombert_topic)
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes, use_latent=True)
                if method == "cisTopic_latent":
                    adata = sc.read_h5ad(file_cistopic)
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes, use_latent=True)
                elif method == "scAGDE_latent":
                    adata = sc.read_h5ad(file_scagde)
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes, use_latent=True)
                elif method == "SCALE_latent":
                    adata = sc.read_h5ad(file_scale)
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes, use_latent=True)
                elif method == "scOpen_latent":
                    # adata = sc.read_h5ad(file_scopen)
                    df_barcodes = pd.read_csv(os.path.join(path_scopen, args.data, "imputed_barcodes.txt"), sep='\t', index_col=0)
                    adata.obsm['latent'] = df_barcodes.values.T
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes, use_latent=True)
                elif method == "scAGDE":
                    adata = sc.read_h5ad(file_scagde)
                    imputed = adata.obsm["impute"]
                    adata = anndata.AnnData(imputed,obs=adata.obs,var=adata[:,adata.var["is_selected"] == 1].var)
                    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
                    sc.pp.log1p(adata)
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes)
                elif method == "SCALE":
                    adata = sc.read_h5ad(file_scale)
                    imputed = adata.obsm["impute"]
                    adata = anndata.AnnData(imputed,obs=adata.obs,var=adata.var)
                    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
                    sc.pp.log1p(adata)
                    # sc.pp.highly_variable_genes(adata, min_mean=0.05, max_mean=1.5, min_disp=.5)
                    # adata = adata[:, adata.var['highly_variable']].copy()
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes)
                elif method == "ChromBERT":
                    adata = sc.read_h5ad(file_chrombert)
                    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
                    sc.pp.log1p(adata)
                    # sc.pp.highly_variable_genes(adata, min_mean=0.05, max_mean=1.5, min_disp=.5)
                    # adata = adata[:, adata.var['highly_variable']].copy()
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes)
                elif method == "ChromBERT_topic":
                    adata = sc.read_h5ad(file_chrombert_topic)
                    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
                    sc.pp.log1p(adata)
                    # sc.pp.highly_variable_genes(adata, min_mean=0.05, max_mean=1.5, min_disp=.5)
                    # adata = adata[:, adata.var['highly_variable']].copy()
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes)
                elif method == "cisTopic":
                    adata = sc.read_h5ad(file_cistopic)
                    adata.X = adata.X * 1e6
                    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
                    sc.pp.log1p(adata)
                    # sc.pp.highly_variable_genes(adata, min_mean=0.05, max_mean=1.5, min_disp=.5)
                    # adata = adata[:, adata.var['highly_variable']].copy()
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes)
                elif method == "scOpen":
                    adata = sc.read_h5ad(file_scopen)
                    adata.X = adata.X.astype("float32")
                    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
                    sc.pp.log1p(adata)
                    # sc.pp.highly_variable_genes(adata, min_mean=0.05, max_mean=1.5, min_disp=.5)
                    # adata = adata[:, adata.var['highly_variable']].copy()
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes)
                elif method == "baseline":
                    adata = sc.read_h5ad(file_baseline)
                    adata.X = adata.X.astype("float32")
                    adata.X = binarize_matrix(adata.X)
                    adata.X = tfidf_transform(adata.X)
                    inner_metric = get_metric(adata, odir=args.odir, oname=oname, resolution=resolution, n_neighbors=n_neighbors, celltypes=celltypes)
                for metric in inner_metric:
                    metric.update({'method': method, "resolution": resolution, "n_neighbors": n_neighbors})
                print(inner_metric)
                list_metrics.extend(inner_metric)
    df_metrics = pd.DataFrame(list_metrics)
    if args.count_cell:
        df_metrics['ratio'] = df_metrics['celltype'].map(dict_ct_to_ratio)
    metrics_file = os.path.join(args.odir, f"metrics_rare_cells.csv")
    if not os.path.exists(metrics_file):
        df_metrics.to_csv(metrics_file, index=False)
    else:
        df_metrics.to_csv(metrics_file, mode='a', header=False, index=False)

if __name__ == "__main__":
    main() 