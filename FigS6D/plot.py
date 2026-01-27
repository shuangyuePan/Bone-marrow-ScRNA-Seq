import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import anndata as ad
import warnings

# Suppress warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# Input directory configuration (using relative path to FigS6D)
data_dir = "../data" 
export_dir = os.path.join(data_dir, "export_scvelo")

# Load data (Assuming all files exist)
X = sc.read_mtx(os.path.join(export_dir, "counts.mtx")).X
adata = ad.AnnData(X=X.transpose().tocsr())

cell_meta = pd.read_csv(os.path.join(export_dir, "metadata.csv"))
adata.obs = cell_meta
adata.obs.index = adata.obs["barcode"]

with open(os.path.join(export_dir, "gene_names.csv"), "r") as f:
    gene_names = f.read().splitlines()
adata.var.index = gene_names

pca = pd.read_csv(os.path.join(export_dir, "pca.csv"))
pca.index = adata.obs.index
adata.obsm["X_pca"] = pca.to_numpy()
adata.obsm["X_umap"] = np.vstack((adata.obs["UMAP_1"].to_numpy(), adata.obs["UMAP_2"].to_numpy())).T
# adata.obsm["X_tsne"] = np.vstack((adata.obs["TSNE_1"].to_numpy(), adata.obs["TSNE_2"].to_numpy())).T

# Load loom files (Assuming loom files exist)
loom_files = [f for f in os.listdir(data_dir) if f.endswith('.loom')]
ldata_list = []
for lf in loom_files:
    path = os.path.join(data_dir, lf)
    ldata = sc.read(path, cache=True)
    cellids = ldata.obs.index.astype(str).tolist()
    samples = [cid.split(":")[0].split("_")[0].lower() for cid in cellids]
    cb = [cid.split(":")[1].rstrip("x") for cid in cellids]
    barcodes = [samples[i] + "_" + cb[i] for i in range(len(cb))]
    ldata.obs.index = barcodes
    ldata.var_names_make_unique()
    ldata_list.append(ldata)

# Merge data (Assuming ldata_list is not empty)
ldata = ldata_list[0]
for i in range(1, len(ldata_list)):
    ldata = ldata.concatenate([ldata_list[i]])

adata = scv.utils.merge(adata, ldata)

# Preprocessing and Velocity calculation
scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

# Plotting Fig.S6D: Velocity Stream on UMAP
output_file = "FigS6D_velocity_stream_umap.pdf"
scv.pl.velocity_embedding_stream(
    adata, 
    basis='umap', 
    color='celltype', 
    size=20, 
    alpha=0.9, 
    dpi=300, 
    save=output_file, 
    legend_loc='on data',
    title=''
)
print(f"Fig.S6D generated: {output_file}")
