########## Run
import scanpy as sc
# import bbknn
import numpy as np
# import scvi
# from rich import print
# from scib_metrics.benchmark import Benchmarker
# from scvi.model.utils import mde
# from scvi_colab import install
from matplotlib import pyplot as plt
import pandas as pd
# !pip install scanorama==1.7.3
import scanorama
# !pip install bbknn
import bbknn
import scvi
from rich import print
from scib_metrics.benchmark import Benchmarker
from scvi.model.utils import mde

## set
sc.settings.verbosity = 3 
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

input_data="/home/data/ssy43/PCa/Step1.seurat.data_for_python.h5ad"
plot_dir='/home/data/ssy43/PCa/Python.vis/'
data_dir='/home/data/ssy43/PCa/'
sc.settings.figdir = plot_dir
batch_ID = "Data.sets"
celltype = "celltype"

### 1.input data
adata = sc.read_h5ad(input_data)
adata
adata.raw = adata  # keep full dimension safe
adata.X = adata.X.tocsr()
adata.layers["counts"] = adata.X.copy()

### 2.Run methods
#### 2.1 Scanorama
# List of adata per batch
batch_cats = adata.obs[batch_ID].cat.categories.tolist()
adata_list = [adata[adata.obs[batch_ID] == b].copy() for b in batch_cats]
# adatas=[adata_list[0],adata_list[1]]
scanorama.integrate_scanpy(adata_list)

adata.obsm["X_scanorama"] = np.zeros((adata.shape[0], adata_list[0].obsm["X_scanorama"].shape[1]))
for i, b in enumerate(batch_cats):
    adata.obsm["X_scanorama"][adata.obs[batch_ID] == b] = adata_list[i].obsm["X_scanorama"]
adata
sc.pp.neighbors(adata, use_rep="X_scanorama")
sc.tl.umap(adata)
adata.obsm["Scanorama"] = adata.obsm["X_umap"]

#### 2.2 BBKNN
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,batch_key=batch_ID)
#sc.pl.highly_variable_genes(adata)
adata.var[adata.var['highly_variable']]
sc.tl.pca(adata, svd_solver='arpack')
# perform bbknn
sc.external.pp.bbknn(adata, batch_key= batch_ID)
# re-compute UMAP
sc.tl.umap(adata)
adata.obsm["BBKNN"] = adata.obsm["X_umap"]

#### 2.3 scVI
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_ID)
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
#Now we train scVI. This should take a couple of minutes on a Colab session
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata)
adata.obsm["scVI"] = mde(adata.obsm["X_scVI"])

#### 2.4 scANVI
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="seurat_annotations",
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=20, n_samples_per_label=100)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
adata.obsm["scANVI"] = mde(adata.obsm["X_scANVI"])

#### 2.5 Combat
# create a new object with lognormalized counts
adata_combat = sc.AnnData(X=adata.raw.X, var=adata.raw.var, obs = adata.obs)
# first store the raw data 
adata_combat.raw = adata_combat
# run combat
sc.pp.combat(adata_combat, key=batch_ID)
sc.pp.highly_variable_genes(adata_combat)
high_var_genes = set(adata.var_names[adata.var['highly_variable']])
adata_combat.var['highly_variable'] = [gene in high_var_genes for gene in adata_combat.var_names]
# Then we run the regular steps of dimensionality reduction on the combat corrected data
sc.pp.pca(adata_combat, n_comps=30, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata_combat)
sc.tl.umap(adata_combat)
adata.obsm["Combat"] = adata_combat.obsm["X_umap"]

### 3. scib
bm = Benchmarker(
    adata,
    batch_key="Data.sets",
    label_key="celltype",
    embedding_obsm_keys=["Unintegrated",'Harmony', 'LIGER', 'CCA', 'RPCA', 'fastMNN', "Connos", "Scanorama", "BBKNN", "scVI", "scANVI", "Combat"],
    n_jobs=-1,
)
bm.benchmark()
bm.plot_results_table(min_max_scale=False)

df = bm.get_results(min_max_scale=False)
print(df)
