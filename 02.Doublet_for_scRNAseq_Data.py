######## in python
import scanpy as sc
adata = sc.read.h5ad(file = '/home/data/ssy43/PCa/all.cells.h5ad')
sc.external.pp.scrublet(adata, expected_doublet_rate = 0.05,threshold=0.25,batch_key = "sampleID")
data_dir="/home/data/ssy43/PCa/Outdata/"
adata.obs[["doublet_score","predicted_doublet", "barcode"]].to_csv('{}adata_scrublet.csv'.format(data_dir))

######## in R
# batch.data = read.csv("./Outdata/adata_scrublet.csv",header=T)
# table(batch.data$predicted_doublet)
# meta.data = merged_seurat@meta.data
# meta.data = left_join(meta.data, batch.data)
# row.names(meta.data) = row.names(merged_seurat@meta.data)
# merged_seurat@meta.data = meta.data
# table(merged_seurat$predicted_doublet)
# seurat.filter <- subset(merged_seurat, subset = predicted_doublet =="False")
