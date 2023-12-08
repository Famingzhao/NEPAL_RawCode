rm(list=ls())
library(Seurat)
# remotes::install_github("satijalab/seurat-data",force = TRUE)
# library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(readr)
library(future)
library(cowplot)
library(scater)
library(data.table)
library(stringr)
library(ggpubr)
library(tidyverse)
library(ggsci)
library(SingleR)
# change the current plan to access parallelization
plan("multisession", workers = 4)
plan()

### 1.input
st.data <- read_rds("./Outdata/Step1.GSM7211257_NEPC_ARPC.rds")
SpatialFeaturePlot(st.data, features = "nCount_Spatial") + 
  theme(legend.position = "right")

### 2.SCT
st.data <- SCTransform(st.data, assay = "Spatial", verbose = FALSE)
st.data <- st.data %>%  RunPCA( assay = "SCT", verbose = FALSE) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  RunUMAP(reduction = "pca", dims = 1:30)
DimPlot(st.data)
SpatialFeaturePlot(st.data, features = "nCount_Spatial") + 
  theme(legend.position = "right")

for (res in c(0.05, 0.2, 0.3, 0.5, 0.8, 1)) {
  # res=0.01
  print(res)
  st.data <- FindClusters(st.data, resolution = res, algorithm = 1, verbose = FALSE)
}

p1.cluster_umap <- plot_grid(ncol = 3, 
                             DimPlot(st.data, reduction = "umap", group.by = "SCT_snn_res.0.05", label = T) & NoAxes(), 
                             DimPlot(st.data, reduction = "umap", group.by = "SCT_snn_res.0.2", label = T)& NoAxes(),
                             DimPlot(st.data, reduction = "umap", group.by = "SCT_snn_res.0.3", label = T)& NoAxes(),
                             DimPlot(st.data, reduction = "umap", group.by = "SCT_snn_res.0.4", label = T) & NoAxes(), 
                             DimPlot(st.data, reduction = "umap", group.by = "SCT_snn_res.0.5", label = T) & NoAxes(),
                             DimPlot(st.data, reduction = "umap", group.by = "SCT_snn_res.0.8", label = T) & NoAxes()
)
p1.cluster_umap
ggsave(plot = p1.cluster_umap, width = 14 ,height = 8,
       filename="./Outplot/Step2.SCT_snn_clusting_PCa.pdf")

st.data$Niche = paste0("Niche_",st.data$SCT_snn_res.0.3)
st.data$Niche = factor(st.data$Niche,levels = levels(as.factor(st.data$Niche)))
Idents(st.data) = "Niche"
saveRDS(st.data,file = "./Outdata/Step2.NEPC_ARPC_pca_Niche.rds")