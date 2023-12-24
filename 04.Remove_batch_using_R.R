rm(list=ls())
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(future)
library(clustree)
library(cowplot)
library(stringr)
library(harmony)
# change the current plan to access parallelization
plan("multisession", workers = 4)
plan()
options(future.globals.maxSize = 10 * 1024^3)

##### 1.input data
seurat.data = read_rds(file = "./Outdata/Step1.seurat.data_after_QC.rds")

##### 2.Standard pipeline
seurat.data <- seurat.data %>% NormalizeData(verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
  ScaleData(verbose = F) %>%
  RunPCA(npcs = 30, verbose = F)

seurat.data <- seurat.data %>% 
  RunUMAP(reduction = "pca", dims = 1:30, verbose = F) 

p1.compare=wrap_plots(ncol = 2,
                      DimPlot(seurat.data, reduction = "pca", group.by = "sampleID")+NoAxes()+ggtitle("Before_PCA"),
                      DimPlot(seurat.data, reduction = "umap", group.by = "sampleID")+NoAxes()+ggtitle("Before_UMAP"),
                      guides = "collect"
)
p1.compare

##### 3. Remove batch using R
scRemoveBatch = function(seurat.data, 
                         batchID = "Data.sets", 
                         methods = "harmony", n.pcs = 30){
  
  library(Seurat)
  cat("##################################\n")
  cat("R code from FamingZhao.\n") 
  cat("Welcome to cite our article <Integrated analysis of single-cell and bulk transcriptomics develops a robust neuroendocrine cell-intrinsic signature to predict prostate cancer progression>.\n")
  cat("##################################\n")
  if(methods %in% c("none")){
    seurat_int <- seurat.data %>% NormalizeData(verbose = F) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
      ScaleData(verbose = F) %>%
      RunPCA(npcs = n.pcs, verbose = F)%>% 
      RunUMAP(reduction = "pca", dims = 1:20, verbose = F)
  }
  
  ### 1. harmony
  if(methods %in% c("harmony")){
    library(harmony)
    seurat.data <- seurat.data %>% NormalizeData(verbose = F) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
      ScaleData(verbose = F) %>%
      RunPCA(npcs = n.pcs, verbose = F)
    
    seurat.data <- seurat.data %>% RunHarmony(batchID, plot_convergence = T)
    
    seurat_int <- seurat.data %>% 
      RunUMAP(reduction = "harmony", dims = 1:20, verbose = F)
  }
  
  ### 2. CCA
  if(methods %in% c("CCA")){
    seurat_list = SplitObject(seurat.data, split.by = batchID)
    rm(seurat.data)
    
    # normalize and identify variable features for each dataset independently
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      return(x)
    })
    
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = seurat_list)
    
    #anchors
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                             anchor.features = features, dims = 1:30)
    rm(seurat_list)
    
    seurat_int <- IntegrateData(anchorset = seurat_anchors, dims = 1:n.pcs) %>% 
      ScaleData() %>% 
      RunPCA(npcs = n.pcs, verbose = FALSE) %>% 
      RunUMAP(reduction = "pca", dims = 1:n.pcs)
  }
  
  ### 3. RPCA
  if(methods %in% c("RPCA")){
    seurat_list = SplitObject(seurat.data, split.by = batchID)
    rm(seurat.data)
    
    # normalize and identify variable features for each dataset independently
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      return(x)
    })
    
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = seurat_list)
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE) %>% 
        RunPCA( features = features, verbose = FALSE)
      return(x)
    })
    
    #Perform integration
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                             anchor.features = features, 
                                             reduction = "rpca", 
                                             k.anchor = 20 # 5 by default. Increasing this parameter to 20 will assist in aligning these populations.
    )
    rm(seurat_list)
    
    seurat_int <- IntegrateData(anchorset = seurat_anchors, dims = 1:n.pcs) %>% 
      ScaleData() %>% 
      RunPCA(npcs = n.pcs, verbose = FALSE) %>% 
      RunUMAP(reduction = "pca", dims = 1:n.pcs)
  }
  
  ### 4. fastMNN
  if(methods %in% c("fastMNN")){
    library(SeuratWrappers)
    seurat_list = SplitObject(seurat.data, split.by = batchID)
    rm(seurat.data)
    
    seurat_list <- lapply(seurat_list, FUN = function(x){
      x = x %>% NormalizeData() %>% 
        FindVariableFeatures()
      return(x)
    })
    
    seurat_int <- RunFastMNN(object.list = seurat_list) %>% 
      RunUMAP(reduction = "mnn", dims = 1:n.pcs)
  }
  
  ### 5. LIGER
  if(methods %in% c("LIGER")){
    library(rliger)
    seurat_int <- seurat.data %>% NormalizeData(verbose = F) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
      ScaleData(verbose = F, do.center = F) %>%
      RunPCA(npcs = n.pcs, verbose = F)
    rm(seurat.data)
    
    seurat_int <- RunOptimizeALS(seurat_int, k = 30, lambda = 5, split.by = batchID) %>% 
      RunQuantileNorm(split.by = batchID)
    
    n.pcs = ncol(seurat_int[["iNMF"]])
    seurat_int <- seurat_int %>% 
      RunUMAP(reduction = "iNMF", dims = 1:n.pcs, verbose = F)
  }
  
  ### 6. conos
  if(methods %in% c("conos")){
    library(conos)
    seurat_list = SplitObject(seurat.data, split.by = batchID)
    rm(seurat.data)
    gc()
    
    panel.preprocessed.seurat <- lapply(seurat_list, function(x){
      x <- NormalizeData(x) %>% 
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
        ScaleData(verbose = FALSE) %>% 
        RunPCA(npcs = n.pcs, verbose = FALSE)
      return(x)
    })
    
    #### 2 Conos
    con <- Conos$new(panel.preprocessed.seurat, n.cores=2)
    space='PCA' # PCA, CPCA,CCA
    
    con$buildGraph(k=30, k.self=5, 
                   space=space,  
                   ncomps=30, 
                   n.odgenes=2000, 
                   matching.method='mNN', 
                   metric='angular', 
                   score.component.variance=TRUE, 
                   verbose=TRUE)
    plotComponentVariance(con, space=space)  
    
    con$findCommunities(method=leiden.community, resolution=1)
    con$embedGraph()
    con$plotGraph(alpha=0.1)
    con$plotPanel(font.size=4)
    plotClusterBarplots(con, legend.height = 0.1)
    seurat_int = con
  }
  
  return(seurat_int)
}

seurat.harmony <- scRemoveBatch(seurat.data = seurat.data, methods = "harmony")
seurat.CCA <- scRemoveBatch(seurat.data = seurat.data, methods = "CCA")
seurat.RPCA <- scRemoveBatch(seurat.data = seurat.data, methods = "RPCA")
seurat.fastMNN <- scRemoveBatch(seurat.data = seurat.data, methods = "fastMNN")
seurat.LIGER <- scRemoveBatch(seurat.data = seurat.data, methods = "LIGER")
conos.data <- scRemoveBatch(seurat.data = seurat.data, methods = "conos")

##### 4. Remove batch using Python, please see python scripts.
