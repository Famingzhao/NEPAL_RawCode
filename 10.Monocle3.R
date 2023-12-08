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
library(monocle3)
plan("multisession", workers = 5)
plan()
options(future.globals.maxSize = 10 * 1024^3)

out.prefix = "./Outplot/"

#### 1.input data
NE.data = read_rds("./Outdata/Step1.NE.cells.rds")

#monocle3
### 1.1 Make the CDS object
data <- GetAssayData(NE.data, assay = 'RNA', slot = 'counts')
cell_metadata <- NE.data@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

### 1.2 Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)

## check batch
# plot_cells(cds, color_cells_by="group", label_cell_groups=FALSE,reduction_method = "PCA")
plot_cells(cds, color_cells_by="celltype.NE", 
           label_cell_groups=FALSE,reduction_method = "PCA")
# plot_cells(cds, color_cells_by="data.sets", label_cell_groups=FALSE,reduction_method = "PCA")

#### Step 2: Remove batch effects with cell alignment
# cds <- align_cds(cds, alignment_group = "sampleID", num_dim = 100)
# plot_cells(cds, color_cells_by="sampleID", label_cell_groups=FALSE,reduction_method = "Aligned")

#### Step 3: Reduce the dimensions using UMAP
# cds <- reduce_dimension(cds, preprocess_method = "Aligned")
# p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="data.sets") + ggtitle('cds.umap')
# p1

#### Step 4: umap from Seurat data
cds.embed <- cds@int_colData$reducedDims$PCA
int.embed <- Embeddings(NE.data, reduction = "bbknn.umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype.NE") +
  plot_cells(cds, reduction_method="UMAP", color_cells_by="group2") 
p2

#### Step 5: Monocle3 Cluster
cds <- cluster_cells(cds, resolution=0.005)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2);p

#### Step 6: Identification trajectory
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "celltype.NE",
           label_groups_by_cluster=T,
           label_leaves=T,
           label_branch_points=T)

#### Step 7: Order
# check# order_cells
FeaturePlot(NE.data,features = "NE_score",
            reduction = "bbknn.umap",pt.size = 0.1)&NoAxes()

cds <- order_cells(cds) 

temp <- RColorBrewer::brewer.pal(11, "Spectral")
temp[6] <- "gold"
rbPal <- colorRampPalette(temp)
(plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
            label_leaves = FALSE,  label_branch_points = F,label_roots = F,
            show_trajectory_graph = T)&NoAxes()&ggtitle("Monocle3"))+
(FeaturePlot(NE.data,features = "NE_score",
             reduction = "bbknn.umap",pt.size = 0.1)&NoAxes())+
  (FeaturePlot(NE.data,features = "CytoTRACE",
               reduction = "bbknn.umap",pt.size = 0.1)&NoAxes()&
     ggplot2::scale_colour_gradientn(name = "Predicted\norder", 
                                     colours = rev(rbPal(50)), 
                                     guide = ggplot2::guide_colourbar(ticks.colour = "black", 
                                                                      ticks.linewidth = 1, frame.colour = "black"), 
                                     breaks = seq(0, 1, 0.2), labels = c("0.0 (More diff.)", 
                                                                         0.2, 0.4, 0.6, 0.8, "1.0 (Less diff.)")) )