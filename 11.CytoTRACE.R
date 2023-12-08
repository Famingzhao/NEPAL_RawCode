rm(list=ls())
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(clustree)
library(cowplot)
library(stringr)
library(CytoTRACE)

#### 1.input data
NE.data = read_rds("./Outdata/Step1.NE.cells.rds")
meta.data = NE.data@meta.data

#### 2. Run iCytoTRACE
expr.data = as.my.matrix(GetAssayData(object = NE.data[["RNA"]],slot = "counts"))
results <- CytoTRACE(mat = expr.data,ncores = 10)
save(results,file= "./Output_Cyto/Step11.results.CytoTRACE.Rdata")

plotCytoGenes(results, numOfGenes = 10,outputDir = "Output_Cyto/")

phe.data = as.character(meta.data$celltype.NE)
names(phe.data) = row.names(meta.data)
plotCytoTRACE(results, phenotype = phe.data,
              emb = as.matrix(NE.data@reductions$bbknn.umap@cell.embeddings),
              gene = "TOP2A",outputDir = "Output_Cyto/")

CytoTRACE.score = read.table("./Output_Cyto/CytoTRACE_plot_table.txt",header = T,row.names = 1)
head(CytoTRACE.score)
CytoTRACE.score = CytoTRACE.score[match(colnames(NE.data),row.names(CytoTRACE.score)),]
NE.data$CytoTRACE = CytoTRACE.score$CytoTRACE
