library(Seurat)
library(dplyr)
library(future)
library(Matrix)
library(ggplot2)
library(readr)

# change the current plan to access parallelization
plan("multisession", workers = 4)
plan()
options(future.globals.maxSize = 100 * 1024^3)

#### 1. input data
seurat.data <- read_rds(file = '/home/data/ssy43/PCa/all.cells.rds')

#### 2. QC
### 2.1 Doublet Linux
batch.data = read.csv("./Outdata/adata_scrublet.csv",header=T)
table(batch.data$predicted_doublet)
meta.data = seurat.data@meta.data
meta.data = left_join(meta.data, batch.data)
row.names(meta.data) = row.names(seurat.data@meta.data)
seurat.data@meta.data = meta.data
table(seurat.data$predicted_doublet)

### 2.2 QC index calculation
sc.QC.caculate = function(seurat.data = seurat.data,org = "human",cellcycle = F){
  cat("##################################\n")
  cat("R code from FamingZhao.\n") 
  cat("Welcome to cite our article <PMID: 38250042.Integrated analysis of single-cell and bulk transcriptomics develops a robust neuroendocrine cell-intrinsic signature to predict prostate cancer progression>.\n")
  cat("##################################\n")
  ###Runs Seurat for QC calculate
  ##Inputs:
  #seurat.data = Seruat object
  #org = "human" or "mouse"
  
  ##Outputs:
  #seurat.data with
  #log10GenesPerUMI
  #Mitochondrial ratio, percent.mt
  #Ribosomal ratio, percent.ribo
  #Erythrocyte ratio, percent.hb
  #Cell cycle, S.Score and G2M.Score
  
  if(org == "human"){
    ### 1.1 Number of genes detected per UMI
    # Add number of genes per UMI for each cell to metadata
    seurat.data$log10GenesPerUMI <- log10(seurat.data$nFeature_RNA) / log10(seurat.data$nCount_RNA)
    
    ### 1.2 QC index calculation
    ## 1.2.1 Mitochondrial ratio
    mito_genes=rownames(seurat.data)[grep("^MT-", rownames(seurat.data))]
    print("mito_genes:")
    cat(mito_genes)
    seurat.data <- PercentageFeatureSet(object = seurat.data, pattern = "^MT-",col.name = "percent.mt")
    
    ## 1.2.2 Ribosomal ratio
    ribo_genes=rownames(seurat.data)[grep("^RP[SL]", rownames(seurat.data),ignore.case = T)]
    print("ribo_genes:")
    cat(ribo_genes)
    
    seurat.data <- PercentageFeatureSet(seurat.data,"^RP[SL]",col.name = "percent.ribo")
    print(summary(seurat.data@meta.data$percent.ribo))
    
    ## 1.2.3 Erythrocyte ratio
    hb_genes <- rownames(seurat.data)[grep("^HB[^(P)]", rownames(seurat.data),ignore.case = T)]
    print("hb_genes")
    cat(hb_genes)
    seurat.data=PercentageFeatureSet(seurat.data, "^HB[^(P)]", col.name = "percent.hb")
    print(summary(seurat.data@meta.data$percent.hb))
    
    if(cellcycle){
      ## 1.2.4 Cell cycle
      # human cell cycle genes
      s.genes=Seurat::cc.genes.updated.2019$s.genes
      print(s.genes)
      g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
      
      seurat.data=CellCycleScoring(object = seurat.data, 
                                     s.features = s.genes, 
                                     g2m.features = g2m.genes, 
                                     set.ident = TRUE)
    }
  }
  if(org == "mouse"){
    ### 1.1 Number of genes detected per UMI
    # Add number of genes per UMI for each cell to metadata
    seurat.data$log10GenesPerUMI <- log10(seurat.data$nFeature_RNA) / log10(seurat.data$nCount_RNA)
    
    ### 1.2 QC index calculation
    ## 1.2.1 Mitochondrial Ratio
    mito_genes=rownames(seurat.data)[grep("^mt-", rownames(seurat.data))]
    print("mito_genes:")
    cat(mito_genes)
    seurat.data <- PercentageFeatureSet(object = seurat.data, pattern = "^mt-",col.name = "percent.mt")
    
    ## 1.2.2 Ribosomal ratio
    ribo_genes=rownames(seurat.data)[grep("^Rp[sl]", rownames(seurat.data),ignore.case = T)]
    print("ribo_genes:")
    cat(ribo_genes)
    
    seurat.data <- PercentageFeatureSet(seurat.data,"^Rp[sl]",col.name = "percent.ribo")
    print(summary(seurat.data@meta.data$percent.ribo))
    
    ## 1.2.3 Erythrocyte ratio
    hb_genes <- rownames(seurat.data)[grep("^Hb[^(p)]", rownames(seurat.data),ignore.case = T)]
    print("hb_genes")
    cat(hb_genes)
    seurat.data=PercentageFeatureSet(seurat.data, "^Hb[^(p)]", col.name = "percent.hb")
    print(summary(seurat.data@meta.data$percent.hb))
    
    if(cellcycle){
      ## 1.2.4 Cell cycle
      load("~/ref_annotation_Geneset/7.CellCycleRdata/scRNAseq_cc.genes.updated.2019_human_mouse.Rdata")
      s.genes=mus_s.genes
      g2m.genes=mus_g2m.genes
      
      seurat.data=CellCycleScoring(object = seurat.data, 
                                     s.features = s.genes, 
                                     g2m.features = g2m.genes, 
                                     set.ident = TRUE)
    }
  }
  return(seurat.data)
}
seurat.data = sc.QC.caculate(seurat.data)

### 2.3 Vis
p.VlnPlot = VlnPlot(seurat.data,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"),
                    ncol = 5,
                    group.by = "sampleID")
p.VlnPlot
ggsave(filename = "./Outplot/01.QC.VlnPlot.pdf",width = 12,height = 4)

### 2.4 fiter
batch.data = read.csv("./Outdata/adata_scrublet.csv",header=T)
table(batch.data$predicted_doublet)

meta.data = seurat.data@meta.data
meta.data = left_join(meta.data, batch.data)
row.names(meta.data) = row.names(seurat.data@meta.data)
seurat.data@meta.data = meta.data
table(seurat.data$predicted_doublet)

seurat.qc <- subset(seurat.data, subset = nFeature_RNA >= 250 & nCount_RNA>= 500 & percent.mt < 20 & predicted_doublet =="False")

### 2.5 high dropout genes
# Extract counts
count.data <- GetAssayData(object = seurat.qc, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- count.data > 0

# Sums all TRUE values and returns TRUE if more than 500 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 500 
table(keep_genes)
filtered_counts <- count.data[keep_genes,]

# Reassign to filtered Seurat object MT+RP+AC
keep.genes = row.names(filtered_counts)[grep("^MT-|^RP[SL]|AC[0-2]", 
                                             row.names(filtered_counts),
                                             invert = T,value = F)]
keep.genes

filtered_counts = filtered_counts[keep.genes,]
dim(filtered_counts)
seurat.qc <- CreateSeuratObject(filtered_counts, 
                                meta.data = seurat.qc@meta.data,
                                min.cells = 0, min.features = 0)

### 2.6 save
saveRDS(seurat.qc, file = "./Outdata/Step1.seurat.data_after_QC.rds")
