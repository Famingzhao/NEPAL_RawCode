library(clustree)
library(cowplot)
library(stringr)
library(gplots)
library(rstatix)
library(ggsci)
library(ComplexUpset)
library(UpSetR)
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(future)
library(ComplexHeatmap)
library(circlize)
# library(SeuratDisk)
# library(SeuratWrappers)
# library(harmony)
# library(rliger)
# library(reshape2)

#### 1. load data
## 1.1 scRNA-seq
seurat.data = read_rds(file = "./Outdata/BBKNN.scRNAseq_data.rds")
DimPlot(seurat.data,raster=FALSE,label = T)
g.colSet = list(cluster.name.main = c(OE = "#3182BDFF", Luminal = "#E6550DFF", NE = "#31A354FF", 
                                      cECs = "#756BB1FF", aECs = "#636363FF", Fib = "#6BAED6FF", MyoFib = "#FD8D3CFF", 
                                      Mye = "#74C476FF", Mast = "#9E9AC8FF", CD4T = "#969696FF", Treg = "#9ECAE1FF", 
                                      CD8T = "#FDAE6BFF", NK = "#A1D99BFF", B = "#BCBDDCFF", `B plasma` = "#BDBDBDFF"),
                loc2 = c("Pri" = "#0072B5FF","CRPC" = "#20854EFF",
                         "mCRPC" = "#E18727FF","NEPC" = "#BC3C29FF",
                         "mHSPC" = "#84B8D7","mmHSPC" = "#6A3D9A"))
scales::show_col(g.colSet$loc2)

if(T){
  text.size = 7
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   #axis.line = element_line(color = "black"),
                   #axis.ticks = element_line(color = "black"),
                   #panel.grid.minor.y = element_blank(),
                   #panel.grid.minor.x = element_blank(),
                   panel.grid=element_blank(),
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size),
                   # panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
                   strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07", 
  )
}
#AverageExpression
aver.data = AverageExpression(seurat.data,group.by = "sampleID")
aver.data =  aver.data.sample$RNA %>% as.data.frame()
aver.meta = data.frame(sampleID=colnames(aver.data.sample))
aver.meta = left_join(aver.meta,select(seurat.data@meta.data,
                                       group,patientID,sampleID,age,data.sets:group))
aver.meta = aver.meta[!duplicated(aver.meta$sampleID),]

seurat.data.list = SplitObject(seurat.data,split.by = "sampleID")

aver.data.ct = lapply(seurat.data.list, function(x){
  as.data.frame(AverageExpression(x,group.by = "celltype.sub")$RNA)
})

## 1.2 gene list
gene.list = read.gmt2list("./gene.list/NE.cells.marker_Ten.gmt")
gene.list = gene.list$gs
input.data = fromList(gene.list)
inpu.gene = names(gene.list)

get.gene = function(DEG.list,n=4){
  tpm = unlist(DEG.list) %>% unique()
  
  data.tmp = matrix(0,nrow = length(tpm) ,ncol = length(DEG.list))
  row.names(data.tmp) = tpm
  colnames(data.tmp) = names(DEG.list)
  for (i in 1:length(DEG.list)) {
    data.tmp[,i] = ifelse(row.names(data.tmp) %in% DEG.list[[i]],1,0)
  }
  data.tmp= as.data.frame(data.tmp)
  data.tmp$sum = as.numeric(rowSums(data.tmp))
  
  row.names(data.tmp)[data.tmp$sum >= n]
}
length(get.gene(gene.list,n=1))
check_genes = get.gene(gene.list,n=4)
all.genes = get.gene(gene.list,n=1) 

#### 2. ComplexUpset Vis
if(F){
  set_size = function(w, h, factor=1.5) {
    s = 1 * factor
    options(
      repr.plot.width=w * s,
      repr.plot.height=h * s,
      repr.plot.res=100 / factor,
      jupyter.plot_mimetypes='image/png',
      jupyter.plot_scale=1
    )
  }
  text.size = 15
  ComplexUpset::upset(input.data, inpu.gene, width_ratio=0.1)+ggtitle("Meta-analysis of NE markers")&
    theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
          axis.ticks = element_line(color = "black"),
          # axis.title = element_text(size = text.size,color ="black"), 
          axis.text = element_text(size=text.size,color = "black"),
          # axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
          #axis.line = element_line(color = "black"),
          #axis.ticks = element_line(color = "black"),
          #panel.grid.minor.y = element_blank(),
          #panel.grid.minor.x = element_blank(),
          # panel.grid=element_blank(), # 去网格线
          # legend.position = legend.position,
          legend.text = element_text(size= text.size),
          legend.title= element_text(size= text.size))
}

#### 3. markers vis
## 3.1 EurUrol.2005
p1 = DotPlot(object = seurat.data, features = gene.list$EurUrol.2005,assay = "RNA",scale = T) + 
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "EurUrol.2005 by cell types")
p1

p2 = DotPlot(object = seurat.data,group.by = "group", 
             features = gene.list$EurUrol.2005,assay = "RNA",scale = T) + 
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "EurUrol.2005 by groups")
p2

## 3.2 CellRep.2018.scRNA
p3 = DotPlot(object = seurat.data, features = gene.list$CellRep.2018.scRNA,
             assay = "RNA",scale = T) + 
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "CellRep.2018.scRNA by cell types")
p3

p4 = DotPlot(object = seurat.data,group.by = "group", 
             features = gene.list$CellRep.2018.scRNA,assay = "RNA",scale = T) + 
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "CellRep.2018.scRNA by groups")
p4

## 3.3 NatMed.2016
p5 = DotPlot(object = seurat.data, features = gene.list$NatMed.2016,
             assay = "RNA",scale = T) + 
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "NatMed.2016 by cell types")
p5

p6 = DotPlot(object = seurat.data,group.by = "group", 
             features = gene.list$NatMed.2016,assay = "RNA",scale = T) + 
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "NatMed.2016 by groups")
p6

## 3.4 HP_NE_neoplasm
p7 = DotPlot(object = seurat.data, features = gene.list$HP_NE_neoplasm,
             assay = "RNA",scale = T) + 
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "HP_NE_neoplasm by cell types")
p7

p8 = DotPlot(object = seurat.data,group.by = "group", 
             features = gene.list$NatMed.2016,assay = "RNA",scale = T) + 
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "HP_NE_neoplasm by groups")
p8

## 3.5 meta analysis
check_genes = get.gene(gene.list,n=4)
check_genes = check_genes[check_genes%in%row.names(seurat.data)]

p9 = DotPlot(object = seurat.data, features = check_genes,
             assay = "RNA",scale = T) + theme_bw()+
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "Overlap >=4")
p9

test = filter(p9$data,id%in%"NE")
test = test[order(test$pct.exp),]
check_genes = test$features.plot[order(test$pct.exp)]

p9 = DotPlot(object = seurat.data, features = check_genes,
             assay = "RNA",scale = T) + theme_bw()+
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "Overlap >=4 (n=61)")
p9

p10 = DotPlot(object = seurat.data,group.by = "group", 
              features = check_genes,assay = "RNA",scale = T) + theme_bw()+
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "Overlap >=4 (n=61)")
p10

p.meta = wrap_plots(p9,p10,guides = "collect",widths = c(1,0.4))
p.meta
ggsave(p.meta,filename = "./Outplot/Step2.Overlap.4.Dot.plot.pdf",
       height = 28,width = 20,units = "cm")

## 3.5 scRNA-seq
p11 = DotPlot(object = seurat.data, features = NE.marker,
              assay = "RNA",scale = T) + 
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "scRNA-seq")
p11

test = filter(p11$data,id%in%"NE")
test = test[order(test$pct.exp),]
NE.marker = test$features.plot[order(test$pct.exp)]

p11 = DotPlot(object = seurat.data, features = NE.marker,dot.min = -10,
              assay = "RNA",scale = T) + 
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "scRNA-seq")
p11

p12 = DotPlot(object = seurat.data,group.by = "group", 
              features = NE.marker,dot.min = -10,
              assay = "RNA",scale = T) + 
  coord_flip() + mytheme + labs(x=NULL,y=NULL,title = "scRNA-seq")
p12

#### 4.heatmap
### 4.1 by cell types
aver.data.ct = AverageExpression(seurat.data)
aver.data.ct = aver.data.ct$RNA %>% as.data.frame()

input.data1 = aver.data.ct[row.names(aver.data.ct)%in%all.genes,levels(seurat.data)]
if(F){
  colnames(input.data1) = paste0(LETTERS[1:15],".",colnames(input.data1))
  
  # color
  col_fun = colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))
  
  # celltype
  anno_group = g.colSet$cluster.name.main
  names(anno_group) = colnames(input.data1)
  
  column_top = HeatmapAnnotation(
    group = colnames(input.data1),
    col = list(group = anno_group
    ),
    border = T)
  
  # left_annotation
  left_annotation =  rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                                    labels = c("UP in Epi",
                                                               "UP in Stroma",
                                                               "UP in Immune",
                                                               "UP in NE"), 
                                                    labels_gp = gpar(col = "white", 
                                                                     fontsize = 8)))
  # plot
  htdf1 <- t(scale(t(input.data1),scale = T,center = T))
  htkm1 = Heatmap(htdf1,
                  km = 4,
                  right_annotation = row_anno,
                  left_annotation = left_annotation,
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  column_title = "Meta−analysis of 1482 NE markers",
                  name = "Z-score",
                  cluster_columns = F,show_row_dend = F,
                  cluster_rows = T,show_row_names = F,
                  row_title = "Gene module",
                  row_names_gp = gpar(fontface = 'italic',fontsize = 10),
                  row_names_side = 'left',
                  border = T,
                  # rect_gp = gpar(col = "white", lwd = 1),
                  column_names_side = 'top',
                  column_names_rot = 45,
                  top_annotation = column_top,
                  column_split = colnames(input.data1),
                  col = col_fun
  )
  pdf(file = "./Outplot/Step2.NE_markers.heatmap.plot.celltypes.pdf",
      height = 6,width = 8)
  htkm1
  dev.off()
  
  str(row_order(htkm1))
  
  gene.module = lapply(row_order(htkm1), function(x){
    htkm1@row_names_param$labels[x]
  })
  str(gene.module)
  names(gene.module) = c("UP in Epi",
                         "UP in Stroma",
                         "UP in Immune",
                         "UP in NE")
  write.list2gmt(gene.module,file = "./Outdata/Step2.Meta.NE.markers.by.celltype.gmt")
 
  mark_gene = c("CHGA", "CHGB", "ENO2", "SYP", "CDKN2C", "SOX2", "EZH2", "NKX2-1", 
                "ASCL1", "SCG2", "SCGN", "TFF2", "TUBB2B", "FOXA2")
  mark_gene = mark_gene[mark_gene%in%rownames(input.data1)]
  gene_pos <- which(rownames(input.data1) %in% mark_gene)
  gene_lab <- rownames(input.data1)[gene_pos]
  
  row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,
                                                   labels = gene_lab))
}

### 4.2 htkm1+htkm2 by group
aver.data = AverageExpression(seurat.data,group.by = "group")
aver.data = aver.data$RNA %>% as.data.frame()

input.data = aver.data[as.character(unlist(gene.module)),levels(seurat.data$group)]
if(F){
  colnames(input.data) = paste0(LETTERS[1:4],".",colnames(input.data))
  
  # color
  col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#0099CC","white","#CC0033"))
  
  # celltype
  anno_group = g.colSet$loc2[levels(seurat.data$group)]
  names(anno_group) = colnames(input.data)
  
  column_top = HeatmapAnnotation(
    group = colnames(input.data),
    col = list(group = anno_group
    ),
    border = T)
  
  # left_annotation
  left_annotation =  rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                                    labels = c(
                                                      "UP in mCRPC",
                                                      "UP in NEPC",
                                                      "UP in Primary",
                                                      "UP in CRPC"), 
                                                    labels_gp = gpar(col = "white", 
                                                                     fontsize = 8)))
  # plot
  htdf <- t(scale(t(input.data),scale = T,center = T))
  htkm2 = Heatmap(htdf,
                  # km = 4,
                  # left_annotation = left_annotation,
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  column_title = "Meta−analysis of 1482 NE markers",
                  name = "Z-score",
                  cluster_columns = F,show_row_dend = F,
                  cluster_rows = F,show_row_names = F,
                  row_title = "Gene module",
                  row_names_gp = gpar(fontface = 'italic',fontsize = 10),
                  row_names_side = 'left',
                  border = T,
                  # right_annotation = row_anno,
                  # rect_gp = gpar(col = "white", lwd = 1),
                  column_names_side = 'top',
                  column_names_rot = 45,
                  top_annotation = column_top,
                  column_split = colnames(input.data),
                  col = col_fun
  )
  pdf(file = "./Outplot/Step2.NE_markers.heatmap.plot.pdf",height = 6,width = 10)
  htkm1 + htkm2
  dev.off()
}

### 4.3 by group
input.data = aver.data[row.names(aver.data)%in%
                         all.genes,levels(seurat.data$group)]
if(F){
  colnames(input.data) = paste0(LETTERS[1:4],".",colnames(input.data))
  
  # color
  col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#0099CC", "white", "#CC0033"))
  
  # celltype
  anno_group = g.colSet$loc2[levels(seurat.data$group)]
  names(anno_group) = colnames(input.data)
  
  column_top = HeatmapAnnotation(
    group = colnames(input.data),
    col = list(group = anno_group
    ),
    border = T)
  
  # left_annotation
  left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                                   labels = c(
                                                     "UP in NEPC",
                                                     "UP in mCRPC",
                                                     "UP in Primary",
                                                     "UP in CRPC"), 
                                                   labels_gp = gpar(col = "white", 
                                                                    fontsize = 8)))
  # plot
  htdf <- t(scale(t(input.data),scale = T,center = T))
  htkm3 = Heatmap(htdf,
                  km = 4,
                  right_annotation = row_anno,
                  left_annotation = left_annotation,
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  column_title = "Meta−analysis of 1482 NE markers",
                  name = "Z-score",
                  cluster_columns = F,show_row_dend = F,
                  cluster_rows = T,show_row_names = F,
                  row_title = "Gene module",
                  row_names_gp = gpar(fontface = 'italic',fontsize = 10),
                  row_names_side = 'left',
                  border = T,
                  # rect_gp = gpar(col = "white", lwd = 1),
                  column_names_side = 'top',
                  column_names_rot = 45,
                  top_annotation = column_top,
                  column_split = colnames(input.data),
                  col = col_fun
  )
  pdf(file = "./Outplot/Step2.NE_markers.heatmap.plot.by.Group.pdf",height = 6,width = 6)
  htkm3
  dev.off()
  
  row_order(htkm3)
  gene.module3 = lapply(row_order(htkm3), function(x){
    htkm3@row_names_param$labels[x]
  })
  str(gene.module3)
  names(gene.module3) = c(
    "UP in NEPC",
    "UP in mCRPC",
    "UP in Primary",
    "UP in CRPC")
  write.list2gmt(gene.module3,file = "./Outdata/Step2.Meta.NE.markers.by.Group.gmt")
  
  mark_gene = c("CHGA", "CHGB", "ENO2", "SYP", "CDKN2C", "SOX2", "EZH2", "NKX2-1", 
                "ASCL1", "SCG2", "SCGN", "TFF2", "TUBB2B", "FOXA2")
  mark_gene = mark_gene[mark_gene%in%rownames(input.data)]
  gene_pos <- which(rownames(input.data) %in% mark_gene)
  gene_lab <- rownames(input.data)[gene_pos]
  
  row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,
                                                   labels = gene_lab))
}

#### 5. ssGSEA scores
(load("Outdata/Step4.SU2C_WCDT_UMSPORE_PCaProfilter_scRNAseq.NE.Rdata"))
library(pROC)
val_dd_list = lapply(names(NE.expr.list), function(x){
  data.sets = data.frame(ID = colnames(NE.expr.list[[x]]))
  data.sets$NE_status = ifelse(NE.phe.list[[x]]$NE_status %in% c("CRPC_NE","NEPC"),
                               1,0)
  data.sets = select(data.sets,"ID","NE_status",everything())
  return(data.sets)
})
names(val_dd_list) = names(NE.expr.list)

CstatisticCI <- function(x) {
  se <- x["S.D."]/sqrt(x["n"])
  Low95 <- x["C Index"] - 1.96*se 
  Upper95 <- x["C Index"] + 1.96*se 
  cbind(x["C Index"], Low95, Upper95) 
}
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

gene.list = read.gmt2list("./gene.list/NE.cells.marker_Ten.gmt")
gene.list = gene.list$gs

### 5.1 Calculate index
for (i in names(NE.expr.list)) {
  pca.data = NE.expr.list[[i]]
  phe.data = val_dd_list[[i]]

  sig_tme <- calculate_sig_score(pdata           = NULL,
                                 eset            = pca.data,
                                 signature       = gene.list,
                                 method          = "ssgsea",
                                 mini_gene_count = 0)
  sig_tme = arrange(sig_tme,Index)
  phe.data = cbind(phe.data,sig_tme[,-c(1,2)])
  val_dd_list[[i]] = phe.data
}
str(val_dd_list)

res.meta = data.frame()
for (i in names(gene.list)) {
  fml = as.formula(paste0("NE_status~",i))
  cc <- data.frame(
    R2=sapply(val_dd_list,function(x){as.numeric(R2(x[,i], x$NE_status))}),
    RMSE=sapply(val_dd_list,function(x){as.numeric(RMSE(x[,i], x$NE_status))}),
    MSE=sapply(val_dd_list,function(x){as.numeric(mean((x[,i]-x$NE_status)^2))}),
    MAE=sapply(val_dd_list,function(x){as.numeric(MAE(x[,i], x$NE_status))}),
    Cindex=sapply(val_dd_list,function(x){as.numeric(CstatisticCI(rcorr.cens(x[,i], x$NE_status))[1])}),
    AUC = sapply(val_dd_list,function(x){round(roc(fml,data=x,aur=TRUE,ci=TRUE)$auc,2)})
  )%>%
    rownames_to_column('ID')
  cc$Gene_sets = i
  res.meta = rbind(res.meta,cc)
}

res.meta_mean <- res.meta%>%
  # filter(ID!='SU2C')%>%
  dplyr::group_by(Gene_sets)%>%
  dplyr::mutate(AUC_mean=mean(AUC,na.rm =T))%>%
  dplyr::mutate(Cindex_mean=mean(Cindex,na.rm =T))%>%
  dplyr::mutate(R2_mean=mean(R2,na.rm =T))%>%
  dplyr::mutate(RMSE_mean=mean(RMSE,na.rm =T))

res.meta %>%
  # filter(ID!='SU2C')%>%
  ggplot(aes(Cindex,reorder(Gene_sets,Cindex)))+
  geom_bar(width = 0.7,stat = 'summary',fun='mean',fill='#0072B5FF')+
  stat_summary(fun = mean,
               geom = "pointrange",size=0.08,col="black",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x))+
  theme_classic()+
  labs(y=NULL)

### 5.2 Vis
res<-roc(NE_status~EurUrol.2005,data=val_dd_list$SU2C,aur=TRUE,
         ci=TRUE,
         smooth=F
)
ggroc(res, color ="red",legacy.axes = TRUE)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw() + 
  ggtitle("PUM1 ROC Curve")+
  theme(plot.title = element_text(hjust = 0.5,size = 16),
        axis.text=element_text(size=12,colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  annotate("text",x=0.75,y=0.25,label=paste("AUC = ", round(res$auc,3)))+
  annotate("text",x=0.75,y=0.20,label=paste("95%CI: ", round(res$ci[1],3),'-',round(res$ci[3],3)))

### scRNA-seq vs. other scores
(load(file = "./Outdata/Step2.Meta.NE_markers_ssGSEA_score.Rdata"))
input.data.scRNA = fre.data[match(val_dd_list$scRNAseq$ID,fre.data$SampleID),]
input.data.scRNA = cbind(input.data.scRNA,dplyr::select(val_dd_list$scRNAseq,3:13))
cor.test(input.data.scRNA$freq,input.data.scRNA$EurUrol.2005)
cor.test(input.data.scRNA$freq,input.data.scRNA$JClinInvest.2019)

plot.data = dplyr::select(input.data.scRNA,freq:BMC.Cancer.2017) %>% t() %>% as.data.frame()
plot.data = batch_cor(gene = "freq",
                      exprSet = plot.data,
                      rownames = row.names(plot.data),
                      method = "p")
plot.data = plot.data[order(plot.data$cor),]
plot.data$mRNAs = factor(plot.data$mRNAs,levels = unique(plot.data$mRNAs))
plot.data$FDR = -log(plot.data$p.value)
range(plot.data$FDR)

draw_stick =  ggplot(plot.data,aes(cor,mRNAs)) +
  geom_segment(aes(xend=0,yend=mRNAs)) +
  geom_point(aes(col=FDR,size=abs(cor))) +
  scale_colour_gradientn(colours=c("#7fc97f","#984ea3"),
                         name = "-log10 (p-value)",limit=c(0,51)) +
  #scale_color_viridis_c(begin = 0.5, end = 1) +
  scale_size_continuous(range =c(0,4))  +
  theme_bw() + mytheme + labs(y=NULL,x="Cor",title = "Cellular fraction vs. Predict scores")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
draw_stick
