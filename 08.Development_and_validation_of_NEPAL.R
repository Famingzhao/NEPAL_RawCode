library(survival)
library(randomForestSRC)
library(glmnet)
# install.packages("plsRcox")
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(StabRvsPair)
library(stringr)
library(rms)
CstatisticCI <- function(x) {
  se <- x["S.D."]/sqrt(x["n"])
  Low95 <- x["C Index"] - 1.96*se 
  Upper95 <- x["C Index"] + 1.96*se 
  cbind(x["C Index"], Low95, Upper95) 
}

#### 1. traning
input.data = NE.expr.list$SU2C[cand.genes,]
phe.data = NE.phe.list$SU2C
val_dd_list = lapply(names(NE.expr.list), function(x){
  data.sets = NE.expr.list[[x]]
  data.sets = data.sets[cand.genes,]
  row.names(data.sets) = cand.genes
  data.sets[is.na(data.sets)] = 0
  
  data.sets <- t(data.sets) %>% as.data.frame()
  data.sets$ID = row.names(data.sets)
  data.sets$NE_status = NE.phe.list[[x]]$NE_status
  data.sets = select(data.sets,"ID","NE_status",everything())
  return(data.sets)
})
names(val_dd_list) = names(NE.expr.list)

result.data <- data.frame()
est_data <- t(input.data) %>% as.data.frame()
est_data$ID = row.names(est_data)
row.names(est_data) ==  NE.phe.list$SU2C$SAMPLE_ID

est_data$NE_status = ifelse(NE.phe.list$SU2C$NE_status == "CRPC_NE",1,0)

table(est_data$NE_status)

est_data = select(est_data,"ID","NE_status",everything())
pre_var <- colnames(est_data)[-c(1:2)]
est_dd <- est_data[,c('NE_status',pre_var)]
dim(est_dd)

rf_nodesize <- 5
seed <- 1

##################################
#### 1.RSF ####
##################################
set.seed(seed)
fit <- rfsrc(NE_status~.,data = est_dd,
             ntree = 1000,
             nodesize = rf_nodesize,
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs.RFS <- lapply(val_dd_list, function(x){cbind(x[,1:2],
                                                RS = predict(RFS.fit, newdata = x)$predicted)})

cc <- data.frame(
  R2=sapply(rs.RFS,function(x){as.numeric(R2(x$RS, x$NE_status))}),
  RMSE=sapply(rs.RFS,function(x){as.numeric(RMSE(x$RS, x$NE_status))}),
  MSE=sapply(rs.RFS,function(x){as.numeric(mean((x$RS-x$NE_status)^2))}),
  MAE=sapply(rs.RFS,function(x){as.numeric(MAE(x$RS, x$NE_status))}),
  Cindex=sapply(rs.RFS,function(x){as.numeric(CstatisticCI(rcorr.cens(x$RS, x$NE_status))[1])}),
  AUC = sapply(rs.RFS,function(x){round(roc(NE_status~RS,data=x,aur=TRUE,ci=TRUE)$auc,2)})
)%>%
  rownames_to_column('ID')

cc$Model <- 'RSF'
result.data <- rbind(result.data, cc)

rid <- var.select(object = rs.RFS, conservative = "high")
feature.gene <- rid$topvars
length(feature.gene)

##################################
#### 2.Enet ####
##################################
x1 <- as.matrix(est_dd[,pre_var])
x2 <- est_dd$NE_status
Enet.res.list = list()
for (alpha in c(0.01, seq(0.01,1,0.1))) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "binomial",alpha=alpha,nfolds = 10)
  
  rs <- lapply(val_dd_list, function(x){
    cbind(x[,1:2],
          RS = as.numeric(predict(fit,type='link',
                                  newx=as.matrix(x[,-c(1,2)]),
                                  s=fit$lambda.min)))})
  Enet.res.list[[paste0("alpha.",alpha)]] = fit
  
  cc <- data.frame(
    R2=sapply(rs,function(x){as.numeric(R2(x$RS, x$NE_status))}),
    RMSE=sapply(rs,function(x){as.numeric(RMSE(x$RS, x$NE_status))}),
    MSE=sapply(rs,function(x){as.numeric(mean((x$RS-x$NE_status)^2))}),
    MAE=sapply(rs,function(x){as.numeric(MAE(x$RS, x$NE_status))}),
    Cindex=sapply(rs,function(x){as.numeric(CstatisticCI(rcorr.cens(x$RS, x$NE_status))[1])}),
    AUC = sapply(rs,function(x){round(roc(NE_status~RS,data=x,aur=TRUE,ci=TRUE)$auc,2)})
  )%>%
    rownames_to_column('ID')
  
  cc$Model <- paste0('Enet','[α=',alpha,']')
  result.data <- rbind(result.data,cc)
}

##################################
#### 3. Ridge ####
##################################
x1 <- as.matrix(est_dd[,pre_var])
x2 <- est_dd$NE_status
set.seed(seed)

fit = cv.glmnet(x1, x2,family = "binomial",alpha=0,nfolds = 10, lambda = NULL)
cvfit = cv.glmnet(x1, x2,
                  nfold = 10, 
                  family = "binomial",
                  type.measure = "class"
)

rs <- lapply(val_dd_list, function(x){
  cbind(x[,1:2],
        RS = as.numeric(predict(fit,type='response',
                                newx=as.matrix(x[,-c(1,2)]),
                                s=fit$lambda.min)))})

cc <- data.frame(
  R2=sapply(rs,function(x){as.numeric(R2(x$RS, x$NE_status))}),
  RMSE=sapply(rs,function(x){as.numeric(RMSE(x$RS, x$NE_status))}),
  MSE=sapply(rs,function(x){as.numeric(mean((x$RS-x$NE_status)^2))}),
  MAE=sapply(rs,function(x){as.numeric(MAE(x$RS, x$NE_status))}),
  Cindex=sapply(rs,function(x){as.numeric(CstatisticCI(rcorr.cens(x$RS, x$NE_status))[1])}),
  AUC = sapply(rs,function(x){round(roc(NE_status~RS,data=x,aur=TRUE,ci=TRUE)$auc,2)})
)%>%
  rownames_to_column('ID')

cc$Model <- "Ridge"
result.data <- rbind(result.data,cc)

Ridge.fit = fit
Ridge.cvfit =  cvfit

##################################
#### 4.superpc####
##################################
data <- list(x=t(est_dd[,-c(1)]),
             y=est_dd$NE_status,
             featurenames=colnames(est_dd)[-c(1)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'regression',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(val_dd_list,function(w){
  test <- list(x=t(w[,-c(1,2)]),
               y=w$NE_status,
               featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,
                        data,
                        test,
                        threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],
                        n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})

cc <- data.frame(
  R2=sapply(rs,function(x){as.numeric(R2(x$RS, x$NE_status))}),
  RMSE=sapply(rs,function(x){as.numeric(RMSE(x$RS, x$NE_status))}),
  MSE=sapply(rs,function(x){as.numeric(mean((x$RS-x$NE_status)^2))}),
  MAE=sapply(rs,function(x){as.numeric(MAE(x$RS, x$NE_status))}),
  Cindex=sapply(rs,function(x){as.numeric(CstatisticCI(rcorr.cens(x$RS, x$NE_status))[1])}),
  AUC = sapply(rs,function(x){round(roc(NE_status~RS,data=x,aur=TRUE,ci=TRUE)$auc,2)})
)%>%
  rownames_to_column('ID')

cc$Model <- "superpc"
result.data <- rbind(result.data,cc)
superpc.cv.fit = cv.fit
superpc.fit = fit
superpc.data = data

##################################
#### 5.GBM ####
##################################
set.seed(seed)
fit <- gbm(formula = NE_status~.,data = est_dd,
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 10)

# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = NE_status~.,data = est_dd,distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 10)
GBM.fit = fit

rs <- lapply(val_dd_list, function(x){
  cbind(x[,1:2],
        RS = as.numeric(predict(GBM.fit,x,n.trees = best,type = 'link')))})
cc <- data.frame(
  R2=sapply(rs,function(x){as.numeric(R2(x$RS, x$NE_status))}),
  RMSE=sapply(rs,function(x){as.numeric(RMSE(x$RS, x$NE_status))}),
  MSE=sapply(rs,function(x){as.numeric(mean((x$RS-x$NE_status)^2))}),
  MAE=sapply(rs,function(x){as.numeric(MAE(x$RS, x$NE_status))}),
  Cindex=sapply(rs,function(x){as.numeric(CstatisticCI(rcorr.cens(x$RS, x$NE_status))[1])}),
  AUC = sapply(rs,function(x){round(roc(NE_status~RS,data=x,aur=TRUE,ci=TRUE)$auc,2)})
)%>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result,cc)

##################################
#### 6.SVM ####
##################################
library(e1071)
fit.svm <- svm(NE_status~., 
               data= est_dd,
               kernel = "polynomial", 
               gamma= 0.1, 
               cost= 1)
save(fit.svm,file = "./Step4.machine_learning/fit.svm.Rdata")

rs <- lapply(val_dd_list, function(x){
  cbind(x[,1:2],
        RS = as.numeric(predict(fit.svm,x)) %>% normalize())})
cc <- data.frame(
  R2=sapply(rs,function(x){as.numeric(R2(x$RS, x$NE_status))}),
  RMSE=sapply(rs,function(x){as.numeric(RMSE(x$RS, x$NE_status))}),
  MSE=sapply(rs,function(x){as.numeric(mean((x$RS-x$NE_status)^2))}),
  MAE=sapply(rs,function(x){as.numeric(MAE(x$RS, x$NE_status))}),
  Cindex=sapply(rs,function(x){as.numeric(CstatisticCI(rcorr.cens(x$RS, x$NE_status))[1])}),
  AUC = sapply(rs,function(x){round(roc(NE_status~RS,data=x,aur=TRUE,ci=TRUE)$auc,2)})
)%>%
  rownames_to_column('ID')
cc$Model <- paste0('SVM')
result.data <- rbind(result.data,cc)

##################################
#### 7.ssGSEA ####
##################################
for (i in names(NE.expr.list)) {
  pca.data = NE.expr.list[[i]]
  phe.data = NE.phe.list[[i]]
  
  sig_tme <- calculate_sig_score(pdata           = NULL,
                                 eset            = pca.data,
                                 signature       = ne.markers.list,
                                 method          = "ssgsea",
                                 mini_gene_count = 5)
  sig_tme = arrange(sig_tme,Index)
  phe.data$NE_UP_DN = sig_tme$UP-sig_tme$DN
  phe.data$NE_UP = sig_tme$UP
  NE.phe.list[[i]] = phe.data
}

## NE_UP
cc <- data.frame(
  R2=sapply(NE.phe.list,function(x){as.numeric(R2(x$NE_UP%>% normalize(), x$NE_status)) }),
  RMSE=sapply(NE.phe.list,function(x){as.numeric(RMSE(x$NE_UP%>% normalize(), x$NE_status)) }),
  MSE=sapply(NE.phe.list,function(x){as.numeric(mean((x$NE_UP%>% normalize()-x$NE_status)^2)) }),
  MAE=sapply(NE.phe.list,function(x){as.numeric(MAE(x$NE_UP%>% normalize(), x$NE_status)) }),
  Cindex=sapply(NE.phe.list,function(x){as.numeric(CstatisticCI(rcorr.cens(x$NE_UP %>% normalize(), x$NE_status))[1])}),
  AUC = sapply(NE.phe.list,function(x){round(roc(NE_status~NE_UP%>% normalize(),data=x,aur=TRUE,ci=TRUE)$auc ,2)})
)%>%
  rownames_to_column('ID')

result.data <- data.frame()
cc$Model <- 'UP_ssGSEA'
result.data <- rbind(result.data, cc)

## NE_UP_DN
cc <- data.frame(
  R2=sapply(NE.phe.list,function(x){as.numeric(R2(x$NE_UP_DN, x$NE_status))}),
  RMSE=sapply(NE.phe.list,function(x){as.numeric(RMSE(x$NE_UP_DN, x$NE_status))}),
  MSE=sapply(NE.phe.list,function(x){as.numeric(mean((x$NE_UP_DN-x$NE_status)^2))}),
  MAE=sapply(NE.phe.list,function(x){as.numeric(MAE(x$NE_UP_DN, x$NE_status))}),
  Cindex=sapply(NE.phe.list,function(x){as.numeric(CstatisticCI(rcorr.cens(x$NE_UP_DN, x$NE_status))[1])}),
  AUC = sapply(NE.phe.list,function(x){round(roc(NE_status~NE_UP_DN,data=x,aur=TRUE,ci=TRUE)$auc,2)})
)%>%
  rownames_to_column('ID')

cc$Model <- 'UP_DN_ssGSEA'
result.data <- rbind(result.data, cc)
result.data 

### Vis
p1 = result.data %>%
  ggplot(aes(AUC,reorder(Model,AUC)))+
  geom_bar(width = 0.7,stat = 'summary',fun='mean',fill='#0072B5FF')+
  stat_summary(fun = mean,
               geom = "pointrange",size=0.08,col="black",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x))+
  theme_classic()+
  labs(y=NULL);p1

p2 = result.data %>%
  ggplot(aes(Cindex,reorder(Model,Cindex)))+
  geom_bar(width = 0.7,stat = 'summary',fun='mean',fill='#0072B5FF')+
  stat_summary(fun = mean,
               geom = "pointrange",size=0.08,col="black",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x))+
  theme_classic()+
  labs(y=NULL);p2

p3 = result.data %>%
  ggplot(aes(R2,reorder(Model,R2)))+
  geom_boxplot(width = 0.7,fill='#E18727FF')+
  stat_summary(fun = mean,
               geom = "pointrange",size=0.08,col="black",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x))+
  theme_classic()+
  labs(y=NULL);p3

## heatmap
input.data = result.data 
input.data$ID = factor(input.data$ID,levels = c("SU2C","PCaProfilter",
                                                "WCDT","UM_SPORE","WCM","scRNAseq"))
if(F){
  dd2 <- input.data[order(input.data$Cindex,decreasing = T),]
  dt <- dplyr::select(dd2,ID,Cindex,Model)
  dt = pivot_wider(dt, names_from = ID,
                   values_from = Cindex) %>% as.data.frame()
  dt = dt %>% column_to_rownames("Model")  
  dt = dt[order.ID,]
  col_ha <- HeatmapAnnotation(which = "col", Cohort = colnames(dt),
                              annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
                              annotation_name_side = "left",
                              col = list(Cohort=c("SU2C"="#bc3b29",
                                                  "PCaProfilter"="#0672b7",
                                                  "WCDT"="#e28729",
                                                  "UM_SPORE"="#208650",
                                                  "scRNAseq"="#7776b0",
                                                  "WCM" = "#6F99ADFF")),
                              annotation_legend_param = list(Cohort=list(title = "Cohort",
                                                                         title_position = "topleft",
                                                                         title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                                         labels_rot = 0,
                                                                         legend_height = unit(1,"cm"),
                                                                         legend_width = unit(5,"mm"),
                                                                         labels_gp = gpar(fontsize = 9,
                                                                                          fontface = "bold"))
                              )
  )
  
  row_ha <- rowAnnotation('Mean Cindex' = anno_barplot(round(rowMeans(dt), 3), bar_width = 1, add_numbers = T,
                                                       labels = c("Mean Cindex"), height = unit(1, "mm"),
                                                       gp = gpar(col = "white", fill = "skyblue1"), numbers_gp = gpar(fontsize = 8),
                                                       axis_param = list(at = c(0, 0.5, 1),
                                                                         labels = c("0", "0.5", "1")),
                                                       width = unit(2.5, "cm")),
                          annotation_name_side = "bottom",
                          annotation_name_gp = gpar(fontsize = 9, fontface = "bold", angle = 90))
  
 cell_fun <- function(j, i, x, y, width, height, fill) {
    grid.text(
      round(dt[i, j], 2), 
      x, y,
      gp = gpar(
        fontsize = 8
      ))
  }
  
  heatmap <- Heatmap(dt,name = " ",
                     heatmap_legend_param = list(title="",title_position = "topleft", labels_rot = 0,
                                                 legend_height = unit(8,"cm"),
                                                 legend_width = unit(5,"mm"),
                                                 labels_gp = gpar(fontsize = 15, fontface = "bold")),
                     border = TRUE,
                     column_split = colnames(dt),
                     column_gap = unit(3, "mm"),
                     show_column_names = F,
                     show_row_names = T,
                     col = colorRamp2(c(0.7,0.85,1), 
                                      c("#4DBBD5B2", "white", "#E64B35B2")),
                     column_title ="", 
                     column_title_side = "top", 
                     row_title_side = "left",
                     row_title_rot = 90, 
                     column_title_gp = gpar(fontsize = 12, fontface = "bold",col = "black"), 
                     cluster_columns =F,
                     cluster_rows = F,
                     column_order = c(colnames(dt)),
                     show_row_dend =F, 
                     cell_fun = cell_fun,
                     top_annotation = col_ha
  )
  draw(heatmap)
}

## ROC plot
input.data = res.sum$PCaProfilter
res<-roc(NE_status~NE_UP_DN_ssGSEA+NE_UP_ssGSEA+alpha.0+alpha.0.9,data=input.data,aur=TRUE,
         ci=TRUE,
         smooth=T
)
p = ggroc(res,legacy.axes = TRUE)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw() + 
  scale_color_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#7876B1FF"))+
  ggtitle("PCaProfilter (NEPC/total = 19/1223)")+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text=element_text(size=10,colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))
p

p+annotate("text",x=0.85,y=0.35,size = 3,label=paste("UP-DN AUC: ", round(res$NE_UP_DN_ssGSEA$auc,3)))+
  annotate("text",x=0.85,y=0.25,size = 3,label=paste("Enet (α=0.01) AUC: ", round(res$alpha.0$auc,3)))+
  annotate("text",x=0.85,y=0.15,size = 3,label=paste("Enet (α=0.9) AUC: ", round(res$alpha.0.9$auc,3)))+
  annotate("text",x=0.85,y=0.05,size = 3,label=paste("UP AUC: ", round(res$NE_UP_ssGSEA$auc,3)))
