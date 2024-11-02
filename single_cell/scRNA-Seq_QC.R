# Library and Load data -----------------
library(Matrix)
library(Seurat)
library(tidyverse)
setwd("/home/chenweiming/Project/HCC_scRNAseq/luo/")

seu = readr::read_rds("/home/bioinformatics/Project/zhluo/HCC_wangfubing/seu_A124.counts_anno.rds")
seu

# QC质控 ----------------------------- 
dir.create('QC')
proj_name <- data.frame(proj_name=rep("QC",ncol(seu)))
rownames(proj_name) <- row.names(seu@meta.data)
seu <- AddMetaData(seu, proj_name)
rm(proj_name)
head(seu@meta.data)
Vln<-VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3, group.by = "proj_name",raster=FALSE)
Vln
ggsave("QC/QC-before.pdf", plot = Vln, width = 7, height = 5.5)

dim(seu)

##对feature进行质控，设置质控标准
minicount = 500 #每个基因在细胞中的表达量
maxicount = 30000 #每个基因在细胞中的表达量
minGene=500  #每个细胞中的最小基因数
maxGene=6000 #每个细胞中的最大基因数
mt=25        #线粒体基因的UMI总数超过mt%的细胞

#filtering 
seu <- subset(seu, subset = percent.mt < mt & nFeature_RNA >= minGene & nFeature_RNA <= maxGene & nCount_RNA >= minicount & nCount_RNA <= maxicount)

##删掉线粒体和核糖体的基因
genes_obtained = names(rowSums(seu@assays$RNA@counts)[rowSums(seu@assays$RNA@counts) >= 3])
mito.genes <- grep('^MT-', genes_obtained, value = TRUE)
ribo.genes <- grep('^RPL|^RPS|^MRPL|^MRPS', genes_obtained, value = TRUE)
genes_obtained = genes_obtained[!genes_obtained %in% c(mito.genes, ribo.genes)]
seu <- seu[rownames(seu) %in% genes_obtained, ]

table(seu@meta.data$Cancer_type)
table(seu@meta.data$Sample)
seu@meta.data$Label<-seu@meta.data$Sample
# 删除 `seu2@meta.data$Label` 列中第一个 `_` 符号及之前的字符
seu@meta.data$Label <- sub("^[^_]*_", "", seu@meta.data$Label)
# 确保 `Label` 列被正确赋值
table(seu@meta.data$Label)
table(seu@meta.data$Cancer_type)

# 只保留HCC
seu = subset(seu, Cancer_type == "HCC")
table(seu@meta.data$Cancer_type)
table(seu@meta.data$Label)

## 载入临床信息----
clin<-read.csv('/home/chenweiming/Project/HCC_scRNAseq/data/Sample_information.CSV')
# 替换列名中的两个连续点为一个下划线，并移除末尾的点
colnames(clin) <- gsub("\\.", "_", colnames(clin))  # 将所有点替换为下划线
colnames(clin) <- gsub("___", "_", colnames(clin))  # 将两个连续下划线替换为一个下划线
colnames(clin) <- gsub("__", "_", colnames(clin))  # 将两个连续下划线替换为一个下划线
colnames(clin) <- gsub("_$", "", colnames(clin))   # 移除末尾的下划线
clin$Sample<-paste(clin$Patient,clin$Cancer_type_short,sep = '_')
unique(clin$Sample)
head(clin)

# 查看CEA缺失情况
cea_value_na = clin$Sample[is.na(clin$CEA_μg_L)]
cea_value_na
# cea_value_na = c("A080_HCC", "A088_HCC", "A120_HCC")
# cea_stage_na = clin$Sample[is.na(clin$CEA_status)]
# cea_stage_na
table(seu$Sample)[cea_value_na]

# 合并几个样本
seu$Sample[seu$Sample %in% c('A074_HCC','A074_HCC_IM')]<-'A074_HCC'
seu$Sample[seu$Sample %in% c('A119_HCC','A119_HCC_IM1','A119_HCC_IM2')]<-'A119_HCC'

# 删掉CEA含量缺失的样本："A019_HH"  "A080_HCC" "A088_HCC" "A120_HCC"
`%!in%` <- Negate(`%in%`)
samples = unique(seu@meta.data$Sample)
samples
samples = samples[samples %!in% cea_value_na]
samples
seu@meta.data$selected_T_F = seu@meta.data$Sample
seu@meta.data$selected_T_F[seu@meta.data$selected_T_F %in% samples] <- "sel_samples"
table(seu@meta.data$selected_T_F)

seu <- subset(seu, selected_T_F == "sel_samples")
seu
library(dplyr)
seu@meta.data <- select(seu@meta.data, -selected_T_F)
head(seu@meta.data)

Vln<-VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, group.by = "proj_name",ncol = 3,raster=FALSE)
Vln
ggsave("QC/QC-after.pdf", plot = Vln, width = 7, height = 5.5)


# 合并临床信息
seu@meta.data <- seu@meta.data %>%
  left_join(clin, by = "Sample")
table(seu$CEA_status)
sum(is.na(seu$CEA_status))
table(seu$Sample[is.na(seu$CEA_μg_L)])

##恢复meta.data行名Barcode
rownames(seu@meta.data)<-colnames(seu)

## 根据中位数分组----
median_value <- median(seu@meta.data$CEA_μg_L, na.rm = TRUE)
seu@meta.data$CEA_Group <- ifelse(seu@meta.data$CEA_μg_L > median_value, 'High_CEA', 'Low_CEA')
head(seu@meta.data)


# 细胞注释：原文章的细胞类型大类 ------------------------------------------------------
seu@meta.data$Prim_clusters = seu@meta.data$clusters
table(seu@meta.data$Prim_clusters)
# Create a function to map clusters to categories based on the prefix
map_clusters_to_category <- function(cluster) {
  if (cluster %in% c("CD4T_09_FOXP3", "CD4T_10_FOXP3_CTLA4", "CD4T_11_FOXP3_STMN1")) {
    return("T reg cells")
  } else if (startsWith(cluster, "CD4T")) {
    return("CD4+ T cells")
  } else if (startsWith(cluster, "CD8T")) {
    return("CD8+ T cells")
  } else if (startsWith(cluster, "B_")) {
    return("B cells")
  } else if (startsWith(cluster, "NK")) {
    return("NK cells")
  } else if (startsWith(cluster, "DC")) {
    return("Dendritic cells")
  } else if (cluster == "MonoDC") {
    return("Dendritic cells")
  } else if (startsWith(cluster, "Mo_")) {
    return("Monocytes")
  } else if (startsWith(cluster, "Mono-like")) {
    return("Monocytes")
  } else if (startsWith(cluster, "Mph")) {
    return("Macrophages")
  } else if (startsWith(cluster, "gdT")) {
    return("γδT cells")
  } else if (startsWith(cluster, "EC")) {
    return("Endothelial cells")
  } else if (startsWith(cluster, "Fb")) {
    return("Fibroblasts")
  } else if (startsWith(cluster, "Mu")) {
    return("Muscle Cells")
  } else if (startsWith(cluster, "Neu")) {
    return("Neutrophils")
  } else if (cluster == "Mast") {
    return("Mast cells")
  } else if (cluster == "Tumor") {
    return("Tumor")
  } else {
    return("Unknown")
  }
}

# Apply the function to create a new column in the metadata
seu@meta.data$Prim_Celltype1 <- sapply(seu@meta.data$Prim_clusters, map_clusters_to_category)
table(seu@meta.data$Prim_Celltype1)
table(seu@meta.data$Cancer_type)

# 删掉A068这个样本，因为这个样本跟其他样本聚类的时候基本没有overlap的区域
seu = seu[, seu@meta.data$Sample != "A068_HCC"] #19057*88202
# 样本A079_HCC，A080_HCC，A102_HCC，A120_HCC分别只有5,6,4,7个细胞，这些细胞的作用不大。
seu = seu[, !seu@meta.data$Sample %in% c("A079_HCC", "A080_HCC", "A102_HCC", "A120_HCC")]

saveRDS(seu, file = "/home/chenweiming/Project/HCC_scRNAseq/luo/data/HCC_scRNA.rds")