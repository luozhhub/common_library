library(Matrix)
library(Seurat)
library(tidyverse)

seu = readr::read_rds("/home/bioinformatics/Project/zhluo/HCC_wangfubing/seu_A124.counts_anno.rds")

table(seu@meta.data$Cancer_type)
table(seu@meta.data$Sample)

seu@meta.data$Label<-seu@meta.data$Sample
# 删除 `seu2@meta.data$Label` 列中第一个 `_` 符号及之前的字符
seu@meta.data$Label <- sub("^[^_]*_", "", seu@meta.data$Label)

# 确保 `Label` 列被正确赋值
table(seu@meta.data$Label)

table(seu@meta.data$Cancer_type)


##提取恶性细胞----
Hep = seu[, HCC@meta.data$Label %in% c('HCC','HCC_IM', 'HCC_IM1', 'HCC_IM2') & HCC@meta.data$clusters %in% 'Tumor']

##创建Seurat对象
Hep<-CreateSeuratObject(counts = Hep@assays$RNA@counts,
                        meta.data = Hep@meta.data)


#筛选高可变基因，并进行数据降维
Hep <- NormalizeData(Hep) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)


##去批次----
library(harmony)
unique(Hep$orig.ident)
##使用harmony进行整合&批次矫正 
Hep <- RunHarmony(Hep,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
#非线性降维，UMAP和tSNE 是Seurat进行非线性降维的方法，首选 UMAP
Hep <- RunUMAP(Hep, reduction = "harmony", dims = 1:30,reduction.name = "umap")

#查看去批次效果
DimPlot(Hep, reduction = "umap",group.by = "orig.ident")

##识别细胞簇(cluster)
Hep <- FindNeighbors(Hep, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.4)  

DimPlot(Hep,reduction = "umap",label=T)+ggtitle('HCC Tumor cell')

save(Hep,file = '~/DATA/luo/tongji/HCC_nature/HepTumor.rdata')

#load( '/home/zhaojingwei//DATA/luo/tongji/HCC_nature/HepTumor.rdata')

##载入临床信息----
clin<-read.csv('/home/zhaojingwei/DATA/luo/tongji/HCC_nature/Nature_HCC_Clin.csv')
clin$Sample<-paste(clin$Patient,clin$Cancer_type_short,sep = '_')
unique(clin$Sample)

unique(Hep$Sample)

##统一Sample信息----
Hep$Sample[Hep$Sample %in% c('A074_HCC','A074_HCC_IM')]<-'A074_HCC'
Hep$Sample[Hep$Sample %in% c('A119_HCC','A119_HCC_IM1','A119_HCC_IM2')]<-'A119_HCC'

Hep@meta.data <- Hep@meta.data %>%
  left_join(clin, by = "Sample")

##恢复meta.data CB----
rownames(Hep@meta.data)<-colnames(Hep)

DimPlot(Hep)

DimPlot(Hep,group.by = 'TNM_stage')
