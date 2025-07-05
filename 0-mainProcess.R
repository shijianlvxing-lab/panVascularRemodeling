#Step0:read data----
#读取10×标准文件
if(T){
  rm(list = ls());gc();
  setwd("D:/科研工作站/生信工作站/3-dataSet/GSE224559")
  library(Seurat)
  dir <- "data/"#设置可读数据目录
  samples <- list.files(dir)
  samples
  sceList <-  lapply(samples,function(samples){
    sce <-  CreateSeuratObject(counts =  Read10X(file.path(dir,samples)),project = samples,min.cells = 3,min.features = 200 ) 
    return(sce)})#逐个创建Seurat，储存为列表
  names(sceList) =  samples
  sce_merged <- merge(sceList[[1]], y= sceList[-1] ,add.cell.ids =  samples)#合并列表中的Seurat对象
  sce_merged <- JoinLayers(sce_merged)#合并层???
  dim(sce_merged[["RNA"]]$counts )#打印细胞数、基因数
  
  
  #给细胞添加样本来源信息
  sampleInfo <- readxl::read_xlsx('meta/sampleInfo.xlsx')#读取样本信息
  sampleName <- sce_merged$orig.ident
  column <- colnames(sampleInfo)#预添加信息列（默认前10项）
  order <- match(sampleName, sampleInfo$sampleID)#每一个样本在样本信息中的次序
  
  for(col in column){
    sce_merged <- AddMetaData(object = sce_merged,metadata = sampleInfo[[col]][order],col.name = col)
  }
  
  saveRDS(sce_merged,"sce_merged.rds")
}



#Step1:QC-----------
#basic QC----
if(T){
  rm(list = ls());gc();
  sce_merged <- readRDS("sce_merged.rds")
  source('D:/科研工作站/生信工作站/0-code/v2/1-qc.R')
  dir.create("./1-qc")
  setwd('1-qc/')
  sce_qc = basic_qc(sce_merged)
  
  setwd('../')
  saveRDS(sce_qc,"sce_qc.rds")
}

#doubleFinder----
if(T){
  rm(list = ls());gc();
  sce_qc <- readRDS("sce_qc.rds")
  source('D:/科研工作站/生信工作站/0-code/v2/3-doubleFinder.R')
  dir.create("./2-df")
  setwd('2-df/')
  sce_df = df_qc(sce_qc)
  setwd('../')
  saveRDS(sce_qc,"sce_df.rds")
} 
