
#单样本doubletFinder分析
run_doublet_finder <- function(sce, doublet_rate=ncol(sce)*8*1e-6, pcs = 1:15) {
  library(DoubletFinder)
  #step0:数据预处理
  sce <- NormalizeData(sce)
  sce <- FindVariableFeatures(sce)
  sce<- ScaleData(sce)
  sce<- RunPCA(sce)
  sce <- FindNeighbors(sce)
  sce <- FindClusters(sce)
  
  #step1:确定最近邻比例pK值
  sweep.res <- paramSweep(sce, PCs = pcs)#参数扫描，不同pN-pK组合下每个细胞的pANN
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.vector() %>% as.numeric()#确定最优pK值，使得所有pN下BCmetric最大的pK值
  #step2:估算异源双细胞的数量
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)#估计同源双细胞的比例
  nExp_poi <- round(doublet_rate * ncol(sce))#估计双细胞总数
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))#估算异源双细胞数量
  #step3：根据每个细胞的最近邻比例和异源双细胞总数和来判断每个细胞是否为异源双细胞
  sce <- doubletFinder(sce, 
                       PCs = pcs, 
                       pN = 0.25, #默认值
                       pK = pK_bcmvn,
                       nExp = nExp_poi.adj, 
                       reuse.pANN = FALSE)
  return(sce)
  
}

#合并的Seurat处理
df_qc <- function(input_sce){
  #step0:分割样本
  library(Seurat)
  library(BiocParallel)
  register(MulticoreParam(workers = 12, progressbar = TRUE))
  
  sce_list <- SplitObject(input_sce, split.by = "orig.ident")
  
  #step1:每个样本单独运行run_doublet_finder  
  sce_list <- lapply(sce_list,run_doublet_finder)
    
  
  #step2:去掉pANN 和DF.classifications列名的后缀，以合并不同的样本
  sce_list <- lapply(sce_list, function(sce) {
      # pANN
      pANN_col <- grep("^pANN_0.25", colnames(sce@meta.data), value = TRUE)
      if (length(pANN_col) > 0) {
        colnames(sce@meta.data)[colnames(sce@meta.data) %in% pANN_col] <- "pANN_0.25"
      }
      
      # DF.classifications
      DF_col <- grep("^DF\\.classifications_0.25", colnames(sce@meta.data), value = TRUE)
      if (length(DF_col) > 0) {
        colnames(sce@meta.data)[colnames(sce@meta.data) %in% DF_col] <- "DF.classifications"
      }
      
      return(sce)
    })
  
    #
    input_sce <- merge(sce_list[[1]], y= sce_list[ -1 ])
    input_sce <- JoinLayers(input_sce)
    write.table(table(input_sce$DF.classifications),'doublet_stat.txt')
    input_sce.filtered <- subset(input_sce, subset = DF.classifications == "Singlet")
    return(input_sce.filtered)
}



