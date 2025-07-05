
#细胞周期计算+作图函数
qc_cc <- function(input_sce){
  input_sce = CellCycleScoring(object = input_sce, 
                               slot = "counts",
                               g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
  p1=VlnPlot(input_sce, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", ncol = 2, pt.size = 0.1)
  p1
  ggsave('vlnplot_cell_cycle.pdf')
  p2= ggplot(input_sce@meta.data,aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+theme_minimal()
  p2
  ggsave('dotplot_cell_cycle.pdf')
  return(input_sce)
}


#基因数、线粒体、核糖体、血红蛋白质控，细胞周期计算
basic_qc <- function(input_sce){
  library(Seurat)
  input_sce <- AddMetaData(input_sce,metadata = PercentageFeatureSet(input_sce, pattern = "^MT-"),col.name = 'percent.mt')
  input_sce <- AddMetaData(input_sce,metadata = PercentageFeatureSet(input_sce, pattern = "^RB[SL]"),col.name = 'percent.rb')
  input_sce <- AddMetaData(input_sce,metadata = PercentageFeatureSet(input_sce, pattern = "^HB[^(p)]"),col.name = 'percent.hb')
  
  feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb", "percent.hb")
  
  p1=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0,ncol = 1) + 
    NoLegend()
  p1 
  ggsave(filename="vlnplot_before.pdf",plot=p1,width = 6,height = 20)
  
  
  input_sce <- subset(input_sce,subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25 & percent.hb <1 & percent.rb<1)
  
  
  p1=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0,ncol = 1) + 
    NoLegend()
  p1 
  ggsave(filename="vlnplot_after.pdf",plot=p1,width = 6,height = 20)
  
  input_sce <- qc_cc(input_sce)
    return (input_sce)

}
