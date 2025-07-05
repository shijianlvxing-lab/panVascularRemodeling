

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


