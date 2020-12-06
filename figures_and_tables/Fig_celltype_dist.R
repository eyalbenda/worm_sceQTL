source("analysis_functions//1_load_data_and_functions.R")

cellPlot = plot_cells(cdsFinBC,color_cells_by = "broad_tissue",label_cell_groups = F)
umap_colorPallete  = ggplot_build(cellPlot)$plot$scales$scales[[1]]$palette.cache
names(umap_colorPallete) = ggplot_build(cellPlot)$plot$scales$scales[[1]]$range$range


ordColor = c(1:5,14,15,8:12,16,18,7,13,6,19,17)

data.frame(nUMI=colSums(exprs(cdsFinBC))) %>% 
  group_by(`Cell Type`=colData(cdsFinBC)$broad_tissue) %>% 
  ggplot(aes(x=`Cell Type`,y=nUMI,color=`Cell Type`,fill = `Cell Type`)) + 
    geom_violin() + 
    scale_y_log10() + 
    theme_cowplot() +   
    scale_color_manual(breaks=names(umap_colorPallete)[ordColor],values=umap_colorPallete[ordColor]) +
    scale_fill_manual(breaks=names(umap_colorPallete)[ordColor],values=umap_colorPallete[ordColor]) + 
    theme(legend.position = "none",axis.text.x = element_text(angle=45,hjust=1))
ggsave(file.path(out_basepath,"Figures/Supplementary Figures/nUMI_per_celltype.png"),dpi=300)
