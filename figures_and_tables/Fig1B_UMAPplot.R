source("analysis_functions//1_load_data_and_functions.R")

cellPlot = plot_cells(cdsFinBC,color_cells_by = "broad_tissue",label_cell_groups = F)
umap_colorPallete  = ggplot_build(cellPlot)$plot$scales$scales[[1]]$palette.cache
names(umap_colorPallete) = ggplot_build(cellPlot)$plot$scales$scales[[1]]$range$range


ordColor = c(1:5,14,15,8:12,16,18,7,13,6,19,17)
cellPlot + theme(legend.position="bottom",legend.spacing = unit(0,"cm"),legend.margin=margin(0,0,0,0),
                 legend.key.size = unit(0, 'lines'),
                 legend.box.margin=margin(0,0,0,0)) + ylab("") + xlab("") + guides(color=guide_legend(ncol=3,title = "Cell type",title.position = "top",keywidth = 0.4,keyheight = 0.4,override.aes = list(size = 5))) + 
  scale_color_manual(breaks=names(umap_colorPallete)[ordColor],values=umap_colorPallete[ordColor])

#ggsave(file.path(out_basepath,"Figures/Figures/Figure 1/UMAPplot.png"),dpi=400,width=4.5,height=6)
# Note: letters in legend needed to be entered manually