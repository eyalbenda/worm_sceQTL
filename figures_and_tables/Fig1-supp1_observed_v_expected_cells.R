source("analysis_functions//1_load_data_and_functions.R")

cellsFrequencies = rep(NA,length(unique(colData(cdsFinBC)$broad_tissue)))
names(cellsFrequencies) = unique(colData(cdsFinBC)$broad_tissue)
cellsFrequencies["XXX"] = 2
### Neurons are 222 at birth + 21 post-embryonic from wormatlas + 1 CEM (male)
cellsFrequencies["Neuron"] = 244
cellsFrequencies["Intestine"] = 20
cellsFrequencies["Body Wall Muscle"] = 94
cellsFrequencies["Somatic Gonad"] = 10
cellsFrequencies["Pharynx and Arcade Cells"] = 51
#### Total number of hypodermis nuclei minus the 22 H,V,T nuclei that are made post L2
cellsFrequencies["Hypodermis"] = 143
#### Need to remove a bunch of L3-L4 cells later
cellsFrequencies["Seam Cells"] = 60
cellsFrequencies["Coelomocytes"] = 6
#### Counted based on wormatlas drawing
cellsFrequencies["Germline"] = 17
cellsFrequencies["Excretory Cells"] = 3
cellsFrequencies["Excretory Gland"] = 2
cellsFrequencies["Glia"] = 51
cellsFrequencies["GLR"] = 6 
cellsFrequencies["Pharyngeal Gland Cells"] = 4 
cellsFrequencies["Sex Myoblast"] = 14
cellsFrequencies["Sphincter and Anal Muscles"] = 4
cellsFrequencies["Vulval Precursor Cells"] = 6 

cellPlot = plot_cells(cdsFinBC,color_cells_by = "broad_tissue",label_cell_groups = F)
umap_colorPallete  = ggplot_build(cellPlot)$plot$scales$scales[[1]]$palette.cache
names(umap_colorPallete) = ggplot_build(cellPlot)$plot$scales$scales[[1]]$range$range


ordColor = c(1:5,14,15,8:12,16,18,7,13,6,19,17)



cellPlotColors = tapply(as.character(ggplot_build(cellPlot)$data[[1]]$colour),colData(cdsFinBC)$broad_tissue,unique)

expected = cellsFrequencies[sort(names(cellsFrequencies))]*55508/  
  sum(cellsFrequencies[sort(names(cellsFrequencies))],na.rm=T) ## na.rm here required due to the Unknown cell cluster having no expected value
obsExpDF = data.frame(observed = table(colData(cdsFinBC)$broad_tissue),
                      expected = expected)
ggplot(obsExpDF,aes(x=expected,y=observed.Freq,color=observed.Var1,label=observed.Var1))+
  geom_point() + 
  theme_classic() + 
  geom_text_repel()  + 
  theme(legend.position = "none") + 
  ylab("Observed number of cells") + 
  xlab("Expected number of cells") + 
  scale_x_log10(limits=c(100,20000)) + 
  scale_y_log10(limits=c(100,20000))
#ggsave(file.path(out_basepath,"Figures/Supplementary Figures/expectedObservedCells.png"),dpi=300,width=6,height=6)


