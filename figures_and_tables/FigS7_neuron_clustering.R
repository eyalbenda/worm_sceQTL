source("analysis_functions/1_load_data_and_functions.R")

plot_cells(neuCDS,color_cells_by = "fine_tissue",label_groups_by_cluster = F)
ggsave(file.path(out_basepath,"Figures/Figure 3/neuronalClustering.png"),dpi=600,width=4,height=4)
