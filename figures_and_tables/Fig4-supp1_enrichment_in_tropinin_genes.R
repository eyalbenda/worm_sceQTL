source("analysis_functions//1_load_data_and_functions.R")
plot_genes_by_group(cdsFinBC,markers=c("tni-3","pat-10","mup-2",
                                       "tni-1","tnt-2","tnt-3","tnc-2",
                                       "tni-4","tnt-4"),group_cells_by = "broad_tissue") + xlab("Cell type") + coord_flip()
ggsave(file.path(out_basepath,"Figures/Supplementary Figures//Figure S6.png"),dpi=300,width=7,height=7)
