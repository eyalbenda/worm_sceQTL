source("analysis_functions/6_single_neuron_eQTL_processing.R")
assertthat::are_equal(colnames(neuCDS),rownames(genoMatrixScEQTL))


plot_sneQTL(genoMatrixScEQTL,neuCDS,neurDF,which(neurDF$gene_short_name=="nlp-21")[1],normalize=F)
ggsave(file.path(out_basepath,"Figures/Figure 3/nlp21_RIC.png"),width=4,height=4,dpi=300)

plot_sneQTL(genoMatrixScEQTL,neuCDS,neurDF,which(neurDF$gene_short_name=="nlp-21")[2],normalize=F)
ggsave(file.path(out_basepath,"Figures/Figure 3/nlp21_RIM.png"),width=4,height=4,dpi=300)
