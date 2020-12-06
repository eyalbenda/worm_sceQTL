source("analysis_functions/1_load_data_and_functions.R")

plotQTLdf(do.call(rbind,lapply(1:length(tissuePeaks),function(x)generatePlotTissue(tissuePeaks[[x]],tissueCis[[x]],dfReturn = T))))

#ggsave(file.path(out_basepath,"Figures/Figure 2/global_map_by_tissuemapping.png"),dpi=600,height=4,width=4)
