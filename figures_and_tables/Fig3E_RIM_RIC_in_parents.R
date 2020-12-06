source("analysis_functions/1_load_data_and_functions.R")

neuCDSparents = readRDS("study_data/neuCDS_parents.rds")

strainRIC = colData(neuCDSparents)$Strain[clusters(neuCDSparents)==26]
strainRIM = colData(neuCDSparents)$Strain[clusters(neuCDSparents)==21]
exprRIM = normalized_counts(neuCDSparents[geneAnnotsCur$wbps_gene_id[geneAnnotsCur$external_gene_id=="nlp-21"],clusters(neuCDSparents) %in% c(21)])
exprRIC = normalized_counts(neuCDSparents[geneAnnotsCur$wbps_gene_id[geneAnnotsCur$external_gene_id=="nlp-21"],clusters(neuCDSparents)==26])



exprAll = normalized_counts(neuCDSparents[geneAnnotsCur$wbps_gene_id[geneAnnotsCur$external_gene_id=="nlp-21"],])
strainAll = colData(neuCDSparents)$Strain


df = rbind(data.frame(neuron="RIM",expression=as.vector(exprRIM),Strain=strainRIM),data.frame(neuron="RIC",expression=as.vector(exprRIC),Strain=strainRIC),
           data.frame(neuron="All",expression=as.vector(exprAll),Strain=strainAll))

ggplot(na.omit(df),aes(x=Strain,y=expression,color=Strain,fill=Strain,group=Strain)) + 
  geom_crossbar(stat="summary", fun="mean") +
  geom_point(size=1, position = position_jitter(w=0.3)) + theme_classic() + facet_grid(~neuron) + scale_color_viridis_d(end=0.6) + scale_fill_viridis_d(end=0.6) +
  theme(legend.position = c(0.12,0.9),legend.background = element_rect(fill="transparent")) + ylab("Normalized Expression") + guides(color = guide_legend(NULL),fill = guide_legend(NULL))
ggsave(file.path(out_basepath,"Figures/Figure 3/nlp21_parents.png"),width=4,height=4,dpi=300)
