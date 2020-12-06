source("analysis_functions/1_load_data_and_functions.R")

alleQTL = list()
for(tis in names(tissuePeaks)) alleQTL[[tis]] = generatePlotTissue(tissuePeaks[[tis]],tissueCis[[tis]],dfReturn=T)
dfeQTL = reshape2::melt(do.call(rbind,lapply(alleQTL,function(x)table(factor(x$type,levels=c("cis","trans")))))) %>% dplyr::rename(Tissue=Var1,type=Var2,eqtl=value)
fct_order = levels(fct_reorder(dfeQTL$Tissue, dfeQTL$eqtl))
ggplot(dfeQTL ,aes(x=factor(Tissue, levels=fct_order),y=eqtl,color=type,fill=type)) + 
  geom_bar(stat="identity",width=0.4,position = position_dodge()) + 
  theme_cowplot() + scale_color_viridis_d(end=0.6,labels = c(expression(italic("cis")), expression(italic("trans")))) + 
  scale_fill_viridis_d(end=0.6,labels = c(expression(italic("cis")), expression(italic("trans")))) + coord_flip() + 
  theme(axis.text.x = element_text(angle=90,vjust=0.5),legend.position=c(0.63,0.12)) + 
  xlab("Cell type") + ylab(expression("eQTL"))

#ggsave("~/hujidrive/Work/Research/Postdoc/scRNAseq/Manuscript/Figures/Figure 2/eQTL_numbers.png",dpi=600,width=4.32,height=3.96)
