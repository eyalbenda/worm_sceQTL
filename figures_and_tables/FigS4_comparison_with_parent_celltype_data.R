source("analysis_functions//3_comparison_to_parental_dataset.R")

corDF$x = min(corDF$x)

ggplot(dfCisAndParents %>% filter(FDR<0.1),aes(x=log2(norm_foldChange),y=negbin.beta/negbin.se)) + 
  geom_point(size=0.2,aes(color=FDR<0.1)) + 
  theme_classic()  + 
  facet_wrap(~tissue) + 
  scale_colour_manual(values=c(viridis::inferno(10)[2])) + 
  geom_text(data=corDF,aes(x=x,y=y,label=label),parse=T,hjust=0) +
  ylab(expression("standardized "*beta*" from cis-eQTL model")) + 
  xlab(expression(log[2]*"(fold change)")) + 
  theme(legend.position = "none")
ggsave(file.path(out_basepath,"Figures/Supplementary Figures/Previous Data/TissueSpecificAndParent.png"),width=7,height=8.5,dpi=300)

