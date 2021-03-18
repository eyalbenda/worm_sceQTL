source("analysis_functions/6_single_neuron_eQTL_processing.R")

looDF = neurDF %>% filter(panNeuFDR<0.5) %>% summarize(same = sum(sign(looBetaNorm)==sign(negbin.beta)),
                                                       opposite=sum(sign(looBetaNorm)!=sign(negbin.beta)),thresh="<0.5")
looDF = rbind(looDF,neurDF %>% filter(panNeuFDR>0.5) %>% summarize(same = sum(sign(looBetaNorm)==sign(negbin.beta)),
                                                                   opposite=sum(sign(looBetaNorm)!=sign(negbin.beta)),thresh=">0.5"))
looDF = reshape2::melt(looDF)



ggplot(looDF,aes(x=thresh,color=variable,fill=variable,y=value)) + 
  geom_bar(position="dodge",width=0.2,stat = "identity") + theme_classic() + 
  scale_color_viridis_d(end=0.6) + scale_fill_viridis_d(end=0.6) + 
  ggplot2::guides(colour=guide_legend(expression("effect\ndirection")),fill=guide_legend(expression("effect\ndirection"))) + xlab("pan neuronal FDR") + ylab("Number of eQTL") + theme(legend.position=c(0.879,0.85))

#ggsave(file.path(out_basepath,"Figures/Fig3E.png"),width=4,height=4,dpi=300)
#Note: this plot ended up redone by Stefan Zdraljevic manually as a diagram based on the output here
