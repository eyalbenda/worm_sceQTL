source("analysis_functions//4_comparison_to_wholeworm_data.R")

corDF = with(dfCisAndRockmanCis %>% filter(significant),makeCorText(previous,ourdata))
ggplot(dfCisAndRockmanCis %>% filter(significant),aes(x=previous,y=ourdata)) + geom_point(size=0.2,aes(color=significant)) + theme_classic() +
  scale_colour_manual(values=c(viridis::inferno(10)[2]))  + theme(legend.position = "none") + geom_text(parse=T,data=corDF,aes(x=x,y=y,label=label),hjust=0,vjust=0) +
  ylab(expression("standardized "*beta*" from cis-eQTL model")) + xlab(expression("effect size in the bulk study")) + theme(legend.position = "none")
#ggsave(file.path(out_basepath,"Figures/Figure 2//RockmanAndI.png"),width=3,height=3,dpi=450)



dfPlot = rbind(data.frame(x="All tested\ngenes",y=c(sum(dfCisAndRockmanCis$significant|dfCisAndRockmanCis$rockmansig),
                                                    sum(!dfCisAndRockmanCis$significant&!dfCisAndRockmanCis$rockmansig)),leg = c("eQTL in at least\n one dataset","No hits"),pan="a"),
               data.frame(x="All significant\ngenes",y=c(enrichmentInRockmanRep[1,1],enrichmentInRockmanRep[1,2],enrichmentInRockmanRep[2,1]),
                          leg = c("Shared hits","Rockman et al only","Our dataset only"),pan="b"),
               data.frame(x="Mapped in \n1 cell type",y=c(enrichmentInRockmanMoreThanOne[1,2],enrichmentInRockmanMoreThanOne[2,2]),leg=c("Shared hits","Our dataset only"),pan="c"),
               data.frame(x="Mapped in \n>1 cell type",y=enrichmentInRockmanMoreThanOne[,1],leg=c("Shared hits","Our dataset only"),pan="c")
)

dfPlot["relY"] = dfPlot$y
dfPlot$leg = factor(dfPlot$leg,levels=unique(dfPlot$leg)[c(1,2,3,5,4)])
for(curX in dfPlot$x)
{
  dfPlot$relY[dfPlot$x==curX] = dfPlot$relY[dfPlot$x==curX]/sum(dfPlot$relY[dfPlot$x==curX])   
}
unique(dfPlot$leg)
col_pallete = viridisLite::viridis(end=0.6,n=4)[c(1,2,3,3,4)]
ggplot(dfPlot,aes(x=x,y=y,color=leg,fill=leg,pattern=leg,label=y)) + geom_bar_pattern(stat="identity",pattern_color=col_pallete[5],pattern_fill=col_pallete[5],position = "fill",width=0.4) +
  theme_classic() +   scale_fill_manual(breaks = unique(dfPlot$leg),values=col_pallete) + 
  scale_color_manual(breaks = unique(dfPlot$leg),values=col_pallete) + 
  scale_pattern_manual(breaks=unique(dfPlot$leg),values=c("None","None","stripe","None","None")) + facet_grid(~pan,scales="free_x") +
  geom_text(aes(y=relY),size = 3.5,fontface=2, position = position_stack(vjust = .5),color="white") + 
  ylab("proportion of genes") +xlab("") + theme(axis.text.x = element_text(angle=45,hjust=1,color="black"),  
                                                strip.background = element_blank(),
                                                strip.text.x = element_blank())

#ggsave(file.path(out_basepath,"Figures/Figure 2/rockman_and_us.svg"),width=5.8,height=3)
# Note, legend is removed from the final plot, and instead moved to the annotation of each bar for presentability
