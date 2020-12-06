source("analysis_functions//5_hotspot_enrichment_per_tissue.R")

k = 0.3
for(tis in c("Body Wall Muscle","Neuron","Seam Cells","Intestine"))
{
  curPlot = generatePlotTissue(tissuePeaks[[tis]],cisPeaks = tissueCis[[tis]],dfReturn = T,cutoff=0.1)
  curPlot["color"] = rep("black",nrow(curPlot))
  vlines = NULL
  for(hs in names(tissueHotspotEnrichment[[tis]]))
  {
  cur_inHS = rownames(plyr::match_df(curPlot,tissueHotspotEnrichment[[tis]][[hs]]$targets))
  curPlot[cur_inHS,"color"] = viridis::viridis(1,begin=k,end=k)
  k = k + 0.1
  }
  curPlot["isblack"] = (curPlot$color == "black")
  pl = ggplot(curPlot,aes(x=pos,y=start_position,color=color,alpha=color)) + scale_color_identity()+ xlab("QTL position (mb)") + ylab("transcript position (mb)") + 
    scale_x_continuous(breaks = c(5e6,1.5e7),labels = function(x)x/1e6) + 
    scale_y_continuous(breaks = c(5e6,1.5e7),labels = function(x)x/1e6) + 
    vlines + 
    facet_grid(chromosome_name~chr,scales="free") + geom_point(data=curPlot %>% filter(!isblack),alpha=1) + 
    geom_point(data =curPlot %>% filter(isblack),alpha=0.4) + theme_classic() + 
    geom_errorbarh(aes(xmin=ciL,xmax=ciR,y=start_position,alpha=-log10(FDR)),alpha=0.3) + 
    theme(legend.position = "none")
  for(hs in names(tissueHotspotEnrichment[[tis]]))
  {
    hs_spot = unlist(strsplit(hs,"_"))
    pl = pl + 
    geom_vline(data=data.frame(chr=hs_spot[1],xintercept=as.numeric(hs_spot[2])),
               aes(xintercept=xintercept,chr=chr),linetype=3,alpha=0.5)
  }
  print(pl)
#  ggsave(file.path(out_basepath,sprintf("Figures/Supplementary Figures/NeuronalHotspot/%s.png",tis)),dpi=600,width=4,height=4)
}
