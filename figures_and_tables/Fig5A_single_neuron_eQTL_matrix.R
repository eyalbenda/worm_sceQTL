source("analysis_functions/6_single_neuron_eQTL_processing.R")

neuron_df_plot = neurDF %>% arrange(Neuron)
sneqtl_beta = matrix(data = NA,ncol=length(unique(neuron_df_plot$external_gene_id)),nrow=length(unique(neuron_df_plot$Neuron))+3)
rownames(sneqtl_beta) = c("Pan neuronal 10% FDR","Pan neuronal 50% FDR","Pan neuronal sign",unique(neuron_df_plot$Neuron))
colnames(sneqtl_beta) = unique(neuron_df_plot$external_gene_id)
sneqtl_beta[1:2,] = 0
for(i in 1:nrow(neuron_df_plot))
{
  if(neuron_df_plot$panNeuFDR[i]<0.1) sneqtl_beta[1,neuron_df_plot$external_gene_id[i]] = 2
  if(neuron_df_plot$panNeuFDR[i]<0.5) sneqtl_beta[2,neuron_df_plot$external_gene_id[i]] = 2
  sneqtl_beta[3,neuron_df_plot$external_gene_id[i]] = (sign(neuron_df_plot$panNeuBetaNorm[i])==(sign(neuron_df_plot$negbin.beta[i])))*2
  sneqtl_beta[neuron_df_plot$Neuron[i],neuron_df_plot$external_gene_id[i]] = 1
  
}
sneqtl_beta = t(data.frame(t(sneqtl_beta)) %>% arrange(-.[[1]],-.[[2]]))
colnames(sneqtl_beta) = unique(neuron_df_plot$external_gene_id)

conflictingSign = 
  unique(neurDF$external_gene_id[neurDF$transcript %in% names(which(lapply(tapply(sign(neurDF$negbin.beta),neurDF$transcript,table),length)>1))])
sneqtl_beta[3,conflictingSign] = 3
row_split = c(rep("",3),rep("",nrow(sneqtl_beta)-3))
rownames(sneqtl_beta)[1:3] = c("10% FDR","50% FDR","same effect direction")


#pdf(file.path(out_basepath,"Figures/Figure 3/Fig3A.heatmap.pdf"),width=4.5,height=8)
Heatmap(sneqtl_beta,row_names_side = "left",row_names_gp = gpar(fontsize=9),
        cluster_rows = F,
        cluster_columns = F,border = T,
        show_column_names = F,na_col = "white",col = circlize::colorRamp2(c(0,1,2,3),c("white",viridisLite::viridis(6)[4:5],"grey")),
        show_heatmap_legend = F,row_split = row_split)
#dev.off()
