#### This script generates the single cell eQTL matrix. All analysis should go here
#### All plotting should go in dedicated functions

source("analysis_functions/1_load_data_and_functions.R")


## single cell eQTL plotting function of raw data
plot_sneQTL = function(genoMatrix,neuCDS,neurDF,ind,normalize=T)
{
  Neuron = neurDF$Neuron[ind]
  gene = neurDF$transcript[ind]
  gene_id = neurDF$external_gene_id[ind]
  plot_df = rbind(data.frame(sample=Neuron,
                             Count = normalized_counts(neuCDS)[gene,colData(neuCDS)$fine_tissue == Neuron],
                             Genotype = genoMatrix[colData(neuCDS)$fine_tissue == Neuron,neurDF$lmarker[neurDF$transcript==gene][1]]),
                  data.frame(sample = "Other neurons",Count = normalized_counts(neuCDS)[gene,colData(neuCDS)$fine_tissue != Neuron],
                             Genotype = genoMatrix[colData(neuCDS)$fine_tissue != Neuron,neurDF$lmarker[neurDF$transcript==gene][1]]))
  if(!normalize)
    plot_df$Count = c(log(1+exprs(neuCDS)[gene,colData(neuCDS)$fine_tissue == Neuron]),log(1+exprs(neuCDS)[gene,colData(neuCDS)$fine_tissue != Neuron]))
  plot_df$sample = factor(plot_df$sample,levels=c(Neuron,"Other neurons"))
  textDF = data.frame(x=0,y=c(max(plot_df$Count)*1.15,max(plot_df$Count)*1.05),
                      sample = c(Neuron,"Other neurons"),label=c(sprintf("list(beta == %s,p == %s)",signif(neurDF$negbin.beta[ind]/neurDF$negbin.se[ind],2),signif(neurDF$FDR[ind],2)),
                                                                 sprintf("list(beta == %s,p == %s)",signif(neurDF$panNeuBetaNorm[ind],2),signif(neurDF$panNeuFDR[ind],2))))
  ggplot(plot_df,aes(x=Genotype,y=Count,color=sample,alpha=sample)) + geom_text(parse=T,data=textDF,aes(x=x,y=y,label=label),alpha=1,fontface="bold",hjust=0,size=3.3) + ylim(c(0,max(plot_df$Count)*1.15)) + geom_point(size=1)+ geom_point(data=subset(plot_df,sample==Neuron),size=1) +
    scale_color_manual(breaks = c(Neuron,"Other neurons"),values=c(viridisLite::viridis(6)[1],viridisLite::viridis(6)[4])) +  geom_smooth(method='lm', formula= y~x,se=F) +
    scale_alpha_manual(breaks = c(Neuron,"Other neurons"),values=c(1,0.2)) + theme_classic() + theme(legend.position=c(0.77,0.9),legend.title = element_blank()) + ylab(as.expression(bquote(italic(.(gene_id))~expression))) + xlab(expression("        "%<-% N2*"      "*CB4856%->%" "))
}




###### process single-neuron mapping to generate table of siginificant hits ######
#### extract significant hits and turn to nice data frame
neurList = lapply(fineTissueCis,function(x)x[x$FDR<0.1,])
neurList = neurList[sapply(neurList,nrow)>0]
neurList= lapply(neurList,function(x)cbind(x,geneAnnotsCur[match(x$transcript,geneAnnotsCur$wbps_gene_id),]))
neurDF = do.call(rbind,lapply(neurList,function(x){names(x)=names(neurList[[1]]);x}))
neurDF["Neuron"] = gsub("\\.+.*","",rownames(neurDF))
#### add info from pan neuronal mapping
neur_in_pan = match(neurDF$wbps_gene_id,tissueCis$Neuron$transcript)
neurDF = data.frame(neurDF,panNeuBetaNorm = tissueCis$Neuron$negbin.beta[neur_in_pan]/tissueCis$Neuron$negbin.se[neur_in_pan],
                    panNeuFDR = tissueCis$Neuron$FDR[neur_in_pan])

##### Differential expression in each cell with single-neuron eGene and all other neurons #####
neurDF["NeuronVsPan_estimate"] = NA
neurDF["NeuronVsPan_pvalue"] = NA
for(i in 1:nrow(neurDF))
{
  within = which(colData(cdsFinBC)$fine_tissue==neurDF$Neuron[i])
  without = which(colData(cdsFinBC)$broad_tissue=="Neuron"&colData(cdsFinBC)$fine_tissue!=neurDF$Neuron[i])
  cds_subset = cdsFinBC[neurDF$transcript[c(i)],c(within,without)]
  colData(cds_subset)$Neuron = c(rep(1,length(within)),rep(0,length(without)))
  gene_fits <- fit_models(cds_subset, model_formula_str = "~Neuron+Batch")
  fit_coefs <- coefficient_table(gene_fits)  %>% filter(term == "Neuron") %>%
    dplyr::select(gene_short_name, term, q_value, estimate)
  neurDF["NeuronVsPan_estimate"][i,1] = fit_coefs$estimate
  neurDF["NeuronVsPan_pvalue"][i,1] = fit_coefs$q_value
}
#### Enrichment in genes that are replicated across neurons in genes that are upregulated in the eqtl neuron compared to rest ######
fisher.test(cbind(c(sum(neurDF$NeuronVsPan_pvalue<0.05/267&neurDF$panNeuFDR<0.1),
                    sum(neurDF$NeuronVsPan_pvalue<0.05/267&neurDF$panNeuFDR>0.1)),
                  c(sum(neurDF$NeuronVsPan_pvalue>0.05/267&neurDF$panNeuFDR<0.1),
                    sum(neurDF$NeuronVsPan_pvalue>0.05/267&neurDF$panNeuFDR>0.1))))


###### Read leave-one-out analysis ####

looMappingRes = list()
for(neur in unique(neurDF$Neuron)) {
  cisNB = readRDS(sprintf("study_data//LOOmapping/out/tissueCisNB_%s_rds",colData(cdsFinBC)$fine_tissue[match(neur,colData(cdsFinBC)$fine_tissue)]))
  looMappingRes[[neur]]=cisNB
}
looMappingRes=lapply(looMappingRes, 
                      function(x) data.frame(t(sapply(x, function(y)  do.call('c', y[1:8] ))),
                                             stringsAsFactors=F))
looMappingRes=lapply(looMappingRes, function(x) {
  for( i in 3:8) {
    x[,i]=as.numeric(x[,i]) }
  return(x) })
looMappingRes=lapply(looMappingRes, function(x) {
  y=x
  y$FDR=p.adjust(y$negbin.p, method='fdr')
  y=na.omit(y)
  return(y)
})

neurDF["looBetaNorm"] = NA
neurDF["looPval"] = NA
for(i in 1:nrow(neurDF))
{
  neurDF[["looBetaNorm"]][i] = with(looMappingRes[[neurDF$Neuron[i]]],negbin.beta[transcript==neurDF$transcript[i]]/negbin.se[transcript==neurDF$transcript[i]])
  neurDF[["looPval"]][i] = with(looMappingRes[[neurDF$Neuron[i]]],FDR[transcript==neurDF$transcript[i]])
}  




saveRDS(neurDF,"study_data/singleNeuronCiseQTL.rds")

genoMatrixScEQTL = readRDS("study_data/geno_matrix_sceqtl.RDS")
