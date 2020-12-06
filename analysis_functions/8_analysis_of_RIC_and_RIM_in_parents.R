source("analysis_functions/6_single_neuron_eQTL_processing.R")

subsetProcess=function(cds,num_dim=50,k=20)
{cds <- preprocess_cds(cds, num_dim = num_dim)
cds <- reduce_dimension(cds)
cluster_cells(cds,k=k)}
histClust = function(cds,gene)
{
  hist(tapply(normalized_counts(cds)[rowData(cds)$external_gene_id %in% c(gene),],clusters(cds),function(x)sum(x>0)/length(x)))
}
whichClust = function(cds,gene,cutoff=0.1)
{
  which(tapply(normalized_counts(cds)[rowData(cds)$external_gene_id %in% c(gene),],clusters(cds),function(x)sum(x>0)/length(x))>cutoff)
}

cdsL2 = readRDS("study_data/L2monocle_parents.rds")
annot_L2 = read.csv("study_data/parents_L2_classification.csv")
neuCDSparents = subsetProcess(cdsL2[,which(annot_L2$cluster[match(colnames(cdsL2),annot_L2$barcode)]=="Neurons")],num_dim = 100,k = 3)
length(unique(clusters(neuCDSparents)))

## tbh-1;flp-22;cat-1;tdc-1;lin-11;unc-62 converge on 26 (RIC)
## tdc-1 and nmr-2 converge on 21 (RIM), nmr-1 also but very lowly expressed


