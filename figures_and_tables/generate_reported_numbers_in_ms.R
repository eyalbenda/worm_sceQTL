source("analysis_functions/1_load_data_and_functions.R")

###### Cell-type-specific eQTL mapping in a single one-pot experiment #####
## correlation expected vs. observed cell tyes (see generation of figure S1 for more detail)
source("figures_and_tables/FigS1_observed_v_expected_cells.R")
cor.test(table(colData(cdsFinBC)$broad_tissue),
    expected,
    method = "spearman",
    use = "pairwise.complete.obs")
rm(list=ls(all=T))

## eQTL mapping in multiple cell types (see Table_all_eQTL_broad for more detail on aggregating table)
source("figures_and_tables/Table_all_eQTL_broad.R")
source("figures_and_tables/Table_summary_per_celltype.R")
# correlation between cell types
cor.test(stats_per_celltype$`cis-eQTL`+stats_per_celltype$`trans-eQTL`,stats_per_celltype$`No. of cells`)
      # number of cis and trans eQTLs
table(all_eqtl_table$Type)
# number of genes with cis and trans eQTLs
tapply(all_eqtl_table$Transcript,all_eqtl_table$Type,function(x)length(unique(x)))
# number of genes with one or more than one eQTL
ciseQTL = data.frame(all_eqtl_table %>% filter(Type=="cis"))
sum(table(ciseQTL$Transcript)==1)
sum(table(ciseQTL$Transcript)>1)
multi_mapping = names(which(table(ciseQTL$Transcript)>1))
multi_mapping_eQTL = ciseQTL %>% filter(Transcript %in% multi_mapping)
number_of_directions_mapped = tapply(sign(multi_mapping_eQTL$Normalized.effect.size),as.character(multi_mapping_eQTL$Transcript),function(x)length(unique(x)))
# Number of genes with only concordant signs
sum(number_of_directions_mapped==1)
# Number of genes with any discordant signs
sum(number_of_directions_mapped>1)
rm(list=ls(all=T))

## comparison between cell type eQTL and parental data
source("analysis_functions/3_comparison_to_parental_dataset.R") ## See that file for parsing parental differential expression data and merging with our eQTL data
# number of differentially expressed genes in the parents
sum(isDifExpressed)
# number of differential expressed genes with cis eQTL in our dataset
# odds ratio and p-value for enriched overlap
contingency = cbind(c(sum(isDifExpressed&dfCisAndParents$FDR<0.1),
                      sum(isDifExpressed&dfCisAndParents$FDR>0.1)),
                    c(sum(!isDifExpressed&dfCisAndParents$FDR<0.1),
                      sum(!isDifExpressed&dfCisAndParents$FDR>0.1)))
fisher.test(contingency)
# how many have the same sign?
with(dfCisAndParents %>% filter(FDR<0.1&isDifExpressed),sum(sign(log2(norm_foldChange))==sign(negbin.beta)))
# percentage:
with(dfCisAndParents %>% filter(FDR<0.1&isDifExpressed),sum(sign(log2(norm_foldChange))==sign(negbin.beta)))/sum(isDifExpressed&dfCisAndParents$FDR<0.1)
rm(list=ls(all=T))


###### Comparison between single cell and bulk (rockman et al 2010) data #####
source("analysis_functions/4_comparison_to_wholeworm_data.R")
# number of genes with information in both datasets
dim(dfCisAndRockmanCis)[1]
# number of eQTL in rockman dataset
sum(dfCisAndRockmanCis$rockmansig)
# number of shared eQTL
sum(dfCisAndRockmanCis$significant&dfCisAndRockmanCis$rockmansig)
# derive odds ratio
enrichmentInRockmanRep = cbind(c(sum(dfCisAndRockmanCis$significant&dfCisAndRockmanCis$rockmansig),
                                 sum(dfCisAndRockmanCis$significant&!dfCisAndRockmanCis$rockmansig)),
                               c(sum(!dfCisAndRockmanCis$significant&dfCisAndRockmanCis$rockmansig),
                                 sum(!dfCisAndRockmanCis$significant&!dfCisAndRockmanCis$rockmansig)))
fisher.test(enrichmentInRockmanRep)
# percentage of shared from all rockman
sum(dfCisAndRockmanCis$significant&dfCisAndRockmanCis$rockmansig)/sum(dfCisAndRockmanCis$rockmansig)
# percentage of shared for our dataset
sum(dfCisAndRockmanCis$significant&dfCisAndRockmanCis$rockmansig)/sum(dfCisAndRockmanCis$significant)
# correlation between data in rockman for all significant cis eQTLs
corDF = with(dfCisAndRockmanCis %>% filter(significant),makeCorText(previous,ourdata))
# percentage of genes with mutliple cis eqtl replicated in the bulk study
sum(dfCisAndRockmanCis$SumSignificantInOut>1&dfCisAndRockmanCis$rockmansig)/sum(dfCisAndRockmanCis$SumSignificantInOut>1)
# percentage of the genes with only one cis eqtl replicated in the bulk study
sum(dfCisAndRockmanCis$SumSignificantInOut==1&dfCisAndRockmanCis$rockmansig)/sum(dfCisAndRockmanCis$SumSignificantInOut==1)
# fisher test
enrichmentInRockmanMoreThanOne = cbind(c(sum(dfCisAndRockmanCis$SumSignificantInOut>1&dfCisAndRockmanCis$rockmansig),
                                         sum(dfCisAndRockmanCis$SumSignificantInOut>1&!dfCisAndRockmanCis$rockmansig)),
                                       c(sum(dfCisAndRockmanCis$SumSignificantInOut==1&dfCisAndRockmanCis$rockmansig),
                                         sum(dfCisAndRockmanCis$SumSignificantInOut==1&!dfCisAndRockmanCis$rockmansig)))
fisher.test(enrichmentInRockmanMoreThanOne)
rm(list=ls(all=T))

######## Shared and cell-type specific trans-eQTL hotspots ########

source("figures_and_tables/Tables_hotspot_GO_enrichments.R")
# how many linkages fall in the hotspots
sum(tissueHotspotEnrichmentTargets$FDR<0.1)
# get number of eQTLs in FDR 0.1
table(tissueHotspotEnrichmentTargets$Hotspot[tissueHotspotEnrichmentTargets$FDR<0.1])
# neuron and seam cells same eqtl sums up to 31
# rinse and repeat for 0.2 FDR
table(tissueHotspotEnrichmentTargets$Hotspot[tissueHotspotEnrichmentTargets$FDR<0.2])
# neuron and seam cells same eqtl sums up to 42

# enrichment terms reported
tissueHotspotEnrichmentTable %>% filter(Term=="innate immune response")
tissueHotspotEnrichmentTable %>% filter(Term=="myofilament")
tissueHotspotEnrichmentTable %>% filter(Term=="actin cytoskeleton")
tissueHotspotEnrichmentTable %>% filter(Term=="vesicle localization")
tissueHotspotEnrichmentTable %>% filter(Term=="BMP receptor binding")

rm(list=ls(all=T))



######## Cell-specific eQTL effects in the C. elegans nervous system #########
source("analysis_functions/6_single_neuron_eQTL_processing.R")
# the number of neurons
dim(neuCDS)[2]
# number of distinct neurons, minimum and maximum number of cells
length(unique(colData(neuCDS)$fine_tissue))
min(table(colData(neuCDS)$fine_tissue))
max(table(colData(neuCDS)$fine_tissue))
# number of sn-eQTL at FDR<0.1
dim(neurDF)[1]
# number of sn-eQTL genes
length(unique(neurDF$transcript))
# number of sn-eQTL in a single subtype
sum(tapply(neurDF$Neuron,neurDF$transcript,length)==1)
source("analysis_functions/7_single_neuron_GO_enrichment.R")
singleNeuronEnrichment$enrichments %>% filter(Term=="signaling")
singleNeuronEnrichment$enrichments %>% filter(grepl("G protein-coupled receptor signaling",Term))
singleNeuronEnrichment$enrichments %>% filter(Term=="neuropeptide signaling pathway")
# number of sn-eQTL not replicating at 0.1
sum(neurDF$panNeuFDR>0.1)
# number of sn-eQTL not replicating at 0.5
sum(neurDF$panNeuFDR>0.5)
# number of sn-eQTL not replicating at 0.5 and same beta direction without neuron
sum(neurDF$panNeuFDR>0.5&sign(neurDF$looBetaNorm)==sign(neurDF$negbin.beta))
sum(neurDF$panNeuFDR>0.5&sign(neurDF$looBetaNorm)!=sign(neurDF$negbin.beta))
# number of sn-eQTL replicating at 0.5 and same beta direction without neuron
sum(neurDF$panNeuFDR<0.5&sign(neurDF$looBetaNorm)==sign(neurDF$negbin.beta))
sum(neurDF$panNeuFDR<0.5&sign(neurDF$looBetaNorm)!=sign(neurDF$negbin.beta))
