source("analysis_functions/1_load_data_and_functions.R")
singleNeuronCis = readRDS("study_data/singleNeuronCiseQTL.rds")


  


# listMarts(host = "parasite.wormbase.org")
# wormbase = useMart(biomart = "parasite_mart",
#                    host = "https://parasite.wormbase.org",
#                    port = 443)
# wormbase = useDataset(mart = wormbase, dataset = "wbps_gene")
# head(listFilters(wormbase))
# listAttributes(wormbase)[grep("go",(listAttributes(wormbase))$name),]
# listAttributes(wormbase)[grep("entrez",(listAttributes(wormbase))$name),]
# go_term_worm = getBM(attributes= c("wbps_gene_id","go_accession"),
#                      filters="wbps_gene_id", values=unique(geneAnnotsCur$wbps_gene_id), mart=wormbase)
#### go terms retrieved on April 21 2020
go_term_worm = readRDS("study_data/go_term_worm.rds")
go_term_worm_l = split(go_term_worm, f=go_term_worm$go_accession)
go_term_list = lapply(go_term_worm_l, function(x){x$wbps_gene_id})
# Go term list to use with Topgo
go_term_list_inv = inverseList(go_term_list)
genes = names(go_term_list_inv)
# Setup vector
relevant_genes = colnames(isExpressed[,isExpressed["Neuron",]])
all_genes = as.factor(rep(1,length(relevant_genes)))
levels(all_genes) = c("0","1")
names(all_genes) = relevant_genes
# select random genes to test for enrichment
real_list = singleNeuronCis$transcript
# real_list = tissueHotspot[[tis]]$hotspot.linkages[[hs]]$transcript
g1 = factor(as.numeric((names(all_genes))%in% real_list))
names(g1) = names(all_genes)
### make topGo object
GODataBP = new("topGOdata", ontology="BP",allGenes=g1, annot=annFUN.gene2GO, gene2GO=go_term_list_inv)
GODataMF = new("topGOdata", ontology="MF",allGenes=g1, annot=annFUN.gene2GO, gene2GO=go_term_list_inv)
GODataCC = new("topGOdata", ontology="CC",allGenes=g1, annot=annFUN.gene2GO, gene2GO=go_term_list_inv)
## test real genes
test.stat =new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
## Run test
resultFisBP <- runTest(GODataBP, algorithm = "classic", statistic = "fisher")
resultFisMF <- runTest(GODataMF, algorithm = "classic", statistic = "fisher")
resultFisCC <- runTest(GODataCC, algorithm = "classic", statistic = "fisher")
# Table of results ordered by test statistc
allResBP <- GenTable(GODataBP, classic = resultFisBP, orderBy = "weight", ranksOf = "classic", topNodes = 300)
allResMF <- GenTable(GODataMF, classic = resultFisMF, orderBy = "weight", ranksOf = "classic", topNodes = 300)
allResCC <- GenTable(GODataCC, classic = resultFisCC, orderBy = "weight", ranksOf = "classic", topNodes = 300)
combined_table = rbind(data.frame(type = "Biological Process",allResBP),
                       data.frame(type = "Molecular Function",allResMF),
                       data.frame(type = "Cellular Component",allResCC)) %>% mutate(pvalue = as.numeric(classic)) %>% arrange(pvalue) %>% dplyr::select(-classic)
genesBP = sapply(combined_table$GO.ID,function(x){genes = unlist(genesInTerm(object = GODataBP,whichGO = x));genes_rel = genes[genes %in% as.character(real_list)];paste(geneAnnotsCur$external_gene_id[match(genes_rel,geneAnnotsCur$wbps_gene_id)],collapse=";")})
genesMF = sapply(combined_table$GO.ID,function(x){genes = unlist(genesInTerm(object = GODataMF,whichGO = x));genes_rel = genes[genes %in% as.character(real_list)];paste(geneAnnotsCur$external_gene_id[match(genes_rel,geneAnnotsCur$wbps_gene_id)],collapse=";")})
genesCC = sapply(combined_table$GO.ID,function(x){genes = unlist(genesInTerm(object = GODataCC,whichGO = x));genes_rel = genes[genes %in% as.character(real_list)];paste(geneAnnotsCur$external_gene_id[match(genes_rel,geneAnnotsCur$wbps_gene_id)],collapse=";")})
combined_table["genes"] = sapply(combined_table$GO.ID,function(x){genes = unlist(genesInTerm(object = GODataBP,whichGO = x));genes_rel = genes[genes %in% as.character(real_list)];paste(geneAnnotsCur$external_gene_id[match(genes_rel,geneAnnotsCur$wbps_gene_id)],collapse=";")})
combined_table[["genes"]][combined_table$type=="Molecular Function"] = genesMF[combined_table$type=="Molecular Function"]
combined_table[["genes"]][combined_table$type=="Cellular Component"] = genesCC[combined_table$type=="Cellular Component"]
combined_table["FDR"] = p.adjust(combined_table$pvalue,method="fdr")
combined_table = combined_table[order(combined_table$pvalue),]
singleNeuronEnrichment = list(targets = singleNeuronCis,
                              enrichments = combined_table)


saveRDS(singleNeuronEnrichment,"study_data/singleNeuronEnrichment.rds")
