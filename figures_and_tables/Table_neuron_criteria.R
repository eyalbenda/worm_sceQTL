source("analysis_functions/1_load_data_and_functions.R")
source("figures_and_tables/table_formatting.R")
#### This script repeats the neuronal classification from to provide reproducibility
#### since original was done over several iteration
uniqueifyGenes = function(...)
{
  paste(unique(unlist(strsplit(...,";"))),collapse=";")
}
whichClust = function(cds,gene,cutoff=0.1)
{
  if(!gene %in% rowData(cds)$external_gene_id)
    return(integer(0))
  which(tapply(normalized_counts(cds)[rowData(cds)$external_gene_id %in% c(gene),],clusters(cds),function(x)sum(x>0)/length(x))>cutoff)
}
intersectGenes = function(cds,genes,cutoff = 0.1)
{
  i = 0
  inter = clusters(cds)
  while(i<length(genes))
  {
    clusts = whichClust(cds,genes[i+1],cutoff)
    print(sprintf("gene %s has clusters %s above cutoff",genes[i+1],paste(clusts,collapse=" ")))
    inter = intersect(inter,clusts)
    i = i+1
  }
  print(inter)
  i = 0
  geneOut=""
  if(length(inter)==0)
    return(NULL);
  for(clust in 1:length(inter))
  {
    while(i<length(genes))
    {
      if(inter[clust] %in% whichClust(cds,genes[i+1],cutoff))
        geneOut = paste0(geneOut,genes[i+1],";")
      i = i+1
    }
    names(inter)[clust] = geneOut
  }
  inter
}



cellmarkers = readxl::read_excel("external_data/Supplementary_Tables_190611.xlsx",sheet = 11,skip = 10)
neuronmarker = cellmarkers %>% filter(UMAP=="Neurons")
neuronmarkerParsed = sapply(strsplit(gsub(",","",gsub("\r\n",", ",neuronmarker$`Marker genes and annotation notes`))," "),function(x)x[x %in% geneAnnotsCur$external_gene_id])

ourAnnotation = as.numeric(unique(clusters(neuCDS)))
names(ourAnnotation)=colData(neuCDS)$fine_tissue[match(unique(clusters(neuCDS)),clusters(neuCDS))]
ourAnnotation = c(ourAnnotation,"AWB"=94,"M1"=95)
neuronAnnot = data.frame(Cluster=sort(ourAnnotation),Neuron=names(sort(ourAnnotation)),genes=NA,stringsAsFactors = F)

programmaticParse= sapply(neuronmarkerParsed,function(x)intersectGenes(neuCDS,x,0.05))
names(programmaticParse) = neuronmarker$`Cell type`
programmaticParse
singleSupportingCao = programmaticParse[which((sapply(programmaticParse,length)>0))]
CaoNeuron = rep(names(singleSupportingCao),times=sapply(singleSupportingCao,length))
singleSupportingCao2 = singleSupportingCao
names(singleSupportingCao2) = NULL
CaoNeuronDF = data.frame(Cluster = unlist(singleSupportingCao2),
                         Neuron = CaoNeuron,
                         Genes = names(unlist(singleSupportingCao2)),stringsAsFactors = F)
CaoNeuronDFMatch = CaoNeuronDF[paste(CaoNeuronDF$Cluster,CaoNeuronDF$Neuron) %in% paste(neuronAnnot$Cluster,neuronAnnot$Neuron),]
neuronAnnot[match(CaoNeuronDFMatch$Cluster,neuronAnnot$Cluster),3] = CaoNeuronDFMatch$Genes

cellmarkersCENGeN =readxl::read_xlsx("external_data//CENgen.xlsx",sheet = 3,skip=2)
neuronmarkerParsedCENGeN = apply(cellmarkersCENGeN[,-c(1:4)],1,function(x)paste(x[!is.na(x)]))
cellmarkersCENGeN = cellmarkersCENGeN %>% filter(sapply(neuronmarkerParsedCENGeN,length)>0)
neuronmarkerParsedCENGeN = apply(cellmarkersCENGeN[,-c(1:4)],1,function(x)paste(x[!is.na(x)]))

programmaticParseCENGeN= sapply(neuronmarkerParsedCENGeN,function(x)intersectGenes(neuCDS,x,0.05))
CENGeNunique = rep(cellmarkersCENGeN$`Neuron Class`,times=sapply(programmaticParseCENGeN,length))
CENGeNuniqueParsed = unlist(programmaticParseCENGeN[sapply(programmaticParseCENGeN,length)>0])
CENGeNuniqueParsedDF = data.frame(cluster=CENGeNuniqueParsed,neuron=CENGeNunique,genes=names(CENGeNuniqueParsed),stringsAsFactors = F)
CENGeNuniqueParsedDFMatch = CENGeNuniqueParsedDF[paste(CENGeNuniqueParsedDF$cluster,CENGeNuniqueParsedDF$neuron) %in% paste(neuronAnnot$Cluster,neuronAnnot$Neuron),]
neuronAnnot[match(CENGeNuniqueParsedDFMatch$cluster,neuronAnnot$Cluster),3] = CENGeNuniqueParsedDFMatch$genes
neuronAnnot[1:10,]
## markers for all of these neurons support same clusters
for(clust in c(84,88))
{
  curGenes = uniqueifyGenes(as.character(unlist(CENGeNuniqueParsedDF$genes[CENGeNuniqueParsedDF$cluster==84],
                                               CENGeNuniqueParsedDF$genes[CENGeNuniqueParsedDF$cluster==88]))) 
  neuronAnnot[clust,2] = "ALM_PLM_AVM_PVM"
  neuronAnnot[clust,1] = clust
  neuronAnnot[clust,3] = curGenes
}
## New neuron in our dataset, CEM
neuronAnnot[neuronAnnot[,2]=="CEM",3] = "cwp-1"
## OLL and URY genes (except for prdf-1 which doesn't exist) both support cluster 23
neuronAnnot[neuronAnnot[,2]=="OLL_URY",3] = c("osm-5;sox-2;hecw-1;ric-4;ceh-32;osm-6;vab-3;ser-2;str-94;F35F11.3;hsd-2;unc-86;vab-3;pros-1;pdrf-1;glr-4;glr-5;tol-1;lim-7")
### CeNGEN markers find  IL2_DV and 61 IL2_LR, unique gene egas-3 (DV) and 	egas-4 (LR) support the assignment
neuronAnnot[neuronAnnot[,2]=="IL2_DV",3] = "cha-1;cho-1;klp-6;tba-6;unc-86;egas-3"
neuronAnnot[neuronAnnot[,2]=="IL2_LR",3] = "cha-1;cho-1;klp-6;tba-6;unc-86;unc-39;egas-4"
neuronAnnot[neuronAnnot[,2]=="ASH",3] = "osm-10;T23B3.6;srm-4"
neuronAnnot[neuronAnnot[,1] %in% c(6,46),3] = uniqueifyGenes(CENGeNuniqueParsedDFMatch$genes[CENGeNuniqueParsedDFMatch$neuron=="VA"])
neuronAnnot[neuronAnnot[,2]=="AIA",3] = "unc-39;glc-3;flp-2;"
neuronAnnot[neuronAnnot[,2]=="ADA",3] = "C10C5.7;flp-26;T05A8.3;flp-7;flp-16"
neuronAnnot[neuronAnnot[,2]=="RMG",3] = "nlp-46;flp-14;flp-21;C44B11.4;R13A1.5"
neuronAnnot[neuronAnnot[,2]=="ASI_ASJ",3] = "cha-1;cho-1;trx-1;ssu-1;ins-6;ins-3;R102.2"
neuronAnnot[neuronAnnot[,2]=="MI",3] ="ceh-45;nlp-3;glr-8;alr-1;ceh-34;R07E3.2;flp-6;C25G6.4"
neuronAnnot[neuronAnnot[,2]=="RME",3] = "ceh-27;unc-46;unc-47;ceh-32;avr-15;slt-1;snf-11"
neuronAnnot[neuronAnnot[,2]=="Unknown_ACh_3",3] = "nlp-5;R11D1.12;lag-1;F31B12.4;flp-18"
neuronAnnot[neuronAnnot[,2]=="BDU",3] = "gcy-35;flp-12lahr-1;C12D5.5;K02B12.9"
neuronAnnot[neuronAnnot[,2]=="Unknown_glut_2",3] = "vab-7;C09G5.13;ZC196.9;mab-23;C50D2.6"
neuronAnnot[neuronAnnot[,2]=="DA",3] = "lron-14;hecw-1;lin-39;twk-40;C13A2.12;R10E11.5"
neuronAnnot[neuronAnnot[,2]=="PDE",3] = "F09E10.7;F59F3.6;dat-1;cat-2;flp-33"
neuronAnnot[neuronAnnot[,2]=="I2_I3",3] = uniqueifyGenes(names(unlist(programmaticParseCENGeN[cellmarkersCENGeN$`Neuron Class` %in% c("I3")])),names(unlist(programmaticParseCENGeN[cellmarkersCENGeN$`Neuron Class` %in% c("I2")])))
neuronAnnot[neuronAnnot[,2]=="Unknown_touch",3] = "mec-12;mec-7"
neuronAnnot[neuronAnnot[,2]=="AUA",3] = "flp-8;ceh-6;C54C6.5;best-23;E02A10.4"
neuronAnnot[neuronAnnot[,2]=="AWB",3] = "srd-23;srab-24;odr-1;daf-11;sox-2"
neuronAnnot[neuronAnnot[,2]=="M1",3] = "ceh-34;flp-5;nlp-13;flp-26;ceh-32;flr-2;glr-2"

neuronAnnotOut = data.frame(Neuron = unique(colData(neuCDS)$fine_tissue),Genes = neuronAnnot[match(unique(colData(neuCDS)$fine_tissue),neuronAnnot$Neuron),3])

neuronAnnotOut = neuronAnnotOut %>% mutate("No. of cells" = table(colData(neuCDS)$fine_tissue)[neuronAnnotOut$Neuron]) %>% arrange(Neuron)

neuron_criteria_table = make_excel_workbook(table_list = list("single neuron annotation" = neuronAnnotOut))
openxlsx::saveWorkbook(neuron_criteria_table,file = file.path(out_basepath,"Tables/Neuron_classification_criteria.xlsx"),overwrite=T)


saveRDS(neuronAnnotOut,"study_data/neuronAnnotationTable.rds")
