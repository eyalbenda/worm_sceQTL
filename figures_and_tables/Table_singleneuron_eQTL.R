source("analysis_functions/7_single_neuron_GO_enrichment.R")
source("figures_and_tables/table_formatting.R")
#### Generate all tissue eQTL table

curNeuron = singleNeuronEnrichment$enrichments %>% filter(FDR<0.05)
colnames(curNeuron) = c("Ontology Type","GO.ID","Term","Annotated in tissue","Observed in eQTL","Expected","pvalue","Genes","FDR")
curTargets  = neurDF %>% 
  mutate(normEffect = negbin.beta/negbin.se,Type="cis") %>%
  dplyr::select(c(`Neuron`="Neuron",`Transcript` = "transcript",`Normalized effect size` = "normEffect",
                  "FDR" = "FDR",`Pan neuronal normalized effect size` = "panNeuBetaNorm",
                  `pan neuronal FDR` = "panNeuFDR",`all neurons but eQTL neuron normalized effect size` ="looBetaNorm",`all neurons but eQTL neuron FDR` = "looPval", 
                  `Transcript chromosome` = "chromosome_name",
                  `Start` = "start_position",`End` = "end_position",
                  `Strand` = "strand", `external_gene_id`,
                  `external_transcript_id`,`wbps_transcript_id`))

wb = make_excel_workbook(list("single-neuronal subtype eQTL" = curTargets,
                              "GO enrichments in sn-eQTL" = curNeuron)) 
openxlsx::saveWorkbook(wb,file=file.path(out_basepath,"Tables/singleNeuroneQTL.xlsx"),overwrite=T)
# wb <- createWorkbook()
# addWorksheet(wb, "Enriched GO terms in hotspots")
# header_st <- createStyle(textDecoration = "Bold")
# writeData(wb,"Enriched GO terms in hotspots","Table S#. Enriched GO terms in tissue specific hotspots",startCol = 1,startRow = 1,headerStyle = header_st)
# writeData(wb,"Enriched GO terms in hotspots",tissueHotspotEnrichmentTable,startRow = 2,headerStyle = header_st)
# openxlsx::saveWorkbook(wb,file="~/hujidrive/Work/Research/Postdoc/scRNAseq/Manuscript/Tables/tissue_hotspots_table.xlsx",overwrite=T)
# 
# addWorksheet(wb, "Tissue-specific hotspot targets")
# writeData(wb,"Tissue-specific hotspot targets","Table S#. Targets of tissue-specific hotspots",startCol = 1,startRow = 1,headerStyle = header_st)
# writeData(wb,"Tissue-specific hotspot targets",tissueHotspotEnrichmentTargets,headerStyle = header_st,startRow = 2)
# openxlsx::saveWorkbook(wb,file="~/hujidrive/Work/Research/Postdoc/scRNAseq/Manuscript/Tables/tissue_hotspots_targets.xlsx",overwrite=T)
