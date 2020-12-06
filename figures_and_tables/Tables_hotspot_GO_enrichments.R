source("analysis_functions/5_hotspot_enrichment_per_tissue.R")
source("figures_and_tables/table_formatting.R")
#### Generate all tissue eQTL table

tissueHotspotEnrichmentTable = NULL
tissueHotspotEnrichmentTargets = NULL
for(tis in names(tissueHotspotEnrichment))
{
  for(hs in names(tissueHotspotEnrichment[[tis]]))
  {
    curEnrichments = tissueHotspotEnrichment[[tis]][[hs]]$enrichments %>% filter(FDR<0.05)
    colnames(curEnrichments) = c("Ontology Type","GO.ID","Term","Annotated in tissue","Observed in eQTL","Expected","pvalue","Genes","FDR")
    curTargets  = tissueHotspotEnrichment[[tis]][[hs]]$targets %>% 
      mutate(normEffect = negbin.beta/negbin.se) %>%
      dplyr::select(c(`Type` = "type",`QTL chromosome` = "chr",     
                      `QTL peak position` = "pos",`CI_left` = "ciL",
                      `CI_right` = "ciR",`Transcript` = "transcript",`Normalized effect size` = "normEffect",
                      "FDR" = "FDR","Transcript chromosome" = "tchr",
                      `Start` = "start_position",`End` = "end_position",
                      `Strand` = "strand", `external_gene_id`,
                      `external_transcript_id`,`wbps_transcript_id`))
    tissueHotspotEnrichmentTargets = rbind(tissueHotspotEnrichmentTargets,
                                           data.frame(Hotspot=sprintf("%s:%s",tis,hs),curTargets,stringsAsFactors = F)
                                           )
    
  if(nrow(curEnrichments)>0)
  {
    
   tissueHotspotEnrichmentTable = rbind(tissueHotspotEnrichmentTable,
                                           data.frame(Hotspot=sprintf("%s:%s",tis,hs),curEnrichments,stringsAsFactors = F)
    )
  }
}
}

wb = make_excel_workbook(list("Enrichment GO terms in hotspots" = tissueHotspotEnrichmentTable,
                              "Tissue-specific hotspot targets" = tissueHotspotEnrichmentTargets)) 
#openxlsx::saveWorkbook(wb,file=file.path(out_basepath,"Tables/tissue_hotspots_targets.xlsx"),overwrite=T)
