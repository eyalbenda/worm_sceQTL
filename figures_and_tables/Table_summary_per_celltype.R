source("analysis_functions///1_load_data_and_functions.R")
source("figures_and_tables//table_formatting.R")
#### Read in genotype counts
refCounts = readRDS("study_data/refCounts_backgroundRemovedstrict.rds")
altCounts = readRDS("study_data/altCounts_backgroundRemovedstrict.rds")


stats_per_celltype = data.frame(colData(cdsFinBC)) %>% dplyr::group_by(`Cell-type` = broad_tissue) %>% summarize(`No. of cells` = n())
stats_per_celltype["Median No. UMI"] = data.frame(a=colSums(exprs(cdsFinBC))) %>% group_by(colData(cdsFinBC)$broad_tissue) %>% summarize(med=median(a)) %>% pull(med) 
stats_per_celltype["Median No. reads in SNVs"] = data.frame(a = colSums((altCounts[,colData(cdsFinBC)$barcode] + refCounts[,colData(cdsFinBC)$barcode]))) %>% 
  group_by(colData(cdsFinBC)$broad_tissue) %>% summarize(med = median(a)) %>% pull(med)
stats_per_celltype["Median No. SNVs with reads"] = data.frame(a = colSums((altCounts[,colData(cdsFinBC)$barcode]>0|refCounts[,colData(cdsFinBC)$barcode]>0))) %>% 
  group_by(colData(cdsFinBC)$broad_tissue) %>% summarize(med = median(a)) %>% pull(med)
stats_per_celltype["cis-eQTL"] = sapply(stats_per_celltype$`Cell-type`,function(x)
  generatePlotTissue(tissuePeaks[[x]],tissueCis[[x]],dfReturn = T,0.1) %>% filter(type=="cis") %>% nrow)
stats_per_celltype["trans-eQTL"] = sapply(stats_per_celltype$`Cell-type`,function(x)
  generatePlotTissue(tissuePeaks[[x]],tissueCis[[x]],dfReturn = T,0.1) %>% filter(type=="trans") %>% nrow)


#### Estimate number of individuals (number of cells with correlation > 0.9) ####
# genotypes_all = readRDS(path_to_genos_matrix) # Too heavy to be included with the package, please write eyal.bendavid@mail.huji.ac.il and we'll share it (sorry!!!)
# over90 = NULL
# for(celltype in stats_per_celltype$`Cell-type`)
# {
#   print(celltype)
#   tmp = cor(t(genotypes_all[colData(cdsFinBC)$broad_tissue==celltype,]))
#   diag(tmp) = NA
#   over90[celltype] = sum(apply(tmp,2,max,na.rm=T)>0.9)
#   rm(tmp)
# }

# stats_per_celltype["Estimated proportion of collisions"] = round(over90[match(stats_per_celltype$`Cell-type`,names(over90))]/stats_per_celltype$`No. of cells`,2)

all_cell_type_table = make_excel_workbook(list("Cell-type statistics"=stats_per_celltype))
openxlsx::saveWorkbook(all_cell_type_table,file=file.path(out_basepath,"Tables/all_cell_type_table.xlsx"),overwrite=T)

