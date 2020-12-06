source("analysis_functions//1_load_data_and_functions.R")
source("figures_and_tables//table_formatting.R")
#### Generate all tissue eQTL table
all_eqtl_table = NULL
for(tis in names(tissuePeaks))
{
  curDF = format_eqtl_table(tissuePeaksObject = tissuePeaks[[tis]],tissueCisObject = tissueCis[[tis]])
  all_eqtl_table = rbind(all_eqtl_table,data.frame(tissue=tis,curDF,stringsAsFactors = F))
}
all_eqtl_table_wb = make_excel_workbook(list("eQTL in cell types"=all_eqtl_table))
#openxlsx::saveWorkbook(all_eqtl_table_wb,file=file.path(out_basepath,"Tables/all_eqtl_table.xlsx"),overwrite=T)
