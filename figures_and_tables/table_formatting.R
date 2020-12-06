require(openxlsx)
if(!exists("generatePlotTissue"))
  warning("WARNING: functions will not work without loading function generatePlotTissue from 1_load_data_and_functions.R")
format_eqtl_table = function(tissuePeaksObject,tissueCisObject)
{
generatePlotTissue(tissuePeaksObject,tissueCisObject,dfReturn=T)%>% mutate(normEffect = negbin.beta/negbin.se) %>%
  dplyr::select(c(`Type` = "type",`QTL chromosome` = "chr",     
                  `QTL position` = "pos",`CI_left` = "ciL",
                  `CI_right` = "ciR",`Transcript` = "transcript",`Normalized effect size` = "normEffect",
                  "FDR" = "FDR","Transcript chromosome" = "tchr",
                  `Start` = "start_position",`End` = "end_position",
                  `Strand` = "strand", `external_gene_id`,
                  `external_transcript_id`,`wbps_transcript_id`))
}

### table list should have name = data.frame, name will e used for sheet name and header (row 1).
make_excel_workbook = function(table_list,file=NULL)
{
  wb <- createWorkbook()
  for(table in 1:length(table_list))
  {
    name = names(table_list)[table]
    df = table_list[[table]]
    # add Excel sheet
    addWorksheet(wb, name)
    # make first row
    addStyle(wb,sheet = name, style = createStyle(textDecoration = "Bold", halign = "left"),rows = 1,cols=1)
    header_st <- createStyle(textDecoration = "Bold",halign = "center")
    # Write data with header style defined above
    addStyle(wb,sheet = name, style = createStyle(halign = "center"),rows = 2:(nrow(df)+1),cols=1:ncol(df),gridExpand=T)
    writeData(wb, name, df,startRow = 2, headerStyle = header_st)
    setColWidths(wb,name,1:ncol(df),widths = "auto")
    setColWidths(wb,name,1,widths = (max(stringr::str_length(df[,1]))+2))
    writeData(wb,name,sprintf("Table ##. %s",name))
  }
  # save to .xlsx file
  if(is.null(file))
  {
    return(wb)
  } else
  {
    saveWorkbook(wb, file,overwrite=T)
  }
}
