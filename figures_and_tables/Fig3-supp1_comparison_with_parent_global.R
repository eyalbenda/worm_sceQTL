source("analysis_functions//3_comparison_to_parental_dataset.R")


##### Plot cis eQTL against global parent differential expression ####

parents_global_differential_expression = readRDS("study_data//L2parents_DEsingle_output.rds.gz")
cis_and_parents_global = data.frame(ourdata = combined_normalized_beta[match(rownames(parents_global_differential_expression),names(combined_normalized_beta))],
                                    parents_global_differential_expression)
cis_and_parents_global = cis_and_parents_global[!is.na(cis_and_parents_global$ourdata),]
cis_and_parents_global["significant"] = rownames(cis_and_parents_global) %in% significant_at_least_in_one
corDF = with(cis_and_parents_global %>% filter(significant),makeCorText(log2(norm_foldChange),ourdata))

ggplot(cis_and_parents_global,aes(x=log2(norm_foldChange),y=ourdata)) + geom_point(size=0.2,aes(color=significant)) + theme_classic()  + 
  scale_colour_manual(values=c(viridis::inferno(10)[2])) + geom_text(parse=T,data=corDF,aes(x=x,y=y,label=label),hjust=0) + ylab(expression("standardized "*beta*" from cis-eQTL model")) + xlab(expression(log[2]*"(fold change)")) + theme(legend.position = "none")
#ggsave(file.path(out_basepath,"Figures/Supplementary Figures/Previous Data/ParentAndI.png"),width=4,height=4,dpi=300)
