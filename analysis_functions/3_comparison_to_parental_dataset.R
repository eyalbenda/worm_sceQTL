source("analysis_functions/2_comparison_to_other_datasets.R")

dfCisAndParents = data.frame()
corDF = data.frame()
for(tis in names(tissueCis))
{
  parent_diff_exp = readRDS((paste("study_data/parent_differential_expression//difExp_",tis,".rds.gz",sep="")))
  cur_beta_norm = all_betas_nofilter[match(rownames(parent_diff_exp),rownames(all_betas)),tis]/all_SE_nofilter[match(rownames(parent_diff_exp),rownames(all_betas)),tis]
  cur_beta_norm = cur_beta_norm[is.finite(cur_beta_norm)]
  if(length(cur_beta_norm)<2)
    next;
  parent_diff_exp = parent_diff_exp[is.finite(log2(parent_diff_exp$norm_foldChange)),]
  parent_diff_exp = parent_diff_exp[match(intersect(names(cur_beta_norm),rownames(parent_diff_exp)),rownames(parent_diff_exp)),]
  colnames(tissueCis[[tis]])[3] = "theta.dispersion"
  # plot(cur_beta_norm[match(rownames(parent_diff_exp),names(cur_beta_norm))],log2(parent_diff_exp$norm_foldChange))  
  curX = tissueCis[[tis]][match(rownames(parent_diff_exp),tissueCis[[tis]]$transcript),]
  dfCisAndParents = rbind(dfCisAndParents,data.frame(tissue = tis,curX,
                                                     parent_diff_exp))
  signif=curX$FDR<0.1
  if(sum(signif)==0)
    next;
  corDF = rbind(corDF,
                data.frame(tissue=tis,
                           makeCorText(log2(parent_diff_exp$norm_foldChange)[signif],curX$negbin.beta[signif]/curX$negbin.se[signif],mult = 0.8)))
  rm(curX)
}

isDifExpressed = (dfCisAndParents$pvalue.adj.FDR<0.1)&abs(log(base=2,dfCisAndParents$norm_foldChange))>1
fisher.test(cbind(c(sum(isDifExpressed&dfCisAndParents$FDR<0.1),
                    sum(isDifExpressed&dfCisAndParents$FDR>0.1)),
                  c(sum(!isDifExpressed&dfCisAndParents$FDR<0.1),
                    sum(!isDifExpressed&dfCisAndParents$FDR>0.1))))
with(dfCisAndParents %>% filter(pvalue.adj.FDR<0.1&FDR<0.1&abs(log(base=2,norm_foldChange))>1),sum(sign(log2(norm_foldChange))!=sign(negbin.beta)))
