#### This script generates the single cell eQTL matrix. All analysis should go here
#### All plotting should go in dedicated functions

source("analysis_functions/1_load_data_and_functions.R")

######## Extract effect sizes for all cis- eQTL ######
### get all betas
## generate empty matrices and copy dimnames
celltypeCis = tissueCis
all_betas = all_betas_nofilter = matrix(nrow=length(unique(unlist(sapply(celltypeCis,rownames)))),ncol=length(celltypeCis))
all_SE = all_SE_nofilter = matrix(nrow=length(unique(unlist(sapply(celltypeCis,rownames)))),ncol=length(celltypeCis))
colnames(all_betas) = colnames(all_SE) = colnames(all_betas_nofilter) = colnames(all_SE_nofilter) = names(celltypeCis)
rownames(all_betas) = rownames(all_SE) = rownames(all_betas_nofilter) = rownames(all_SE_nofilter) = unique(unlist(sapply(celltypeCis,rownames)))
thresh = 0.1
#### Get all genes that are significant at 10% fdr at least in one celltype
significant_at_least_in_one = unique(unlist(sapply(celltypeCis,function(x)rownames(x)[x$FDR<thresh])))

## populate beta and SE tables
for(tis in names(celltypeCis))
{
  all_betas[rownames(celltypeCis[[tis]])[celltypeCis[[tis]]$FDR<thresh],tis] = celltypeCis[[tis]]$negbin.beta[celltypeCis[[tis]]$FDR<thresh]
  all_SE[rownames(celltypeCis[[tis]])[celltypeCis[[tis]]$FDR<thresh],tis] = celltypeCis[[tis]]$negbin.se[celltypeCis[[tis]]$FDR<thresh]
  all_betas_nofilter[rownames(celltypeCis[[tis]]),tis] = celltypeCis[[tis]]$negbin.beta
  all_SE_nofilter[rownames(celltypeCis[[tis]]),tis] = celltypeCis[[tis]]$negbin.se
}


### Normalize betas by standard error
assertthat::are_equal(is.na(all_betas),is.na(all_SE))
all_normalized_betas = all_betas/all_SE

### combined normalized betas using Stouffer's method (weighing by sample size)
# number of cells per celltype (ordered based on normalized beta matrix)
n_cells_per_celltype = table(colData(cdsFinBC)$broad_tissue)[colnames(all_normalized_betas)]

# make empty matrix
combined_normalized_beta = rep(NA,nrow(all_normalized_betas))
names(combined_normalized_beta) = rownames(all_normalized_betas)
# populate by normalizing all non missing betas across all cell-types
for(i in 1:nrow(all_normalized_betas))
{
  curSubset = !is.na(all_normalized_betas[i,])
  combined_normalized_beta[i] = sum(n_cells_per_celltype[curSubset]*all_normalized_betas[i,curSubset])/sqrt(sum(n_cells_per_celltype[curSubset]^2))
}
