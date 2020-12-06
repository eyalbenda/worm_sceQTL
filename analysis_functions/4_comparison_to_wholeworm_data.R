source("analysis_functions/2_comparison_to_other_datasets.R")
######## Comparison to whole worm eQTL data from Rockman et al 2010 ####

### Read rockman mapping data
rockman_mapping = readRDS("external_data/rockman_mapping.rds")
rockman_effect_sizes = readRDS("external_data/rockman_effect_sizes.rds")

### convert to data frame, create a "gene" column and filter only those who are in our dataset (7001 genes total)
rockman_effect_sizes = as.data.frame(rockman_effect_sizes) %>% 
  rownames_to_column(var="gene") 
rockman_effect_sizes$gene[stringr::str_count(rockman_effect_sizes$gene,"\\.")>1] = sub("^([^\\.]*\\.[^\\.]*).*", "\\1", 
                                                                                       
                                                                                       rockman_effect_sizes$gene[stringr::str_count(rockman_effect_sizes$gene,"\\.")>1])
rockman_effect_sizes = rockman_effect_sizes[!stringr::str_ends(rockman_effect_sizes$gene,"[b-z]"),] 
rockman_effect_sizes$gene[stringr::str_ends(rockman_effect_sizes$gene,"[a-z]")] = substr(rockman_effect_sizes$gene[stringr::str_ends(rockman_effect_sizes$gene,"[a-z]")],1,nchar(rockman_effect_sizes$gene[stringr::str_ends(rockman_effect_sizes$gene,"[a-z]")])-1)



rockman_effect_sizes = rockman_effect_sizes %>%
  filter(gene %in% geneAnnotsCur$wormbase_gseq[match(rownames(all_betas),geneAnnotsCur$wbps_gene_id)])



rockman_effect_sizes = rockman_effect_sizes[match(unique(rockman_effect_sizes$gene),rockman_effect_sizes$gene),]



rownames(rockman_effect_sizes) = rockman_effect_sizes$gene
# remove gene column, no longer needed
rockman_effect_sizes = rockman_effect_sizes[,-1]
### Get variant positions from column names
variant_positions = do.call(rbind,strsplit(gsub("X","",colnames(rockman_effect_sizes)),"_"))
### remap numbers to roman numerals (ugly)
variant_positions[variant_positions[,1]=="",1] = "X";variant_positions[variant_positions[,1]=="1",1] = "I";variant_positions[variant_positions[,1]=="2",1] = "II";variant_positions[variant_positions[,1]=="3",1] = "III";variant_positions[variant_positions[,1]=="4",1] = "IV";variant_positions[variant_positions[,1]=="5",1] = "V"
### find the nearest marker to each gene - 
### for query: attach a sorted gene annotation matrix and make GRanges object with genes start end 
### for subject: use variant position we extracted
cis_marker = nearest(x=with(geneAnnotsCur[match(rownames(rockman_effect_sizes),geneAnnotsCur$wormbase_gseq),],
                            GRanges(seqnames = chromosome_name,ranges = IRanges(start=start_position,end=end_position))),
                     subject = GRanges(seqnames = variant_positions[,1],IRanges(as.numeric(variant_positions[,2]),as.numeric(variant_positions[,2]))))
### Finally get the effect size for the cis marker for each gene
effectsize_cis_rockman = rep(NA,nrow(rockman_effect_sizes));
for(i in 1:nrow(rockman_effect_sizes))
  effectsize_cis_rockman[i] = rockman_effect_sizes[i,cis_marker[i]]
### We've already subset the rockman data to include only genes found in our dataset, so a simple match will do it
combined_normalized_beta_with_rockman = combined_normalized_beta[geneAnnotsCur$wbps_gene_id[match(rownames(rockman_effect_sizes),geneAnnotsCur$wormbase_gseq)]]
assertthat::noNA(combined_normalized_beta_with_rockman)

#### Make comparison tables
## Note: thresh is set by 7_comparison_to_ther_datasets.R
get_lod_score = function(r,n){
  return(-n * (log(1- r^2) / (2 *log(10))))
}
LODcutoff = 4.2 ## Derived from 100 permutations of the data to be a 10% FDR threshold

dfCisAndRockmanCis = data.frame(previous = -effectsize_cis_rockman,
                                ourdata = combined_normalized_beta_with_rockman,
                                significant=names(combined_normalized_beta_with_rockman) %in% significant_at_least_in_one)
dfCisAndRockmanCis = dfCisAndRockmanCis %>% rownames_to_column("transcript") 
dfCisAndRockmanCis = dfCisAndRockmanCis %>% mutate(rockmansig = get_lod_score(previous,n=208)>LODcutoff)
dfCisAndRockmanCis["PercentExpressed"] = rowSums(exprs(cdsFinBC[dfCisAndRockmanCis$transcript,])>0)/ncol(cdsFinBC)
dfCisAndRockmanCis["SumExpressed"] = log(rowSums(normalized_counts(cdsFinBC[dfCisAndRockmanCis$transcript,])))
# count number of hits below FDR threshold across cell types for each gene in the joint analysis
all_genes_fdr = do.call(rbind,lapply(celltypeCis,function(x)data.frame(transcript=x$transcript,FDR=x$FDR,stringsAsFactors = F))) %>% filter(transcript %in% dfCisAndRockmanCis$transcript)
dfCisAndRockmanCis["SumSignificantInOut"] = sapply(dfCisAndRockmanCis$transcript,function(x)sum(all_genes_fdr$FDR[all_genes_fdr$transcript==x]<thresh))
# significant in more datasets : significant in rockman, significant in our (at least ones), and the weighted-Z combined effect size has the same direction
sig_in_rockman_and_our=dfCisAndRockmanCis$transcript[which(dfCisAndRockmanCis$rockmansig&
                                                             dfCisAndRockmanCis$significant&
                                                             (sign(dfCisAndRockmanCis$previous)==sign(dfCisAndRockmanCis$ourdata)))]

### For downstream analyses, set discordant eQTL to false to not count them as successful replication
dfCisAndRockmanCis$significant[dfCisAndRockmanCis$rockmansig&dfCisAndRockmanCis$significant&(sign(dfCisAndRockmanCis$previous)!=sign(dfCisAndRockmanCis$ourdata))] = F
dfCisAndRockmanCis$SumSignificantInOut[!dfCisAndRockmanCis$significant] = 0
concordant_index = match(sig_in_rockman_and_our,dfCisAndRockmanCis$transcript)

enrichmentInRockmanRep = cbind(c(sum(dfCisAndRockmanCis$significant&dfCisAndRockmanCis$rockmansig),
                                 sum(dfCisAndRockmanCis$significant&!dfCisAndRockmanCis$rockmansig)),
                               c(sum(!dfCisAndRockmanCis$significant&dfCisAndRockmanCis$rockmansig),
                                 sum(!dfCisAndRockmanCis$significant&!dfCisAndRockmanCis$rockmansig)))
fsTest = fisher.test(enrichmentInRockmanRep)

enrichmentInRockmanMoreThanOne = cbind(c(sum(dfCisAndRockmanCis$SumSignificantInOut>1&dfCisAndRockmanCis$rockmansig),
                                         sum(dfCisAndRockmanCis$SumSignificantInOut>1&!dfCisAndRockmanCis$rockmansig)),
                                       c(sum(dfCisAndRockmanCis$SumSignificantInOut==1&dfCisAndRockmanCis$rockmansig),
                                         sum(dfCisAndRockmanCis$SumSignificantInOut==1&!dfCisAndRockmanCis$rockmansig)))

