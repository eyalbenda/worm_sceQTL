## Define outoput directory for plots and tables
out_basepath = ""


if(!"devtools" %in% installed.packages())
  install.packages("devtools")
if(!"ggpattern" %in% installed.packages())
  remotes::install_github("coolbutuseless/ggpattern")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!"monocle3" %in% installed.packages())
{
  BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                         'limma', 'S4Vectors', 'SingleCellExperiment',
                         'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
  devtools::install_github('cole-trapnell-lab/leidenbase')
  devtools::install_github('cole-trapnell-lab/monocle3')
}
cran.packages =  c("tidyverse","ggrepel","DescTools","openxlsx","cowplot","foreach","parallel","doParallel","pbapply")
bioc.packages =  c("ComplexHeatmap","DEGreport","topGO","biomaRt","GEOquery","Biostrings")
new.packages <- cran.packages[!(cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

new.packages.bioc <- bioc.packages[!(bioc.packages %in% installed.packages()[,"Package"])]
if(length(new.packages.bioc)) BiocManager::install(new.packages.bioc)

require(ComplexHeatmap)
require(ggrepel)
require(ggpattern)
require(monocle3)
require(DescTools)
require(DEGreport)
require(tidyverse)
require(openxlsx)
require(cowplot)
library(topGO)
library(biomaRt)
require(foreach)
require(parallel)
require(doParallel)
require(pbapply)
require(GEOquery)
require(Biostrings)

##### Load files #####
#### Load gene annotation
geneAnnotsCur = readRDS("study_data//geneAnnotationMatrix.rds.gz")
#### Load global monocle object
cdsFinBC = readRDS("study_data/backgroundCorrectedMonocle3_V4.rds")


#### Load eQTL results
jointPeaks = readRDS("study_data/Feb21//jointModel_globalPeaks.RDS")
tissueCis = readRDS("study_data/Feb21//perTissue_cisOnly.RDS")
tissueHotspot = readRDS("study_data/Apr9/perTissue_hotspots.RDS")
tissuePeaks = readRDS("study_data/Feb21//perTissue_globalPeaks.RDS")
jointCis = readRDS("study_data/Feb21//jointModel_cisOnly.RDS")
isExpressed = readRDS("study_data/Feb21//expressed_in_20_cells_broad_tissue.RDS")
covariates = readRDS("study_data/Feb21//covariates.RDS")
assertthat::are_equal(covariates$barcode,colData(cdsFinBC)$barcode)
fineTissueCis = readRDS("study_data/fineTissue_20_cutoff.rds")

neuCDS = readRDS("study_data/neuronalCDS_V2.rds")

##### Functions #####
#### plot eQTL
plotQTLdf=function(plot_df)ggplot(plot_df,aes(x=pos,y=start_position)) + xlab("QTL position (mb)") + ylab("transcript position (mb)") + scale_x_continuous(breaks = c(5e6,1.5e7),labels = function(x)x/1e6) + scale_y_continuous(breaks = c(5e6,1.5e7),labels = function(x)x/1e6) + facet_grid(chromosome_name~chr,scales="free") + geom_point(size=0.4) + theme_classic() + geom_errorbarh(aes(xmin=ciL,xmax=ciR,y=start_position,alpha=-log10(FDR)),alpha=0.4)

#### Generate eQTL plot for the joint eQTL calls
generatePlot = function(jointPeaks,cisPeaks = NULL,dfReturn=F,cutoff=0.1)
{
  colnames(jointPeaks)[colnames(jointPeaks)=="peak.marker"] = "marker"
  colnames(cisPeaks)[colnames(cisPeaks)=="lmarker"] = "cisMarker"
  colnames(cisPeaks)[colnames(cisPeaks)=="LLik"] = "LOD"
  jointPeaks = data.frame(jointPeaks,chr=gsub("_+.*","",jointPeaks$marker),pos=as.numeric(gsub(".*_","",jointPeaks$marker)),
                          ciL = as.numeric(gsub(".*_","",jointPeaks$CI.l)),ciR = as.numeric(gsub(".*_","",jointPeaks$CI.r)),type="trans")
  jointPeaksGeneAnnots = geneAnnotsCur[match(jointPeaks$transcript,geneAnnotsCur$wbps_gene_id),] 
  jointPeaks = jointPeaks[jointPeaksGeneAnnots$chromosome_name!="MtDNA",]
  jointPeaksGeneAnnots = jointPeaksGeneAnnots[jointPeaksGeneAnnots$chromosome_name!="MtDNA",]
  
  plot_df = data.frame(jointPeaks,jointPeaksGeneAnnots) %>% filter(FDR<cutoff)
  plot_df$chromosome_name = factor(plot_df$chromosome_name,levels=c("X","V","IV","III","II","I"))
  if(!is.null(cisPeaks))
  {
    cisPeaks = cisPeaks %>% filter(FDR<cutoff)
    cisPeaksGeneAnnots = jointPeaksGeneAnnots = geneAnnotsCur[match(cisPeaks$transcript,geneAnnotsCur$wbps_gene_id),]
    cisPeaks = data.frame(ntissues=cisPeaks$ntissues,chrom=cisPeaksGeneAnnots$chromosome_name,transcript=cisPeaks$transcript,
                          marker=cisPeaks$cisMarker,CI.l=cisPeaksGeneAnnots$start_position,CI.r=cisPeaksGeneAnnots$start_position,LOD=cisPeaks$LOD,FDR=cisPeaks$FDR,
                          chr=cisPeaksGeneAnnots$chromosome_name,pos=cisPeaksGeneAnnots$start_position,tchr=cisPeaksGeneAnnots$chromosome_name,tpos=cisPeaksGeneAnnots$start_position,chr.1=cisPeaksGeneAnnots$chromosome_name,pos.1=cisPeaksGeneAnnots$start_position,
                          ciL=cisPeaksGeneAnnots$start_position,ciR=cisPeaksGeneAnnots$start_position,type="cis")
    toRemove = (plot_df$chrom==plot_df$tchr) & sapply(1:nrow(plot_df),function(x)(c(plot_df$CI.l[x]-1e6,plot_df$CI.r[x]+1e6) %overlaps% c(plot_df$start_position[x],plot_df$end_position[x])))
    plot_df = rbind(plot_df[!toRemove,],data.frame(cisPeaks,cisPeaksGeneAnnots))
  }
  if(dfReturn)
  {
    plot_df 
  }else
  {
    plotQTLdf(plot_df)
  }
}


### Generate eQTL plot for tissue data
generatePlotTissue = function(jointPeaks,cisPeaks = NULL,dfReturn=F,cutoff=0.1)
{
  colnames(jointPeaks)[colnames(jointPeaks)=="peak.marker"] = "marker"
  colnames(cisPeaks)[colnames(cisPeaks)=="lmarker"] = "cisMarker"
  colnames(cisPeaks)[colnames(cisPeaks)=="LLik"] = "LOD"
  jointPeaks = data.frame(jointPeaks,chr=gsub("_+.*","",jointPeaks$marker),pos=as.numeric(gsub(".*_","",jointPeaks$marker)),
                          ciL = as.numeric(gsub(".*_","",jointPeaks$CI.l)),ciR = as.numeric(gsub(".*_","",jointPeaks$CI.r)),type="trans")
  jointPeaksGeneAnnots = geneAnnotsCur[match(jointPeaks$transcript,geneAnnotsCur$wbps_gene_id),] 
  jointPeaks = jointPeaks[jointPeaksGeneAnnots$chromosome_name!="MtDNA",]
  jointPeaksGeneAnnots = jointPeaksGeneAnnots[jointPeaksGeneAnnots$chromosome_name!="MtDNA",]
  
  plot_df = data.frame(jointPeaks,jointPeaksGeneAnnots) %>% filter(FDR<cutoff)
  plot_df$chromosome_name = factor(plot_df$chromosome_name,levels=c("X","V","IV","III","II","I"))
  if(!is.null(cisPeaks))
  {
    cisPeaks = cisPeaks %>% filter(FDR<cutoff)
    cisPeaksGeneAnnots = jointPeaksGeneAnnots = geneAnnotsCur[match(cisPeaks$transcript,geneAnnotsCur$wbps_gene_id),]
    cisPeaks = data.frame(chrom=cisPeaksGeneAnnots$chromosome_name,transcript=cisPeaks$transcript,
                          marker=cisPeaks$cisMarker,CI.l=cisPeaksGeneAnnots$start_position,CI.r=cisPeaksGeneAnnots$start_position,LOD=cisPeaks$LOD,FDR=cisPeaks$FDR,
                          chr=cisPeaksGeneAnnots$chromosome_name,pos=cisPeaksGeneAnnots$start_position,tchr=cisPeaksGeneAnnots$chromosome_name,tpos=cisPeaksGeneAnnots$start_position,
                          negbin.beta=cisPeaks$negbin.beta,negbin.se=cisPeaks$negbin.se,negbin.p=cisPeaks$negbin.p,chr.1=cisPeaksGeneAnnots$chromosome_name,pos.1=cisPeaksGeneAnnots$start_position,ciL=cisPeaksGeneAnnots$start_position,ciR=cisPeaksGeneAnnots$start_position,type="cis")
    toRemove = (plot_df$chrom==plot_df$tchr) & sapply(1:nrow(plot_df),function(x)(c(plot_df$CI.l[x]-1e6,plot_df$CI.r[x]+1e6) %overlaps% c(plot_df$start_position[x],plot_df$end_position[x])))
    plot_df = rbind(plot_df[!toRemove,],data.frame(cisPeaks,cisPeaksGeneAnnots))
  }
  if(dfReturn)
  {
    plot_df 
  }else
  {
    plotQTLdf(plot_df)
  }
}

#### Plot genotypes in a chromosome 
plotGenos = function(chr,barcode)
{
  if(!exists(genotypes)|!exists(refCounts)|!exists(altCounts))
    stop("genotypes,refCounts and altCounts must be loaded")
  if(chr!="all")
  {
    plotDF = data.frame(hmm=genotypes[covariates$barcode==barcode,grep(sprintf("^%s_",chr),colnames(genotypes))],
                        alt=altCounts[grep(sprintf("^%s_",chr),rownames(altCounts)),barcode],ref=refCounts[grep(sprintf("^%s_",chr),rownames(refCounts)),barcode])
    plotDF = data.frame(chr = chr,position = as.numeric(gsub("+.*_","",rownames(plotDF))),plotDF)
    pl = ggplot(plotDF,aes(x=position,y=hmm)) + geom_line() + ylim(c(-0.03,1.03)) + theme_classic() + geom_hline(yintercept=0.5,linetype=3) + geom_linerange(data=plotDF %>% filter(ref>0),aes(x=position,ymax=0.5,ymin=0)) + geom_linerange(data=plotDF %>% filter(alt>0),aes(x=position,ymax=1,ymin=0.5)) + 
      geom_text(data=plotDF %>% filter(alt>0),aes(x=position,y=1.03,label=alt))+ geom_text(data=plotDF %>% filter(ref>0),aes(x=position,y=-0.03,label=ref)) + scale_x_continuous(breaks = c(5e6,1e7,1.5e7,2e7),labels = function(x)x/1e6)
  } else
  {
    plotDFAll = NULL
    for(chr in c("I","II","III","IV","V","X"))
    {
      plotDF = data.frame(hmm=genotypes[covariates$barcode==barcode,grep(sprintf("^%s_",chr),colnames(genotypes))],
                          alt=altCounts[grep(sprintf("^%s_",chr),rownames(altCounts)),barcode],ref=refCounts[grep(sprintf("^%s_",chr),rownames(refCounts)),barcode])
      plotDF = data.frame(chr = chr,position = as.numeric(gsub("+.*_","",rownames(plotDF))),plotDF)
      plotDFAll = rbind(plotDFAll,plotDF)
    }
    pl = ggplot(plotDFAll,aes(x=position,y=hmm)) + geom_line() + ylim(c(-0.03,1.03)) + theme_classic() + geom_hline(yintercept=0.5,linetype=3) + geom_linerange(data=plotDFAll %>% filter(ref>0),aes(x=position,ymax=0.5,ymin=0),alpha=0.5) + geom_linerange(data=plotDFAll %>% filter(alt>0),aes(x=position,ymax=1,ymin=0.5),alpha=0.5) + 
      geom_text(data=plotDFAll %>% filter(alt>0),aes(x=position,y=1.03,label=alt))+ scale_x_continuous(breaks = c(5e6,1e7,1.5e7,2e7),labels = function(x)x/1e6) + facet_grid(~chr,scales="free_x")+geom_text(data=plotDFAll %>% filter(ref>0),aes(x=position,y=-0.03,label=ref))
  }
  pl
}
makeCorText = function(x,y,top=20,mult=0.85)
{
  corPval = format.pval(signif(cor.test(x,y,method="spearman")$p.value),2)
  if(!grepl("<",corPval)) 
  {
    corPval = paste("p == ",corPval)
  } else
  {
    corPval = paste("p ",corPval)
  }
  rbind(data.frame(x=min(x[is.finite(x)]),y=top,label=paste("rho == ",signif(cor(x,y,method="spearman"),2))),
        data.frame(x=min(x[is.finite(x)]),y=top*mult,label=corPval))
}

makeFSText = function(fs,top=20,mult=0.85)
{
  fsPval = format.pval(signif(fs$p.value,2))
  if(!grepl("<",fsPval)) 
  {
    fsPval = paste("p == ",fsPval)
  } else
  {
    fsPval = paste("p ",fsPval)
  }
  rbind(data.frame(x=0.1,y=top,label=paste("OR == ",signif(fs$estimate,2))),
        data.frame(x=0.1,y=top*mult,label=fsPval))
}

