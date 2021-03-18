source("analysis_functions/1_load_data_and_functions.R")

plotGenos = function(chr,barcode)
{
  if(chr!="all")
  {
    plotDF = data.frame(hmm=genotypes[rownames(genotypes)==barcode,grep(sprintf("^%s_",chr),colnames(genotypes))],
                        alt=altCounts[grep(sprintf("^%s_",chr),rownames(altCounts)),barcode],ref=refCounts[grep(sprintf("^%s_",chr),rownames(refCounts)),barcode])
    plotDF = data.frame(chr = chr,position = as.numeric(gsub("+.*_","",rownames(plotDF))),plotDF)
    pl = ggplot(plotDF,aes(x=position,y=hmm)) + geom_line() + ylim(c(-0.03,1.03)) + theme_classic() + geom_hline(yintercept=0.5,linetype=3) + geom_linerange(data=plotDF %>% filter(ref>0),aes(x=position,ymax=0.5,ymin=0)) + geom_linerange(data=plotDF %>% filter(alt>0),aes(x=position,ymax=1,ymin=0.5)) + 
      geom_text(data=plotDF %>% filter(alt>0),aes(x=position,y=1.03,label=alt))+ geom_text(data=plotDF %>% filter(ref>0),aes(x=position,y=-0.03,label=ref)) + scale_x_continuous(breaks = c(5e6,1e7,1.5e7,2e7),labels = function(x)x/1e6)
  } else
  {
    plotDFAll = NULL
    for(chr in c("I","II","III","IV","V","X"))
    {
      plotDF = data.frame(hmm=genotypes[rownames(genotypes)==barcode,grep(sprintf("^%s_",chr),colnames(genotypes))],
                          alt=altCounts[grep(sprintf("^%s_",chr),rownames(altCounts)),barcode],ref=refCounts[grep(sprintf("^%s_",chr),rownames(refCounts)),barcode])
      plotDF = data.frame(chr = chr,position = as.numeric(gsub("+.*_","",rownames(plotDF))),plotDF)
      plotDFAll = rbind(plotDFAll,plotDF)
    }
    pl = ggplot(plotDFAll,aes(x=position,y=hmm)) + 
      geom_line() + 
      geom_segment(data=plotDFAll %>% filter(ref>0),aes(xend=position,x=position,y=hmm-0.01,color=as.factor(ref),yend=hmm-0.1)) +
      geom_segment(data=plotDFAll %>% filter(alt>0),aes(xend=position,x=position,y=hmm+0.01,color=as.factor(alt),yend=hmm+0.1)) + 
      ylim(c(-0.03,1.03)) + 
      theme_classic() + 
      geom_hline(yintercept=0.5,linetype=3) +  
      scale_x_continuous(breaks = c(5e6,1e7,1.5e7,2e7),labels = function(x)x/1e6) + 
      facet_grid(~chr,scales="free_x") + 
      ylab("HMM\n output probabilities") +
      xlab("Chromosomal Position (Mb)") +
      scale_color_viridis_d(end=0.9) +
      labs(color="SNV Counts")
  }
  pl
}
refCounts = readRDS("study_data/refCounts_backgroundRemovedstrict.rds")
altCounts = readRDS("study_data/altCounts_backgroundRemovedstrict.rds")
genotypes = readRDS("study_data/Feb21/geno_matrix_figureS2.RDS") ## Full genotype matrix can't fit on github :/ so this is prechosen to have the median cell at slot 11
refCounts = refCounts[match(colnames(genotypes),rownames(refCounts)),]
altCounts = altCounts[match(colnames(genotypes),rownames(altCounts)),]


plotGenos(chr="all",barcode=rownames(genotypes)[11])
#ggsave(file.path(out_basepath,"Figures/Supplementary Figures/FigS2_hmm.png"),dpi=600,width=8,height=4)
