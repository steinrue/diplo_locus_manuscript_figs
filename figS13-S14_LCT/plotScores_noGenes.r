#This sample script plot scan scores with tracks of filters and gene
##Load required packages
library(data.table)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(latex2exp)


##plot the stat scores of the designated region
### use default palette of blue --> red
plotManChrs <- function(dt, left, right, statname, colorstat,
                mask=FALSE, bed=NULL, color_midpoint=0, palname='Greys'){
  # colnames(dt) <- c('Chr','physPos','stat')
  p <- ggplot(dt)
  chr = unique(dt$Chr)[1]
  if(colorstat != FALSE){
    p <- p  + geom_point(aes(x=physPos, y=stat, color=colorscale), size=1.5) + 
              scale_x_continuous(name=paste('Positions on chromosome', ch, '(Mbp)'), labels=function(x){return(as.character(x/1e6))}, limits=c(left, right))
    ### add colorscale #ffffbf change midpoint from white-ish yellow to gray
    p <- p + scale_color_gradient2(name=colorstat, high="#940015", mid="#888888", low="#202584", midpoint=color_midpoint, na.value="#969696")
  }
  else{
    p <- p  + geom_point(aes(x=physPos, y=stat), color='#666666', size=1.5) + 
              scale_x_continuous(name=paste('Positions on chromosome', ch, '(Mbp)'), labels=function(x){return(as.character(x/1e6))}, limits=c(left, right))
  }
  ###adding mask (bedfile: chr start end feature)
  if(mask){
    if(nrow(bed)>0){
      yMin=min(dt$stat)
      width = (max(dt$stat)-yMin)/15
      ntype = length(unique(bed$type))
    if(ntype < 3){
      mycols = c("#a3a3a3","#939393")[1:ntype]
    }else{
      mycols = brewer.pal(ntype+1, palname)[2:(ntype+1)]
      names(mycols) = sort(unique(bed$type))
    }
    p <- p + geom_rect(data=bed, aes(xmin=start, xmax=end, ymin=yMin-width, ymax=yMin, fill=type) )+ 
      scale_fill_manual('Filters:', values = mycols)  + 
      theme(legend.position = 'top', legend.direction = 'horizontal', legend.title=element_text(size=15) )
    }
  }
  ###finishing touches
  ####statname can be customized (instead of a commandline argument) e.g. statname=expression(italic('B')[2])
  p <- p + ylab(statname) + #xlim(left, right) + 
    theme(panel.border=element_blank(), panel.grid.minor.y=element_blank(), panel.background=element_blank()
      )
  return(p)
}


##only plot the scores
plotSingleChrom <- function(ch, inputname, outputname, statname){
  DT = data.table(read.table(inputname, header=FALSE, sep="\t", comment="#"))
  colnames(DT) <- c('locus', 'ongrid_s2hat', 'ongrid_maxLogLikelihood', 's2hat', 'maxLogLikelihood', 'MLR', 'chi2_p')
  DT[, c("Chr", "physPos", "rsID", "placeholder4", "placeholder5")] <- tstrsplit(DT$locus, "_", fixed = TRUE)
  ### format
  DT$physPos <- as.numeric(DT$physPos)

  ###get peak[physPos > left & physPos < right]
  peak=DT

  ###plot the statistics
  if (statname == "MLR"){
    p1 <- plotManChrs(peak[,.(Chr, physPos, stat=MLR, colorscale=s2hat)], 
                      left=min(DT$physPos), right=max(DT$physPos), statname=TeX(r"($log_{10}LR_{max}$)"),
                      colorstat=FALSE) #TeX(r"($\hat{s}_{AA}$)")
    # label the top SNP
    top <- peak[which.max(MLR)]
    p1 <- p1 + geom_label_repel(data=top, aes(x=physPos, y=MLR, label=rsID), arrow=arrow(length=unit(4, "pt")),
                                  force=20, direction="both", min.segment.length=unit(8, "pt"), segment.linetype=1, #
                                  fill="white", size=6, box.padding=0.5)
    # annotate corresponding p value?
    neglogPval_labs = c(0, 3, 5, 8, 10, 12, 15)
    p1 <- p1 + scale_y_continuous(name=TeX(r"($log_{10}LR_{max}$)"),
                                  sec.axis=sec_axis(trans=~., name=TeX(r"($-log_{10}p_{\chi^2(1)}$)"),
                                                    breaks=neglogPval_to_LRT(neglogPval_labs),
                                                    labels=as.character(neglogPval_labs))
                                  )
    # need to rotate y axis label
    p1 <- p1 + theme_minimal() + theme(axis.title.y.left=element_text(angle=90),
                                       axis.title.y.right=element_text(angle=270, color='gray')
                                       )
  }else if (statname == "pval"){
    p1 <- plotManChrs(peak[,.(Chr, physPos, stat=-log10(chi2_p), colorscale=s2hat)],
                      left=min(DT$physPos), right=max(DT$physPos), statname=TeX(r"($-log_{10}p_{\chi^2(1)}$)"),
                      colorstat=FALSE ) #TeX(r"($\hat{s}_{AA}$)")
    # need to rotate y axis label
    p1 <- p1 + theme_minimal() + theme(axis.title.y=element_text(angle=90))
    # maybe label the top SNPs?
    tops <- peak[order(MLR, decreasing=TRUE)]
    tops$logP <- -log10(tops$chi2_p)
    tops <- tops[logP > 12]
    print(tops)
    p1 <- p1 + geom_label_repel(data=tops, aes(x=physPos, y=logP, label=rsID), arrow=arrow(length=unit(4, "pt")),
                                  force=20, direction="both", min.segment.length=unit(8, "pt"), segment.linetype=1, #
                                  fill="white", size=5, box.padding=0.5, ylim=c(9,NA))
  }

  ## draw y=0
  p1 <- p1 + geom_hline(yintercept=0, color='black', linewidth=0.5)

  ## add bonferroni line
  p1 <- p1 + geom_hline(yintercept=-log10(0.05/nrow(DT)), color='black', linewidth=0.5, linetype='dashed')
  
  ###final touches, do not remove the x-axis
  p1 <- p1 + xlab(paste('Positions on chromosome', ch) ) +
    theme(panel.grid.major.x=element_line(linewidth=0.3, color='#cccccc'), 
      legend.position = 'top', legend.direction = 'horizontal', 
      legend.key.width = unit(50, "pt"),
      legend.text=element_text(size=16), legend.title=element_text(size=16), 
      axis.title.y=element_text(size=18, hjust=0.5, vjust=0.5, margin=margin(r=0)), #, angle=0
      axis.text.y=element_text(size=15), axis.line=element_line(color='black', linewidth=0.5),   #
      axis.title.x=element_text(size=18, hjust=0.5), axis.text.x=element_text(size=15, color='black'),
      plot.margin=margin(b=5, unit='pt')
      )

  ggsave(outputname, plot=p1, width=14, height=7, units='in', dpi=500) 
  
  if(interactive())  return(DT)
}


if(!interactive()) {
    args <-commandArgs(trailingOnly=TRUE)
    print(args)
    inputname=args[1]
    ch=as.integer(args[2])
    statname=args[3]
    outputname=args[4]
    plotSingleChrom(ch, inputname, outputname, statname=statname)
}
