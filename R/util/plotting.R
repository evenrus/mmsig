## Plotting functions for mmsig

rotatedAxisElementText = function(angle,position='x'){
    angle     = angle[1]; 
    position  = position[1]
    positions = list(x=0,y=90,top=180,right=270)
    if(!position %in% names(positions))
        stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
    if(!is.numeric(angle))
        stop("'angle' must be numeric",call.=FALSE)
    rads  = (angle - positions[[ position ]])*pi/180
    hjust = 0.5*(1 - sin(rads))
    vjust = 0.5*(1 + cos(rads))
    element_text(angle=angle,vjust=vjust,hjust=hjust)
}

scale_fill_sigs <- function(...){
    ggplot2:::manual_scale('fill', 
                           values = setNames(c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(1, "Set3")),
                                             c("SBS1", "SBS2", "SBS5", "SBS8", "SBS9", "SBS13", "SBS18", "SBS-MM1", "SBS84")), 
                           ...)
}

plot_signatures = function(sig, filename = NULL){
  "
  Plotting function for mmSig package
  filename: path to save plot. If NULL (default), returns plot to environment
  "
  if(is.null(filename)){
      par(mfrow=c(1,1), mar=c(7,5,5,10), xpd=T)
      barplot(as.matrix(t(sig[,-ncol(sig)])), 
              col = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")), 
              las=2, cex.names = 0.7, names.arg = sig$sampleID, space=rep(0, nrow(sig)), border = NA)
      legend("topright",legend=(colnames(sig)[-ncol(sig)]),bty="n", pch=15, 
            col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
            cex=1, pt.cex=1, inset=c(-0.15,0.0),x.intersp = 1,y.intersp = 1)
  } else{
      pdf(filename, width = 12, height = 7)
      barplot(as.matrix(t(sig[,-ncol(sig)])), 
              col = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")), 
              las=2, cex.names = 0.7, names.arg = sig$sampleID, space=rep(0, nrow(sig)), border = NA)
      legend("topright",legend=(colnames(sig)[-ncol(sig)]),bty="n", pch=15, 
            col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
            cex=1, pt.cex=1, inset=c(-0.15,0.0),x.intersp = 1,y.intersp = 1)
      dev.off()
  }
}

bootSigsPlot <- function(mutSigsSummary){
    ggplot(mutSigsSummary, aes(signature, estimate, fill = signature))+
        geom_bar(position="dodge", stat="identity")+
        geom_errorbar(data = mutSigsSummary, mapping = aes(x = signature, ymin = CI025, ymax = CI975))+
        facet_grid(~sample)+
        theme(text = element_text(size = 12),
              axis.text.x = rotatedAxisElementText(90, 'top'),
              axis.title.x = element_blank(),
              strip.background = element_blank(),
              legend.position = 'none')+
        scale_fill_sigs()+
        labs(y = 'Relative contribution')
}
