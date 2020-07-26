#' Plot signature contributions with 95 % CI
#'
#' @param mutSigsSummary data frame of bootstrapping output in the format from mm_fit_signatures
#'
#' @return plots showing mutational signature estimates with 95 % CI by sample
#' @export
#'
bootSigsPlot <- function(mutSigsSummary){

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
                           values = setNames(c(RColorBrewer::brewer.pal(8, "Dark2"), "lightblue", "#FB8072"),
                                             c("SBS1", "SBS2", "SBS5", "SBS8", "SBS9", "SBS13", "SBS18", "SBS-MM1", "SBS35", "SBS84")),
                           ...)
  }

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
