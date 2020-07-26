#' Plotting function for mmsig package
#'
#' @param sig data frame of mutational signature estimates in the format from mm_fit_signatures
#' @param sig_order order of signatures to plot in stacked bar chart, from top to bottom
#' @param samples boolean whether or not to plot sample names on the x axis
#'
#' @return relative contribution of mutational signatures in each sample
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 labs
#' @export
#'
plot_signatures = function(sig,
                           sig_order = c("SBS1", "SBS2", "SBS13", "SBS5", "SBS8", "SBS9", "SBS18", "SBS-MM1", "SBS35"),
                           samples = FALSE){

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

  sigPlot <- sig %>%
    rownames_to_column(var = "sample") %>%
    dplyr::select(-mutations) %>%
    melt(id.var = "sample", variable.name = "SBS", value.name = "prop") %>%
    mutate(SBS = factor(SBS, levels = sig_order)) %>%
    ggplot(aes(sample, prop, fill = SBS)) +
    geom_col(width=1)+
    scale_fill_sigs()+
    scale_y_continuous(expand = c(0,0))+
    theme_bw()+
    labs(x = "Sample",
         y = "Relative contribution",
         fill = "Signature")+
    theme(text = element_text(size = 10, color = "black"),
          axis.text = element_text(color = "black"),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.title.x = element_blank())

  if(samples){
    sigPlot +
      theme(axis.text.x = rotatedAxisElementText(90, "top"))
  } else{
    sigPlot +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
}
