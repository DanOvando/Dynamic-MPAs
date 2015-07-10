npb_surface_plot_fun <- function(PlotData, Theme)
{
  npb_surface_plot <- (ggplot(PlotData,aes(Intercept,FinalReserve))+
                      geom_raster(aes(fill=NPB),interpolate = T)+facet_wrap(~Species,scales='fixed')+
                      scale_fill_gradient2(name = 'NPB',low = 'red', mid = 'yellow', high = 'green', midpoint = 0) + 
                      scale_y_continuous(labels = percent, name = '% in Reserve') + 
                      Theme + 
                      theme(axis.line = element_blank()))
  return(npb_surface_plot)
}