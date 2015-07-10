
yield_surface_plot_fun<- function(PlotData, Theme)
{
  
  yield_surface_plot <- (ggplot(PlotData,aes(Intercept,FinalReserve))+
                        geom_raster(aes(fill=net_yields),interpolate = T)+facet_wrap(~Species,scales='fixed')+
                        scale_fill_gradient(name = 'Cumulative Yields',low = 'red', high = 'green') + 
                        scale_y_continuous(labels = percent, name = '% in Reserve') + 
                        Theme + 
                        theme(axis.line = element_blank()))
  return(yield_surface_plot)
}