priceinc_surface_plot_fun <- function(PlotData, Theme)
{
  
  
  priceinc_surface_plot <- (ggplot(PlotData,aes(Intercept,ReserveSize))+
                              geom_raster(aes(fill=pmin(100,PriceInc)/100),interpolate = T)+facet_wrap(~Species,scales='fixed')+
                              scale_fill_gradient2(name = '% Price Gain Needed', labels = percent,low = 'green', mid = 'yellow', high = 'red', midpoint = .25) + 
                              scale_y_continuous(labels = percent, name = '% in Reserve') + 
                              Theme + 
                              theme(axis.line = element_blank()))
  
  return(priceinc_surface_plot)
  
}