loan_surface_plot_fun <- function(PlotData, Theme)
{
  
  loan_surface_plot <- (ggplot(PlotData,aes(Intercept,ReserveSize))+
                       geom_raster(aes(fill=pmin(100,MaxInterestRate)/100),interpolate = T)+facet_wrap(~Species,scales='fixed')+
                       scale_fill_gradient2(name = 'Max Loan Rate', labels = percent,low = 'red', mid = 'yellow', high = 'green', midpoint = .25) + 
                       scale_y_continuous(labels = percent, name = '% in Reserve') + 
                       Theme + 
                       theme(axis.line = element_blank()))
}