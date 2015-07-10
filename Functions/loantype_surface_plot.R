loantype_surface_plot_fun<- function(PlotData, Theme)
{
  
  loantype_surface_plot <- (ggplot(PlotData,aes(Intercept,ReserveSize))+
                           geom_raster(aes(fill=LoanType))+facet_wrap(~Species,scales='fixed')+
                           scale_fill_manual(values = c('green','red','steelblue2'),name = element_blank()) + 
                           scale_y_continuous(labels = percent, name = '% in Reserve') + 
                           Theme + 
                           theme(axis.line = element_blank()))
  
  return(loantype_surface_plot)
}