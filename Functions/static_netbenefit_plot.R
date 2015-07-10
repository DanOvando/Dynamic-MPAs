static_netbenefit_plot_fun<- function(PlotData, Theme)
{
  
  
  static_netbenefit_plot<- (ggplot(data = PlotData  ,
                                  aes(x = ReserveSize,y = TimeToNPB,fill = TimeToNPB)) +
                             geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
                             facet_wrap(~Species,scales = 'free_y') + 
                             geom_hline(aes(yintercept = 0), linetype = 'longdash') + 
                             scale_fill_gradient(low = 'green', high = 'red') +
                             #                  scale_fill_brewer(palette = 'Dark2') +
                             xlab('Size of Reserve') + 
                             ylab('Years to Net Benefit') + 
                             scale_x_continuous(labels = percent) + 
                             Theme + 
                             theme(legend.position = 'none'))
  
  return(static_netbenefit_plot)
}