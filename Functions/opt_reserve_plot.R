opt_reserve_plot_fun<- function(OptRun, Theme) 
{
  
  opt_reserve_plot<- (ggplot(data = OptRun ,aes(x = Year, y = CurrentReserve))+geom_line(size = 1.2)+
                      geom_point(aes(fill = NPB), shape = 21) + 
                      facet_wrap(~Species,scales  = 'fixed') +
                      scale_fill_gradient2(name = 'NPB',low = 'red', mid = 'yellow', high = 'green',midpoint = 1) + 
                      SimpleTheme + 
                      geom_hline(aes(yintercept = 0, linetype = 'longdash')) + 
                      xlab('Year') + 
                      scale_y_continuous(labels = percent) + 
                      theme(legend.position = 'right'))
  
  return(opt_reserve_plot)
}
