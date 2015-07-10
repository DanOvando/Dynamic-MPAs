static__maxinterest_plot_fun <- function(PlotData, Theme)
{
  
  
  static__maxinterest_plot <- (ggplot(data = PlotData  ,
                                    aes(x = ReserveSize,y = MaxInterestRate/100,fill = LoanType)) +
                               geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
                               facet_wrap(~Species, scales = 'free_y') + 
                               scale_fill_manual(values = c('green','red','steelblue2'),name = element_blank()) + 
                               xlab('Size of Reserve') + 
                               ylab('Maximum Interest Rate') + 
                               scale_x_continuous(labels = percent) + 
                               scale_y_continuous(labels = percent) + 
                               Theme )
return(static__maxinterest_plot)
}