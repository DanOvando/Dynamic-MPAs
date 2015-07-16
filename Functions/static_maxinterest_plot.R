static__maxinterest_plot_fun <- function(PlotData, Theme)
{
  
  PlotData$MaxInterestRate[PlotData$MaxInterestRate >= 25] <- 25
  
  Breaks <- seq(0,25,by=5)
  
  Labels <- c(paste(seq(0,20,by = 5),'%',sep = ''), '>25%')
  
  static__maxinterest_plot <- (ggplot(data = PlotData  ,
                                    aes(x = ReserveSize,y = MaxInterestRate/100,fill = LoanType)) +
                               geom_bar(stat = 'identity', position = 'dodge', color = 'black', alpha = 0.85) +
                               facet_wrap(~Species, scales = 'fixed') + 
                               scale_fill_manual(values = c('green','red','steelblue2'),name = element_blank()) + 
                               xlab('Size of Reserve') + 
                               ylab('Maximum Interest Rate') + 
                               scale_x_continuous(labels = percent) + 
                               scale_y_continuous(labels = Labels) + 
                               Theme )
return(static__maxinterest_plot)
}