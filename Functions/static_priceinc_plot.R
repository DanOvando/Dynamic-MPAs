static__priceinc_plot_fun <- function(PlotData, Theme)
{
  
  
  static__priceinc_plot<- (ggplot(data = PlotData  ,
                                 aes(x = ReserveSize,y = PriceInc/100,fill = PriceInc)) +
                            geom_bar(stat = 'identity', position = 'dodge', color = 'black', alpha = 0.85) +
                            facet_wrap(~Species, scales = 'fixed') + 
                            geom_hline(aes(yintercept = 0), linetype = 'longdash') + 
                            scale_fill_gradient(low = 'green', high = 'red') +
                            #                  scale_fill_brewer(palette = 'Dark2') +
                            xlab('Size of Reserve') + 
                            ylab('Price Increase Needed') + 
                            scale_x_continuous(labels = percent) + 
                            scale_y_continuous(labels = percent) + 
                            Theme + 
                            theme(legend.position = 'none'))
 return(static__priceinc_plot) 
}