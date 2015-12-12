static_priceinc_plot_fun <- function(PlotData, Theme)
{
  
  PlotData$PriceInc <- pmin(50,PlotData$PriceInc)
  
  Breaks <- seq(0,50,by = 10)
  
  Labels <- c(paste(seq(0,40,by = 10),'%',sep = ''), '>50%')
  
  static__priceinc_plot<- (ggplot(data = PlotData  ,
                                                                    aes(x = ReserveSize,y = PriceInc),fill = 'grey76') +
#                                   aes(x = ReserveSize,y = PriceInc,fill = PriceInc)) +
                             
                             geom_bar(stat = 'identity', position = 'dodge', color = 'black', alpha = 0.85) +
                             facet_wrap(~Species, scales = 'fixed') + 
                             geom_hline(aes(yintercept = 0), linetype = 'longdash') + 
#                              scale_fill_gradient(guide = F,low = 'green', high = 'red') +
                             #                  scale_fill_brewer(palette = 'Dark2') +
                             xlab('Size of Reserve') + 
                             ylab('Price Increase Needed') + 
                             scale_x_continuous(labels = percent) + 
                             scale_y_continuous(labels = Labels, breaks = Breaks) + 
                             Theme + 
                             theme(legend.position = 'none',strip.text.x = element_text(size = 7)))
  return(static__priceinc_plot) 
}