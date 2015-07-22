static_comp_plot_fun <- function(PlotData,Theme,DiscRate)
{
  Breaks <- seq(10,30,by=5)
  
  PlotData <- PlotData %>%
    gather('VerbDiscountRate','NPB',which(colnames(PlotData) %in% c('FinalNPB','FinalNB')))
  
  PlotData$DiscountRate <- 0
  
  PlotData$DiscountRate[PlotData$VerbDiscountRate == 'FinalNPB'] <- DiscRate
  
  PlotData$TimeToBen <- PlotData$TimeToNPB
  
  PlotData$TimeToBen[PlotData$DiscountRate == 0] <- PlotData$TimeToNB[PlotData$DiscountRate == 0]
  
  PlotData$TimeToBen[ PlotData$TimeToBen>30] <- 30
  
  static_comp_plot <- (ggplot(data = PlotData  ,
                            aes(x = ReserveSize,y = NPB,fill = TimeToBen, group = DiscountRate)) +
                       geom_bar(stat = 'identity', position = 'dodge', color = 'black', alpha = 0.85) +
                       facet_wrap(~Species,scales = 'fixed') + 
                       geom_hline(aes(yintercept = 0), linetype = 'longdash') + 
                       scale_fill_gradient2(labels = c(seq(10,25,by=5),'>30'), breaks = Breaks, low = 'green', mid = 'yellow', high = 'red', midpoint = 15, name = 'Years to Positive NPB') +
                       #                  scale_fill_brewer(palette = 'Dark2') +
                       xlab('Size of Reserve') + 
                       ylab('NPB)') +
                       scale_x_continuous(labels = percent) + 
                       Theme + 
                       theme(legend.position = 'top'))
  return(static_comp_plot) 
}