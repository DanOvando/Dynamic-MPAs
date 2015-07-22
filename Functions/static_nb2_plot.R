static_NB2_plot_fun <- function(PlotData,Theme)
{
  
  Breaks <- seq(10,30,by=5)
  
  PlotData$TimeToNPB[ PlotData$TimeToNPB>30] <- 30
  
  static_NB2_plot <- (ggplot(data = PlotData  ,
                            aes(x = ReserveSize,y = FinalNB,fill = FinalNPB)) +
                       geom_bar(stat = 'identity', position = 'dodge', color = 'black', alpha = 0.85) +
                       facet_wrap(~Species,scales = 'fixed') + 
                       geom_hline(aes(yintercept = 0), linetype = 'longdash') + 
                       scale_fill_gradient2(name = 'NPB r=0.2',low = 'red', mid = 'yellow', high = 'green', midpoint = 0) +
                       #                         scale_fill_gradient2(labels = c(seq(10,25,by=5),'>30'), breaks = Breaks, low = 'green', mid = 'yellow', high = 'red', midpoint = 15, name = 'Years to Positive NPB') +
                       #                  scale_fill_brewer(palette = 'Dark2') +
                       xlab('Size of Reserve') + 
                       ylab('NPB r=0') +
                       scale_x_continuous(labels = percent) + 
                       Theme + 
                       theme(legend.position = 'top'))
  return(static_NB2_plot) 
}