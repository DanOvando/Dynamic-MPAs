static_NPB_plot_fun <- function(PlotData,Theme)
{
  
  Breaks <- seq(10,30,by=5)
  
  PlotData$TimeToNPB[ PlotData$TimeToNPB>30] <- 30
  
  PlotData$FinalNPB[PlotData$ReserveSize == 0] <- 0
  
  leg.title = expression(paste('% ',Delta,' Cumu. Yields', sep = ''), size = 12)
  
  
  static_NPB_plot <- (ggplot(data = PlotData  ,
                             aes(x = ReserveSize,y = FinalNPB,fill = Final_PercNB)) +
                        geom_bar(stat = 'identity', position = 'dodge', color = 'black', alpha = 0.85) +
                        geom_point(aes(size = TimeToNPB), shape = 21, fill = 'grey71', alpha = 0.7) +
                        scale_size_continuous(breaks = Breaks, labels = c(seq(10,25,by=5),'>30'),
                                              name = 'Years to NPB>0', range = c(1,5)) + 
                        facet_wrap(~Species,scales = 'fixed') + 
                        geom_hline(aes(yintercept = 0), linetype = 'longdash') + 
                        scale_fill_gradient(low = 'blue', high = 'green', name = leg.title,
                                            labels = percent) +
#                         scale_fill_gradient2(labels = c(seq(10,25,by=5),'>30'), breaks = Breaks, low = 'green', mid = 'yellow', high = 'red', midpoint = 15, name = 'Years to Positive NPB') +
                        #                  scale_fill_brewer(palette = 'Dark2') +
                        xlab('Size of Reserve') + 
                        ylab('NPB') +
                        scale_x_continuous(labels = percent) + 
                        Theme + 
                        theme(legend.position = 'right', legend.text = element_text(size = 8),
                              legend.title = element_text(size = 8,face = 'plain'), strip.text = element_text(size = 8)))
  return(static_NPB_plot) 
}