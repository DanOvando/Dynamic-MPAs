static_NPB_plot_fun <- function(PlotData,Theme)
{
  
  
  static_NPB_plot <- (ggplot(data = PlotData  ,
                           aes(x = FinalReserve,y = NPB,fill = NPB)) +
                      geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
                      facet_wrap(~Species,scales = 'fixed') + 
                      geom_hline(aes(yintercept = 0), linetype = 'longdash') + 
                      scale_fill_gradient2(low = 'red', mid = 'yellow', high = 'darkgreen', midpoint = 0, name = 'Relative Yields') +
                      #                  scale_fill_brewer(palette = 'Dark2') +
                      xlab('Size of Reserve') + 
                      scale_x_continuous(labels = percent) + 
                      Theme + 
                      theme(legend.position = 'none'))
 return(static_NPB_plot) 
}