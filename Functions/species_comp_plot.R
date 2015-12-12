species_comp_plot_fun <- function( PlotData, Theme)
{
  species_comp_plot<- (ggplot(data = PlotData, aes(Species1,Species2,fill = round(Percent_of_best_Utility,2)),color='white') + 
                       geom_tile(color = 'black') + 
                      geom_text(aes(label = ReserveStrategy), size = 3) +
                       scale_fill_gradient(name = '% of Opt. Utility',labels = percent, low = 'blue', high = 'green') + 
                       Theme + 
                       xlab('Reserve Used') + 
                       ylab('Species Managed') + 
                       theme(axis.line = element_blank(),axis.ticks = element_blank(),
                             axis.text.x = element_text(angle = 35, vjust = 0.9,hjust = 0.9)
                      ,legend.title = element_text(size = 6, face = 'plain'),
                      legend.text = element_text(size = 6),
                      legend.position = 'right'))

return(species_comp_plot)
}