coral_static_NPB_plot_fun <- function(PlotData, Theme)
{
  Breaks <- seq(10, 30, by = 5)

  PlotData$TimeToNPB[PlotData$TimeToNPB > 30] <- 30

  PlotData$FinalNPB[PlotData$ReserveSize == 0] <- 0

  leg.title = expression(paste('% ', Delta, ' Cumu. Yields', sep = ''), size = 12)


  static_NPB_plot <- (
    ggplot(data = PlotData  ,
           aes(
             x = ReserveSize, y = 1000*FinalNPB, fill = Final_Perc_Biomass
           )) +
      geom_bar(
        stat = 'identity',
        position = 'dodge',
        color = 'black',
        alpha = 0.85
      ) +
      geom_point(
        aes(size = TimeToNPB),
        shape = 21,
        fill = 'grey71',
        alpha = 0.7
      ) +
      scale_size_continuous(
        breaks = Breaks,
        labels = c(seq(10, 25, by = 5), '>30'),
        name = 'Years to Benefits',
        range = c(1, 5)
      ) +
      facet_grid(.~Species, scales = 'fixed') +
      geom_hline(aes(yintercept = 0), linetype = 'longdash') +
      scale_fill_gradient(
        low = 'blue',
        high = 'green',
        name = '% Increase in Fish',
        labels = percent
      ) +
      #                         scale_fill_gradient2(labels = c(seq(10,25,by=5),'>30'), breaks = Breaks, low = 'green', mid = 'yellow', high = 'red', midpoint = 15, name = 'Years to Positive NPB') +
      #                  scale_fill_brewer(palette = 'Dark2') +
      scale_x_continuous(name ='% of Fishery in Marine Reserve', labels = percent) +
      scale_y_continuous(name = 'Net Community Benefit',labels = dollar) +
      Theme +
      theme(
        legend.position = 'right',
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8, face = 'plain'),
        strip.text = element_text(size = 8)
      )
  )
  return(static_NPB_plot)
}