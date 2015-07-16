opt_biomass_plot_fun<- function(OptRun,Theme)
{
  opt_biomass_plot<- (ggplot(data = OptRun ,aes(x = Year, y = Biomass))+geom_line(size = 1.2)+
                    geom_point(aes(fill = NPB,size = CurrentReserve), shape = 21, alpha = 0.85) + 
                    facet_wrap(~Species,scales = 'fixed') +
                    scale_fill_gradient2(name = 'NPB',low = 'red', mid = 'yellow', high = 'green',midpoint = 0)
                  + geom_hline(aes(yintercept = 0), linetype = 'longdash', color = 'grey2') + 
                    Theme + 
                    scale_size_continuous(labels = percent, name = '% Reserve') + 
                    theme(legend.position = 'right'))
  
  return(opt_biomass_plot)  

}