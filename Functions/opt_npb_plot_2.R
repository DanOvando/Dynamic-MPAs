opt_npb_plot_fun2<- function(OptRun,OptStaticRun,Theme)
{
  
  leg.title = expression(paste('% ',Delta,' Yields', sep = ''), size = 12)
  
  opt_npb_plot<- (ggplot(data = OptRun ,aes(x = Year, y = NPB))+geom_line(size = 1.5)+
                    geom_point(aes(fill = s_Balance-1,size = CurrentReserve), shape = 21, alpha = 0.85) +
                    geom_line(data = OptStaticRun,aes(Year,NPB),color = 'red', alpha = 0.5, linetype = 'dashed',size = 1.5) +  
#                     geom_point(data = OptStaticRun,aes(Year,NPB,size = CurrentReserve), alpha = 0.7,shape = 21) +  
                    facet_wrap(~Species,scales = 'free_y') +
                    scale_fill_gradient2(labels = percent, name = leg.title,low = 'red', mid = 'yellow', high = 'green',midpoint = 0)
                  + geom_hline(aes(yintercept = 0), linetype = 'longdash', color = 'grey2') + 
                    Theme + 
                    scale_size_continuous(labels = percent, name = '% Reserve', range = c(1,6)) + 
                    theme(legend.position = 'right', text = element_text(size = 12), 
                          strip.text.x = element_text(size = 8), legend.title= 
                            element_text(size = 8, face = 'plain'),legend.text= 
  element_text(size = 8)))
  
  return(opt_npb_plot)
}