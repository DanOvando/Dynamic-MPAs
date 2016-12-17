opt_npb_plot_fun<- function(OptRun,Theme)
{

  opt_npb_plot<- (ggplot(data = OptRun ,aes(x = Year, y = NPB))+geom_line(size = 1.5)+
   geom_point(aes(fill = s_Balance-1,size = CurrentReserve), shape = 21, alpha = 0.85) +
   facet_wrap(~Species,scales = 'free') +
   scale_fill_gradient2(labels = percent, name = '% Change from SQ Yields',low = 'red', mid = 'yellow', high = 'green',midpoint = 0)
 + geom_hline(aes(yintercept = 0), linetype = 'longdash', color = 'grey2') +
   Theme +
   scale_size_continuous(labels = percent, name = '% Reserve', range = c(2,8)) +
   theme(legend.position = 'right'))

  return(opt_npb_plot)
}