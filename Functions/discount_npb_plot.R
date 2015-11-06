discount_npb_plot_fun<- function(OptRun,Theme)
{
  
  discount_npb_plot<- (ggplot(data = OptRun ,aes(x = Year, y = NPB))+geom_line(aes(group = DiscountRate),size = 1, alpha = 0.7)+
                         geom_point(aes(size = CurrentReserve,fill = factor(DiscountRate)), shape = 21, alpha = 0.4) + 
                         #                          geom_point(aes(shape = factor(DiscountRate),size = CurrentReserve),fill = 'grey95', alpha = 0.6) + 
                         
                         #                          scale_shape_discrete(name = 'Discount Rate') +
                         #                          scale_fill_manual(name = 'Discount Rate',values = c('green','yellow','red')) +
                         geom_hline(aes(yintercept = 0), linetype = 'longdash', color = 'grey2') + 
                         facet_wrap(~Species) +
                         xlim(c(0,30))+
                         ylim(NA,3000)+
                         facet_grid(.~DiscountRate, labeller = label_bquote(r: .(x))) +
                         Theme + 
                         scale_size_continuous(labels = percent, name = '% Reserve', range = c(2,6)) +
                         #                          scale_fill_manual(guide = FALSE,values = c('green','yellow','red')) +
                         scale_fill_manual(guide = FALSE,values = c('green','red')) +
                         
                         theme(legend.position = 'top'))
  
  return(discount_npb_plot)
}