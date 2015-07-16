discount_npb_plot_fun<- function(OptRun,Theme)
{
  
  discount_npb_plot<- (ggplot(data = OptRun ,aes(x = Year, y = NPB))+geom_line(aes(group = DiscountRate),size = 1.2)+
                         geom_point(aes(fill = factor(DiscountRate),size = CurrentReserve), alpha = 0.8, shape = 21) + 
                         scale_fill_manual(name = 'Discount Rate',values = c('blue','yellow','green')) +
                         geom_hline(aes(yintercept = 0), linetype = 'longdash', color = 'grey2') + 
                         Theme + 
                         scale_size_continuous(labels = percent, name = '% Reserve') + 
                         theme(legend.position = 'right'))
  
  return(discount_npb_plot)
}