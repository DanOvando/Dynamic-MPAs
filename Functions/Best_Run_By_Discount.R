Best_Run_By_Discount <- function(Discount,Runs,Alpha)
{
  Fleet$YieldDiscount <- Discount
  show(Discount)
  
  Processed <- ProcessNuts(ReserveResults = Runs,Fleet = Fleet,Scale_Yields = T,long_term_objective = Alpha)
  
  DiscRun <- Processed$ReserveResults
  
  DiscRun$DiscountRate <- Discount
  
  return(BestRun = filter(DiscRun, BestRun == T))
  
}