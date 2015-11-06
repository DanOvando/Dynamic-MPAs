Target_Run_By_Discount <- function(Discount,Runs,npb_focus, LoanDiscount, target_reserve)
{
  Fleet$YieldDiscount <- Discount
  show(Discount)
  
  Processed <- ProcessNuts(ReserveResults = Runs,Fleet = Fleet,Scale_Yields = T,npb_focus = npb_focus)
  
  DiscRun <- Processed$ReserveResults
  
  DiscRun$DiscountRate <- Discount
  
  return(BestRun = filter(DiscRun, FinalReserve == target_reserve))
  
}