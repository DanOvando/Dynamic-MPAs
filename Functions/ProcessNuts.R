ProcessNuts <- function(ReserveResults,Fleet,Scale_Yields,long_term_objective)
{
  
  if (Scale_Yields == T)
  {
    YieldFactor<- ReserveResults$Yield/ReserveResults$SQYield
    
    ReserveResults$SQYield <- 100
    
    ReserveResults$Yield <- ReserveResults$SQYield * YieldFactor
    
    ReserveResults$YieldBalance <- ReserveResults$Yield - ReserveResults$SQYield
    
  }
  
  ReserveResults$YieldBalance <- ReserveResults$Yield - ReserveResults$SQYield
  
  ReserveResults <- ReserveResults %>%
    group_by(Run) %>% 
    arrange(Year) %>%
    mutate(PresentYield=Yield*(1+Fleet$YieldDisc)^-(Year-1),PresentBalance=(Yield-SQYield)*(1+Fleet$YieldDisc)^-(Year-1)) %>%
    mutate(NPB=cumsum(PresentBalance),NPY=cumsum(PresentYield),RequestedLoan = sum(PresentBalance[YieldBalance<0])) %>%
    ungroup() %>%
    group_by(Species) %>% 
    mutate(s_Yield = Yield/max(Yield,na.rm=T),s_PresentYield=s_Yield*(1+Fleet$YieldDisc)^-(Year-1)) %>%
    ungroup() %>%
    group_by(Run) %>%
    arrange(Year) %>%
    mutate(net_yields = cumsum(Yield), s_NPY=cumsum(s_PresentYield),s_Balance=(Yield/SQYield),
           s_NPB=cumsum(s_Balance*(1+Fleet$YieldDisc)^-(Year-1))) %>%
    ungroup() %>%
    group_by(Species) %>%
    mutate(Utility = long_term_objective * (net_yields / max(net_yields)) + (1 - long_term_objective) * (NPY/max(NPY)) ) %>%
    ungroup()
  
  
  
  OptimalRun<- filter(ReserveResults,Year==max(Year)) %>%
    group_by(Species) %>%
    summarize(OptNPB = Run[NPB == max(NPB)][1], OptEQ = Run[net_yields == max(net_yields)][1], OptUtility = Run[Utility == max(Utility)][1])
  
  ReserveResults$Scenario <- with(ReserveResults,paste(FinalReserve,Intercept,Slope,sep = '-'))
  
  UnifiedRuns<- filter(ReserveResults,Year==max(Year)) %>%
    group_by(Scenario) %>%
    summarize(TotalUtility = sum(Utility))
  
  BestUnifiedScenario <- UnifiedRuns$Scenario[UnifiedRuns$TotalUtility == max(UnifiedRuns$TotalUtility)]
  
  ReserveResults$BestUnifiedRun <- ReserveResults$Scenario %in% BestUnifiedScenario
  
  ReserveResults$BestRun <- ReserveResults$Run %in% OptimalRun$OptUtility
  
  ResSummary <- ReserveResults %>%
    group_by(Run) %>%
    summarize(TimeToNPB = which(NPB >= 0)[1],
              Species = unique(Species), 
              BestRun = unique(BestRun), 
              FinalNPB = last(NPB),
              FinalUtility = last(Utility), 
              ReserveSize = mean(FinalReserve), 
              Intercept = mean(Intercept), 
              Slope = mean(Slope),
              NegativeYields = Discount(pmin(0,YieldBalance),0,length(YieldBalance))$NPV,
              StatusQuoYields = Discount(SQYield * (YieldBalance <= 0),0,length(YieldBalance))$NPV,
              PriceInc = 100 * (StatusQuoYields / (StatusQuoYields + NegativeYields) - 1),
              AvailableSurplus = Discount(pmax(0, YieldBalance),0,length(YieldBalance))$NPV,
              LoanNegativeYields = Discount(pmin(0,YieldBalance[1:LoanTime]),0,length(YieldBalance[1:LoanTime]))$NPV,
              LoanStatusQuoYields = Discount((SQYield * (YieldBalance <= 0))[1:LoanTime],0,length(YieldBalance[1:LoanTime]))$NPV,
              LoanAvailableSurplus = Discount(pmax(0, YieldBalance[1:LoanTime]),0,length(YieldBalance[1:LoanTime]))$NPV)
  
  RunNames <- unique(ResSummary$Run)
  for (i in seq_len(length(RunNames)))
  {
    Where <- ResSummary$Run ==  RunNames[i]
    ResSummary$MaxInterestRate[i] <- 100 * exp( optim(-4, FindMaxInterestRate, LoanTime = LoanTime, LoanAmount = -ResSummary$LoanNegativeYields[Where],
                                                      Surplus = ResSummary$LoanAvailableSurplus[Where],lower = -10, upper = 10,method = 'Brent')$par)
  }
  
  ResSummary$TimeToNPB[is.na(ResSummary$TimeToNPB)] <- max(ReserveResults$Year)
  
  SpeciesList<- unique(ReserveResults$Species)
  
  species_comparison<- as.data.frame(matrix(NA,nrow=length(SpeciesList)*length(SpeciesList),ncol=7))
  
  colnames(species_comparison)<- c('Species1','Species2','NPB','Utility','TimeToNPB','PrinceInc','MaxInterestRate')
  
  cc <- 0
  
  ResSummary$RunDescription<- paste(ResSummary$ReserveSize,ResSummary$Intercept,ResSummary$Slope,sep='-')
  
  for (i in 1:length(SpeciesList))
  {
    Where<- ResSummary$Species == SpeciesList[i] & ResSummary$BestRun == T
    
    species_opt_res<- ResSummary$RunDescription[Where]
    
    for (j in 1:length(SpeciesList))
    {
      cc <- cc+1
      
      Where2<- ResSummary$Species == SpeciesList[j] & ResSummary$RunDescription == species_opt_res
      
      species_comparison[cc,]<- data.frame(SpeciesList[i], SpeciesList[j], 
                                           ResSummary$FinalNPB[Where2],
                                           ResSummary$FinalUtility[Where2],
                                           ResSummary$TimeToNPB[Where2], 
                                           ResSummary$PriceInc[Where2], 
                                           ResSummary$MaxInterestRate[Where2],stringsAsFactors = F)
      
    }
  }
  
  species_comparison<- species_comparison %>%
    ungroup() %>%
    group_by(Species2) %>%
    mutate(Percent_of_best_NPB = sign(max(NPB)) * (NPB / max(NPB) - 1 ), 
           Percent_of_best_Utility = sign(max(Utility)) * (Utility / max(Utility)) ) %>%
    arrange(Species2)
  
  return(list(ReserveResults = ReserveResults, species_comparison = species_comparison,ResSummary = ResSummary))
}
