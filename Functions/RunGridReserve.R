#' Function to run one iteration of Species, MPA, FIshing, etc. combination
#' 
#' \code{RunGridReserve} runs each iteration of NUTS model 
#' defined by the RunMatrix
#' @param r the run to use
#' @param RunMatrix a matrix with rows for unique runs
#' and columns for parameter values to be used
#' in a particular run
#' @param BasePatches the default patch structure
#' @param Populations a list of the equilibrium pop
#' for each species
#' @param BatchFolder the folder to store outputs
#' @param TimeToRun the number of years to run
#' 
RunGridReserve <- function(r,RunMatrix,BasePatches,Populations,BatchFolder,TimeToRun)
{
  
  # Set up Run ---------------------------------------------------------------
  Run <- RunMatrix[r,]
  
  for (d in 1:dim(Run)[2]) { # Unclear actually what this does
    eval(parse(text = paste(colnames(Run)[d],'<- I(Run[d])',sep = '')))
  }
  
  store_vars <- c('Species','Run','discount_rate','CurrentReserve','FinalReserve','Intercept','Slope','f','Frate','FracNTZ','Year','Yield','Biomass','Numbers','SQYield','SQNumbers','SQBiomass')
  
  TotalStorage<- as.data.frame(matrix(NA,nrow= TimeToRun,ncol=17)) #storage space
  
  colnames(TotalStorage) <- store_vars 
  
  Patches<- BasePatches #revert patches to basepatches
  
  short_species<- gsub(' ','',Species,fixed = T)
  
  lh<- eval(parse(text=paste('Populations$',short_species,sep=''))) #store life history
  
  EQPopulation<- lh$UnfishedPopulation
  
  UnfishedPopulation<- EQPopulation$FinalNumAtAge #Unfished numbers at age
  
  Fmsy<- lh$Fmsy

  if (is.null(lh$FatTarget) | is.numeric(FLevel)) # Find F/Fmsy that acchieves target B/Bmsy
  {
    FatTarget<- optimize(log(2*Fmsy$par),f=FindReferencePoint,Target='BvBmsy',TargetValue=FLevel,lower=-10,upper=4,Species=Species,lh=lh,UnfishedPopulation=UnfishedPopulation,Patches=Patches) #Find F that results in target B/Bmsy
    
    FatTarget$par<- exp(FatTarget$minimum)
  }
  if (is.null(lh$FatTarget)==F & is.numeric(FLevel)==F) {FatTarget<- lh$FatTarget}
  
  # Fish population to target B/Bmsy
  BatTargetPopulation<- GrowPopulation(UnfishedPopulation, rep(FatTarget$par,NumPatches),'EQ',0,'BatTarget Run',Species=Species,lh=lh,Patches=Patches,FigureFolder=FigureFolder)
  
  EvalTime<- OptTime #Time span to evaluate results on
  
  RunTime<- 'Custom'
  
  FScenarios<- c(FatTarget$par) #load in fishing scenarios
  
  BaseConditions<- as.data.frame(matrix(NA,nrow= length(FScenarios),ncol=3))
  colnames(BaseConditions)<- c('Yield','Biomass','Numbers')
  
  OptimalConditions<- as.data.frame(matrix(NA,nrow= length(FScenarios),ncol=3))
  colnames(OptimalConditions)<- c('Yield','Biomass','Numbers')
  
  BasePop<- GrowPopulation(EQPopulation$FinalNumAtAge,FatTarget$par,'EQ',0,paste('FvFmsy is',round(FatTarget$par/Fmsy$par,2)),Species=Species,lh=lh,Patches=Patches,FigureFolder=FigureFolder)
  
  StartPop<- BasePop$FinalNumAtAge
  
  # Calculate base conditions ----
  
  BaseConditions$Yield<- round(BasePop$Performance$Yields[length(BasePop$Performance$Yields)],2)
  
  BaseConditions$Biomass<- round(sum(lh$WeightAtAge %*% BasePop$FinalNumAtAge),2)
  
  BaseConditions$Numbers<- round(sum(BasePop$FinalNumAtAge),2)
  
  TempPop<- StartPop
  
  Patches<- BasePatches
  
  # Run Simulation ----------------------------------------------------------
  
  MPAParams<- c(Run$Intercept,Run$Slope)
  
  for (y in 1:TimeToRun)
  {
    
    CurrentMPA<- MPAFunction(MPAParams,y,Run$ReserveSize,EvalTime=TimeToRun) #set MPA in current year
    
    if (y==1)
    {
      CurrentMPA<- 0
    }
    
    Patches<- AssignNTZ(CurrentMPA,ReservePosition,BasePatches) #assign MPA
    
    FTemp<- FatTarget$par 
    
    FVec<- MoveFleet(FTemp,CurrentMPA,FleetSpill,0) # Distribute fishing effort
    # Grow population
    PassPop<- GrowPopulation(TempPop,FVec,1,0,'TEST',Species=Species,lh=lh,Patches=Patches,FigureFolder=FigureFolder) #grow the new population
    
    TempPop<- PassPop$FinalNumAtAge

    TotalStorage[y,]<- data.frame(as.character(Species),r,as.numeric(DiscountRate),CurrentMPA,Run$ReserveSize,Run$Intercept,Run$Slope,FLevel,FatTarget$par,CurrentMPA,y,PassPop$Performance$MeanYield,PassPop$Performance$MeanBiomass,PassPop$Performance$MeanNumbers, BaseConditions$Yield,BaseConditions$Numbers, 
                                  BaseConditions$Biomass,stringsAsFactors=F)
  } #Close TimeToRun loop   
  write.table(paste(round(100*(r/dim(RunMatrix)[1])), '% Done With Nuts Runs',sep=''), file = 'NutsProgress.txt', append = TRUE, sep = ";", dec = ".", row.names = FALSE, col.names = FALSE)
  
  return(TotalStorage)
} #Close function






