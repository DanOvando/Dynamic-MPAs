# Function to run one iteration of Species, MPA, FIshing, etc. combination --------------------------------------------------------------

RunGridReserve <- function(r,RunMatrix,BasePatches,Populations,BatchFolder,TimeToRun)
{
  
  # Set up Run ---------------------------------------------------------------
  show(r)
  Run <- RunMatrix[r,]
  
  for (d in 1:dim(Run)[2]) {
    eval(parse(text = paste(colnames(Run)[d],'<- I(Run[d])',sep = '')))
  }
  
  TotalStorage<- as.data.frame(matrix(NA,nrow= TimeToRun,ncol=16))
  
  colnames(TotalStorage)<- c('Species','Run','CurrentReserve','FinalReserve','Intercept','Slope','f','Frate','FracNTZ','Year','Yield','Biomass','Numbers','SQYield','SQNumbers','SQBiomass')
  
  Patches<- BasePatches
  
  short_species<- gsub(' ','',Species,fixed = T)
  
  lh<- eval(parse(text=paste('Populations$',short_species,sep='')))
  
  EQPopulation<- lh$UnfishedPopulation
  
  UnfishedPopulation<- EQPopulation$FinalNumAtAge #Unfished numbers at age
  
  Fmsy<- lh$Fmsy
  
  if (is.null(lh$FatTarget) | is.numeric(FLevel))
  {
  FatTarget<- optimize(log(2*Fmsy$par),f=FindReferencePoint,Target='BvBmsy',TargetValue=FLevel,lower=-10,upper=4,Species=Species,lh=lh,UnfishedPopulation=UnfishedPopulation,Patches=Patches) #Find F that results in target B/Bmsy
  
  FatTarget$par<- exp(FatTarget$minimum)
  }
  if (is.null(lh$FatTarget)==F & is.numeric(FLevel)==F) {FatTarget<- lh$FatTarget}
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
  
  BaseConditions$Yield<- round(BasePop$Performance$Yields[length(BasePop$Performance$Yields)],2)
  
  BaseConditions$Biomass<- round(sum(lh$WeightAtAge %*% BasePop$FinalNumAtAge),2)
  
  BaseConditions$Numbers<- round(sum(BasePop$FinalNumAtAge),2)
  
  TempPop<- StartPop
  
  Patches<- BasePatches
   
  # Run Simulation ----------------------------------------------------------

  MPAParams<- c(Run$Intercept,Run$Slope)
  
  for (y in 1:TimeToRun)
  {
    
#     CurrentMPA<- MPA[y]

    CurrentMPA<- MPAFunction(MPAParams,y,Run$ReserveSize,EvalTime=TimeToRun)
#     
#     if (y>length(MPA))
#     {
#       #             CurrentMPA<- MPAs[dim(MPAs)[2],m]
#       CurrentMPA<-  OptNTZSize$par
#       
#     }
    
    Patches<- AssignNTZ(CurrentMPA,ReservePosition,BasePatches)
    
    FTemp<- FatTarget$par 
    
    FVec<- MoveFleet(FTemp,CurrentMPA,FleetSpill,0)
    
#     if (ReserveStrategy=='CatchShareEqNTZ'& y>1)
#     {
#       #             FVec<- MoveFleet(FTemp,CurrentMPA,0,0)            
#       FVec<- MoveFleet(Fmsy$par,CurrentMPA,0,0)            
#       
#     }
#     
    PassPop<- GrowPopulation(TempPop,FVec,1,0,'TEST',Species=Species,lh=lh,Patches=Patches,FigureFolder=FigureFolder) #grow the new population
    
    TempPop<- PassPop$FinalNumAtAge
    
    TotalStorage[y,]<- data.frame(as.character(Species),r,CurrentMPA,Run$ReserveSize,Run$Intercept,Run$Slope,FLevel,FatTarget$par,CurrentMPA,y,PassPop$Performance$MeanYield,PassPop$Performance$MeanBiomass,PassPop$Performance$MeanNumbers, BaseConditions$Yield,BaseConditions$Numbers, 
                                  BaseConditions$Biomass,stringsAsFactors=F)
  } #Close TimeToRun loop   
  write.table(paste(round(100*(r/dim(RunMatrix)[1])), '% Done With Nuts Runs',sep=''), file = 'NutsProgress.txt', append = TRUE, sep = ";", dec = ".", row.names = FALSE, col.names = FALSE)
  
  return(TotalStorage)
} #Close function






