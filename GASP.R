
####### G.A.S.P. (General age structured population model... sort of) ######

#Dan Ovando
#4.29/13
#Summary: This suite of functions run a general age structured fishery model. Will prepare a short summary shortly to use the set. 

GrowPopulation<- function(InitialPopulation,FishingPressure,Time,MakePlots,GroupFigName)
{
            # InitialPopulation<- 100
              # FishingPressure<- c(0,.4)           
              # Time=3
              # MakePlots=0
            # GroupFigName= 'Test'
       
       
  RelativePatchSizes <- Patches$PatchSizes/sum(Patches$PatchSizes)
  # InitialPopulation<- 1000
  # FishingPressure<- c(0,0)
  # Time='EQ'
  # MakePlots=1
  # GroupFigName= 'Test'
  
          # InitialPopulation<- UnfishedPopulation
          # FishingPressure<- 4*Fmsy$par
           # Time='EQ'
         # MakePlots=0
         # GroupFigName= 'DIE'
  
    # InitialPopulation<-  TempPop
    # FishingPressure<- FTemp
    # Time='EQ'
    # MakePlots=1
    # GroupFigName= 'MOVE'    
  
  
  EQMarker<- 0 # Blank marker to start population loop, leave at 0
  if (Time=='EQ')
  {
    TempTime<- 2000
    # StopTime<- TempTime-1
  }
  if (Time !='EQ')
  {
    TempTime<- Time+2 #TempTime is used to pass the amount of time of the population run
    StopTime<- Time+1
  }
  
  if (is.vector(FishingPressure)) #Reformat the fishing pressure to be by patch over time
  {
    if (length(FishingPressure)==1) #If there's one constant fishing pressure
    {
      FishingPressure<- rep(FishingPressure,(TempTime+1)*NumPatches)
    }
    else  #If you have one year of patch specific fishing pressure
    {
      FishingPressure<- rep(FishingPressure,TempTime+1)
    }
    FishingPressure<- matrix(FishingPressure,nrow=TempTime+1,ncol=NumPatches,byrow=T)
    
  }
    
  PopArray<- array(0,c(TempTime+1,lh$MaxAge,NumPatches)) #blank storage for population array
  
  PopArray[1,,]<- InitialPopulation #Seed population array with initial population
  
  WeightArray<- array(0,c(TempTime+1,lh$MaxAge,NumPatches)) #blank storage for population array
  
  WeightArray[1,,]<- PopArray[1,,] * WeightAtAge
  
  StoreFishedWeight<- as.data.frame(matrix(0,nrow=TempTime+1,ncol=NumPatches)) #Blank for fishing yields
  
  t<- 1 #leave at 1
  
  NumberOfEggs<- Eggs(PopArray[1,,],FecundityAtAge,MaturityAtAge) #Produce eggs by patch    

  while(t<TempTime) #Loop as long as you want
  {    
    
        
    RelocatedEggs<- MoveEggs(NumberOfEggs,'Common') #Move eggs
    
    RecruitsPerPatch<- Recruits(RelocatedEggs,lh$DDForm,'Random',c(lh$RecDevMean,lh$RecDevSTD)) #Calculate recruits

    SelectivityAtAge<- FishingSelectivity(LengthAtAge,Fleet$s50,Fleet$s95,1000)	#Calculate selectivity curve
    
    FishingAtAge<- DistFishingAtAge(FishingPressure[t,],SelectivityAtAge) #Calculate fishing mortality at age
    
    Survival<- exp(-(FishingAtAge+lh$m)) #calculate survival from each age class to the next
        
    FishingYields<- GoFish(FishingAtAge,Survival,PopArray[t,,])     
        
    StoreFishedWeight[t,]<- FishingYields$Total 

    # PopArray[t,,]<- MoveAdults(PopArray[t,,],lh$Range,lh$MoveType)
    
    PlusGroup<- PopArray[t,lh$MaxAge,]*Survival[lh$MaxAge,] #The numbers in the plus group that survive
    	    
    PopArray[t+1,2:lh$MaxAge,]<- PopArray[t,1:(lh$MaxAge-1),] * Survival[1:(lh$MaxAge-1),] #Grow to the next age class
    
    PopArray[t+1,lh$MaxAge,]<- PlusGroup #Add in the survivors of the plus group
    
    PopArray[t+1,,]<- MoveAdults(PopArray[t+1,,],lh$Range,lh$MoveType)

    NumberOfEggs<- Eggs(PopArray[t+1,,],FecundityAtAge,MaturityAtAge) #Produce eggs by patch    

    PopArray[t+1,1,]<- RecruitsPerPatch #Add recruits to population 
    
    PopArray[PopArray<1]<- 0

	# show(MoveAdults(PopArray[t+1,,],lh$Range,lh$MoveType))
   

    WeightArray[t+1,,]<- PopArray[t+1,,] * WeightAtAge
        
     PopChange<- sum(PopTolerance< abs((colSums(WeightArray[t+1,,])-colSums(WeightArray[t,,])))) #Measure change in 

    # PopChange<- sum(PopTolerance< abs((colSums(WeightArray[t+1,,])/colSums(WeightArray[t,,]))-1)) #Measure change in 


    t<- t+1
    if (t>1900 & PopChange>0)
	
    {
    	TempTime<- t
    	StopTime<- t-1
    	Time<- StopTime-1
    	   if (length(dim(PopArray))>=3) #Pain in the ass step because R is a dick
   {
     PopTrajectory<- apply(PopArray,3,rowSums)
     WeightTrajectory<- apply(WeightArray,3,rowSums)
   }
   if (length(dim(PopArray))<3)
   {
     PopTrajectory<- t(as.matrix(colSums(PopArray)))
     WeightTrajectory<- t(as.matrix(colSums(WeightArray)))
   }
    	warning(paste('Population never reached EQ, F is ',FishingPressure))
     quartz()
     matplot(WeightTrajectory,type='l',lwd=4,xlab='Time',ylab='B/Bmsy',col=terrain.colors(2*NumPatches)[1:NumPatches],bty='n')

    }
    
    
    if (Time=='EQ' & PopChange==0 & t>5) #If you're running to EQ, stop once population isn't changing
    {
      StopTime<- t-1
      if (StopTime<2)
      {
        show(paste('Time is',t))
        show(FishingPressure[t,])
      }
      t<- TempTime+1
      Time<- StopTime-1
      PopArray<- PopArray[1:StopTime,,]
      WeightArray<- WeightArray[1:StopTime,,]
      StoreFishedWeight<- StoreFishedWeight[1:StopTime,]
    }
    
  } #close time loop
  
  
  PopArray<- (PopArray[2:StopTime,,])
  WeightArray<- WeightArray[2:(StopTime),,]
  StoreFishedWeight<- StoreFishedWeight[2:(StopTime),]
  
  
  if (length(dim(PopArray))>=3) #Pain in the ass step because R is a dick
  {
    PopTrajectory<- apply(PopArray,3,rowSums)
    WeightTrajectory<- apply(WeightArray,3,rowSums)
  }
  if (length(dim(PopArray))<3)
  {
    PopTrajectory<- t(as.matrix(colSums(PopArray)))
    WeightTrajectory<- t(as.matrix(colSums(WeightArray)))
  }
  
  ### Produce performance metrics ###
  
  Performance<- NULL
  
  # Performance$Yields<- rowSums(StoreFishedWeight[2:dim(StoreFishedWeight)[1],])
  Performance$Yields<- rowSums(StoreFishedWeight)
  
  # Performance$DiscYields<- Discount(StoreFishedWeight[2:dim(StoreFishedWeight)[1],],Fleet$YieldDiscount,Time) #Discount yield and calculate NPV
  Performance$DiscYields<- Discount(StoreFishedWeight,Fleet$YieldDiscount,Time) #Discount yield and calculate NPV
  
  # Performance$DiscNumbers<- Discount(PopTrajectory[2:dim(PopTrajectory)[1],],Fleet$BiomassDiscount,Time)  
  Performance$DiscNumbers<- Discount(PopTrajectory,Fleet$BiomassDiscount,Time)  
  
  # Performance$DiscBiomass<- Discount(WeightTrajectory[2:dim(WeightTrajectory)[1],],Fleet$BiomassDiscount,Time)  
  Performance$DiscBiomass<- Discount(WeightTrajectory,Fleet$BiomassDiscount,Time)  
  
  Performance$YieldInstability<- sd((Performance$Yields[2:length(Performance$Yields)]-Performance$Yields[1:(length(Performance$Yields)-1)])/Performance$Yields[1:(length(Performance$Yields)-1)])
  
  TempTotalPop<- rowSums(PopTrajectory)
  Performance$MeanNumbers<- mean(TempTotalPop)
  
  TempTotalBiomass<- rowSums(WeightTrajectory)
  Performance$MeanBiomass<- mean(TempTotalBiomass)
  
  TempTotalYields<- rowSums(StoreFishedWeight)
  Performance$MeanYield<- mean(TempTotalYields)
  
  Performance$MeanYieldChange<- mean((Performance$Yields[2:length(Performance$Yields)]-Performance$Yields[1:(length(Performance$Yields)-1)])/Performance$Yields[1:(length(Performance$Yields)-1)])
  
  if (Time>1)
  {
    FinalNumAtAge<- PopArray[dim(PopArray)[1],,]
    FinalLengthFrequency<- LengthFrequency(PopArray[dim(PopArray)[1],,],LengthAtAge)
  }
  if (Time==1)
  {
    FinalNumAtAge<- PopArray
    FinalLengthFrequency<- LengthFrequency(PopArray,LengthAtAge)
  }
  
  if (lh$CarryingCapacity[1] !=-999)
  {
    Performance$FreqCollapse<- mean(TempTotalPop/sum(lh$CarryingCapacity)<= CollapseThreshold)
  }
  
  
  if (MakePlots==1)
  {   
    FormatFigure(paste(GroupFigName,'Total Numbers (N) Over Time.pdf'))
    matplot(PopTrajectory,type='l',lwd=4,xlab='Time',ylab='Total Numbers',col=terrain.colors(2*NumPatches)[1:NumPatches],bty='n')
    legend('topright',legend=paste('Patch',1:NumPatches),col=terrain.colors(2*NumPatches)[1:NumPatches],lty=1:NumPatches,lwd=4,bty='n')
    dev.off()  
    
    if (lh$Bmsy[1]!=-999)
    {
      FormatFigure(paste(GroupFigName,'BvBmsy Over Time.pdf'))
      matplot(rowSums(WeightTrajectory)/sum(lh$Bmsy),type='l',lwd=4,xlab='Time',ylab='B/Bmsy',col=terrain.colors(2*NumPatches)[1:NumPatches],bty='n')
      legend('topright',legend=paste('Patch',1:NumPatches),col=terrain.colors(2*NumPatches)[1:NumPatches],lty=1:NumPatches,lwd=4,bty='n')
      dev.off()
      
    }    
    FormatFigure(paste(GroupFigName,'Final Length Frequency Bar.pdf'))
    sizemat<- matrix(c(1:NumPatches),nrow=NumPatches,ncol=1)
    sizelay<- layout(mat=sizemat)
    for (n in 1:NumPatches)
    {   
      
     if (sum( FinalLengthFrequency[,n],na.rm=T)>0)
     {
     	FinalLengthFrequency[FinalLengthFrequency==0]<- NA
     }
      
      
      b=hist(FinalLengthFrequency[,n],breaks=15,xlab='Length',ylab='Frequency',col=terrain.colors(2*NumPatches)[n],main=paste('Patch',n),xlim=c(0,lh$Linf*1.2))
    }
    dev.off()
    
    FormatFigure(paste(GroupFigName,'Fishing Yields.pdf'))
    matplot(cbind(StoreFishedWeight,rowSums(StoreFishedWeight)),type='l',lwd=4,xlab='Time',ylab='Landings', col=terrain.colors(2*NumPatches)[1:(NumPatches+1)],bty='n')
    legend('topright',legend=c(paste('Patch',1:NumPatches),'Total'),col=terrain.colors(2*NumPatches)[1:(NumPatches+1)],lty=1,lwd=4,bty='n')
    dev.off()
    
    FormatFigure(paste(GroupFigName,'Biomass Over Time.pdf'))
    matplot(cbind(WeightTrajectory,rowSums(WeightTrajectory)),type='l',lwd=4,xlab='Time',ylab='Biomass', col=terrain.colors(2*NumPatches)[1:(NumPatches+1)],bty='n')
    legend('topright',legend=c(paste('Patch',1:NumPatches),'Total'),col=terrain.colors(2*NumPatches)[1:(NumPatches+1)],lty=1,lwd=4,bty='n')
    dev.off()
    
  } 
  
  return(list(NumatAge=PopArray,FinalNumAtAge=FinalNumAtAge,FishingYields=StoreFishedWeight,TotalPop=PopTrajectory,Performance=Performance)) 
  
#   return(list(NumatAge=PopArray,FinalNumAtAge=FinalNumAtAge,FishingYields=StoreFishedWeight,TotalPop=PopTrajectory,Performance=Performance,FinalLengthFrequency= FinalLengthFrequency)) 
} #Close population growth function

FormatFigure<- function(name)
{
  pdf(file=paste(FigureFolder,name,sep=''),family=Font,pointsize=FontSize)
}

Discount<- function(DataVec,DR,Time)
{
  
  # DataVec<- PopTrajectory
  DiscValues<- DataVec*(1+DR)^-(0:(Time-1))
  NPV<- sum(DiscValues)
  CumNPV<- cumsum(DiscValues)
  return(list(DiscValues= DiscValues,NPV=NPV,CumValues=CumNPV))
}

Length <- function(ages)
{
  Length<-lh$Linf*(1-exp(-1*lh$k*(ages-lh$t0)))
  return(Length)
}

LengthWiError<- function(Lengths,NumDraws) #Produce a distribution 
{
  
  HasNumber<- as.numeric(NumDraws>0)
  NumDraws[NumDraws==0]<- 1
  LenSD<- lh$VBSD*(1+lh$VBErrorSlope*Lengths/lh$Linf)
  LengthAtAgeDist<- matrix(0,nrow<-length(Lengths),ncol=max(1,NumDraws))
  for (l in 1:length(Lengths))
  {
    LengthAtAgeDist[l,]<- HasNumber*(Lengths[l]* rlnorm(NumDraws,mean=lh$VBErrorMean,sd=LenSD[l]))
  }
  
  return(LengthAtAgeDist)
}

LengthFrequency<- function(NumAtAge,MeanLengths)
{
  # NumAtAge<- PopArray[dim(PopArray)[1],,]
  # MeanLengths<- LengthAtAge
  TempNumAtAge<- NumAtAge
  TempNumAtAge[TempNumAtAge==0]<- 1
  FreqMat<- matrix(NA,nrow=max(1,max(colSums(ceiling(TempNumAtAge)))),ncol=NumPatches)
  for (p in 1:dim(NumAtAge)[2])
  {
    for (f in 1:dim(NumAtAge)[1])
    {
      NewLengths<- LengthWiError(MeanLengths[f],ceiling(NumAtAge[f,p]))
#       VecLength<- 1:length(FreqMat[,p])
#       Where<- (is.na(FreqMat[,p]))
#       Where<- VecLength[Where]
#       Where<- VecLength[1]
#       
#        Where<- which(is.na(FreqMat[,p]))[1]
#        FreqMat[Where:((Where-1)+length(NewLengths)),p]<- NewLengths    
     FreqMat[1:length(NewLengths),p]<- NewLengths    

    }
    
  }
  return(FreqMat)	
}

Weight<- function(Lengths,WeightForm)
{
  LengthDist<- LengthWiError(Lengths,1000)
  
  if (WeightForm=='Exponential')
   {
  weight<- (lh$wa*LengthDist^lh$wb)
  }
  if (WeightForm =='Linear')
  {
  	  weight<- (lh$wa+LengthDist*lh$wb)

  }
  if (WeightForm =='Exponential2')
  {
  	  # weight<- exp((lh$wa+log(LengthDist)*lh$wb))
  	  weight<- exp((lh$wa+(LengthDist)*lh$wb))
  }
  if (WeightForm =='WTF')
  {
  	  weight<- 10^((lh$wa+lh$wb*log(LengthDist)))

  }
  
  weight<- rowMeans(weight)
  return(weight)
}  

Fecundity <- function (Data,FecForm)
{ 
  
  if (FecForm=='Length')
  {
  LengthDist<- LengthWiError(Data,1000)
  fecund<- rowMeans(lh$fa*LengthDist^lh$fb)
  }
  if (FecForm=='Weight')
  {
  fecund<- (lh$fa* Data ^lh$fb) #Assume weight is kg
  }
  return(fecund)
}

Maturity <- function(Data,Mode)
{
  
  if (Mode=='Length')
{
  
  
  Data<- Length(Data)
  
  s50<- lh$LengthMa50

  s95<- lh$LengthMa95

  LengthDist<- LengthWiError(Data,1000)
   # mature<-rowMeans(round(1/(1+exp(lh$ma50*LengthDist+lh$ma95)),2))
  
    mature<- rowMeans(round(1/(1+exp(-log(19)*((LengthDist-s50)/(s95-s50)))),2))	
}

if (Mode=='Age')
{
	    s50<- lh$AgeMa50

  		s95<- lh$AgeMa95

	    mature<-(round(1/(1+exp(-log(19)*((Data-s50)/(s95-s50)))),2))	
}
  
   # mature<-rowMeans(round((1+exp(-log(19)*(LengthDist-lh$ma50)/lh$maBeta))^-1,2))	
  return(mature)
}

Eggs <- function(NatAge,Fecund,Mature) 
{  
  eggs <- (Fecund *Mature) %*% (NatAge * lh$SexRatio)  
  return(eggs) 
}

MoveEggs<- function(NumberOfEggs,PoolType)
{
  if(PoolType=='Common')
  {
    eggs<- sum(NumberOfEggs)
    
    if (lh$LarvalChoice==0)
    {
      RelocateEggs<- eggs*(Patches$PatchSizes / sum(Patches$PatchSizes))
    }
    if (lh$LarvalChoice==1)
    {
      RelocateEggs<- eggs*(Patches$PatchSizes / sum(Patches$PatchSizes))*(Patches$HabQuality / sum(Patches$HabQuality))
    }
  }
  return(RelocateEggs)
}

MoveAdults<- function(NatAge,Distance,MoveType) ### Work in progress, place holder for adult movement
{

  if (Distance==0)
  {
    MovedAdults<- NatAge
  }
  
  if (MoveType=='Simple') #Reallllllly basic movement: move X% from patch 1 to patch 2 and vice versa. Distance is betweein 0 and 1
  {
  		     # NatAge<- PopArray[t,,]
  		     # Distance<- lh$Range
  		     # MoveType<- lh$MoveType
    # show(NatAge)
  	RelativePatchSize<- Patches$PatchSizes/sum(Patches$PatchSizes)
  	NumberMoved<- NatAge/matrix((RelativePatchSize*10),dim(NatAge)[1],dim(NatAge)[2],byrow=T)
  	NumberMoved[is.nan(NumberMoved)]<- 0
  	NumberMoved<- NumberMoved *(matrix(rev(RelativePatchSize),dim(NatAge)[1],dim(NatAge)[2],byrow=T))
  	
  	NumberMoved<- pmin(NumberMoved,NatAge)
  	
  	PotentiallyMoved<- colSums(NumberMoved)
  	
  	if (is.numeric(lh$CarryingCapacity))
  	{
  		SpaceAvailable<- lh$CarryingCapacity-colSums(NatAge)
  		
  		WillMoved<- pmin(rev(PotentiallyMoved),SpaceAvailable)
  		
  		ProportionMoved<- WillMove/PotentiallyMoved
  	
  		NumberMoved<- NumberMoved * matrix(rep(ProportionMoved,dim(NumberMoved)[1]),nrow=dim(NumberMoved)[1],ncol=dim(NumberMoved)[2],byrow=T)

  	}
  	
  	NatAge[,1]<- NatAge[,1]+NumberMoved[,2]-NumberMoved[,1]
  	NatAge[,2]<- NatAge[,2]+NumberMoved[,1]-NumberMoved[,2]	
  	MovedAdults<- NatAge
  }

  if (MoveType=='Simple-DO') #Reallllllly basic movement: move X% from patch 1 to patch 2 and vice versa. Distance is betweein 0 and 1
  {
  		     # NatAge<- PopArray[t,,]
  		     # Distance<- lh$Range
  		     # MoveType<- lh$MoveType
    # show(NatAge)  	
  	RelativePatchSize<- Patches$PatchSizes/sum(Patches$PatchSizes)
  	  	
  	NumberMoved<- NatAge/matrix((RelativePatchSize*(1/Distance)),dim(NatAge)[1],dim(NatAge)[2],byrow=T) #Divide the sustem into 1/Distance patches, then determine how many of these hypothetical patches are in each real patch
  	
  	NumberMoved[is.nan(NumberMoved)]<- 0 
  	
  	NumberMoved<- NumberMoved *(matrix((1-RelativePatchSize),dim(NatAge)[1],dim(NatAge)[2],byrow=T)) #Scale the number of adults that decide to settle in the adjacent patch by the size of that patch
  	  	
  	NumberMoved<- pmin(NumberMoved,NatAge) #Prevent more movement than individuals

  	NatAge[,1]<- NatAge[,1]+NumberMoved[,2]-NumberMoved[,1] #Move things around
  	
  	NatAge[,2]<- NatAge[,2]+NumberMoved[,1]-NumberMoved[,2]	#Move things around
  	
  	MovedAdults<- NatAge
  	
  }
  
  if (MoveType=='None')
  {
  	MovedAdults<- NatAge
  }
  
  return(MovedAdults)
}


Recruits<- function(NofEggs,DDType,RecDevForm,RecDevValue) #Calculate recruitment
{
  #This function calculates recruitment incorporating density dependence of choice. So far only works with 'BH' (Beverton-Holt)
  #Eggs= number of eggs per patch
  #DDType= Density dependence type, so far only 'BH' will work
  #RecDevForm= the type of recruitment deviation, can either be random (from a distribution) or known (a fixed multiplier). 
  #RecDevValue= the values of the recruitment deviate. If form is random, this needs to be a vector, with the first value being the mean of the log normal distribution (probably leave at 0), and the second being the standard deviation 
  
  if (RecDevForm=='Random')
  {
    RecDeviate<- rlnorm(1,meanlog=RecDevValue[1],sdlog=RecDevValue[2])  #Calculate log normal recruitment deviate (not patch specific yet)
  }
  
  if (RecDevForm=='Known')
  {
    RecDeviate<- RecDevValue
  }
  
  if (DDType=='BH')
  {
    
    RelativePatchSizes<- Patches$PatchSizes/sum(Patches$PatchSizes)
    
     Alpha<- (4* lh$BH.Steepness*RelativePatchSizes*lh$R0)/(5* lh$BH.Steepness-1)
        
     Beta<- (RelativePatchSizes*lh$R0)*(1-lh$BH.Steepness)/(5* lh$BH.Steepness-1)
    
	# show(paste('Alpha is ',Alpha))

	# show(paste('Beta is ',Beta))

     # Alpha<- (RelativePatchSizes*lh$B0)*((1-lh$BH.Steepness)/(4*lh$BH.Steepness*(RelativePatchSizes*lh$R0)))
    
     # Beta<- (5*lh$BH.Steepness-1)/(4*lh$BH.Steepness*(RelativePatchSizes*lh$R0))
    

    
    if (lh$LarvalChoice==1)
    {
       Recs<- ((Alpha*NofEggs)/(Beta+NofEggs)) * RecDeviate 
      
      # Recs<- (NofEggs/(Alpha+Beta*NofEggs)) * RecDeviate 

    }
    if (lh$LarvalChoice==0)
    {
        Recs<- ((Alpha*NofEggs)/(Beta+NofEggs)) * RecDeviate * Patches$HabQuality 
    
      # Recs<- (NofEggs/(Alpha+Beta*NofEggs)) * RecDeviate * Patches$HabQuality 

    }
    Recs[is.na(Recs)]<- 0
  }
  if (DDType=='RASS-BH')
  {
    RelativePatchSizes<- Patches$PatchSizes/sum(Patches$PatchSizes)
    Alpha<- lh$RassAlpha
    Beta<- lh$RassBeta*RelativePatchSizes    
    if (lh$LarvalChoice==1)
    {
      Recs<- (NofEggs*Alpha)*(1/(1+Beta*NofEggs)) * RecDeviate 
    }
    if (lh$LarvalChoice==0)
    {
      Recs<- (NofEggs*Alpha)*(1/(1+Beta*NofEggs)) * RecDeviate* Patches$HabQuality 
    }
    Recs[is.na(Recs)]<- 0
  }
return(Recs)

}



FishingSelectivity<- function(Lengths,s50,s95,NumDraws)
{
  # Lengths<- LengthAtAge
  # s50=35
  # s95=45
  # NumDraws=1000
  
  LengthDist<- LengthWiError(Lengths,1000)
  
  GroupSelectivity<- 1/(1+exp(-log(19)*((LengthDist-s50)/(s95-s50))))	
  
  MeanSelectivity<- rowMeans(GroupSelectivity)	
  
  return(MeanSelectivity)	
}

DistFishingAtAge<- function(Effort,Selectivity)
{
  # Effort<- c(0.1,0)
  # Selectivity<- TargetSelectivity
  # Population<- PopArray[10,,]
  NofAges<- length(Selectivity)
  TempSelect<- rep(Selectivity,length(Effort))
  TempSelect<- matrix(TempSelect,nrow=NofAges,ncol=length(Effort),byrow=F)	
  TempEffort<- rep(Effort,NofAges)
  TempEffort<- matrix(TempEffort,nrow=NofAges,ncol=length(Effort),byrow=T)
  FishMortalityAtAge<- TempEffort * TempSelect
  return(FishMortalityAtAge)
}


AssignNTZ<- function(NTZSize,CurrentPopulation)
{
  # NTZSize<- 1
  #Patches<- BasePatches
   Patches$MPALocations<<- c(1,0)
  # CurrentPopulation<- TempPop
  RelativePatchSize<- Patches$PatchSizes/sum(Patches$PatchSizes)
  OldPatches<- RelativePatchSize
  Patches$PatchSizes[Patches$MPALocations==1]<<- NTZSize/sum(Patches$MPALocations==1)
  Patches$PatchSizes[Patches$MPALocations==0]<<- (1-NTZSize)/sum(Patches$MPALocations==0)
  SizeChange<- (Patches$PatchSizes/OldPatches)-1
  SizeChange[is.nan(SizeChange)]<- 0
  ChangeSign<- sign(SizeChange)
  PopChange<- matrix(rep(ChangeSign,dim(CurrentPopulation)[1]),dim(CurrentPopulation),byrow=T) * (CurrentPopulation[,SizeChange<=0]*abs(SizeChange[SizeChange<=0]))
  NewPopulation<- CurrentPopulation+PopChange
  return(NewPopulation)
}


GoFish<- function(FishingMortality,TotalSurvive,Population)
{
  FishingProp<- FishingMortality/(FishingMortality+lh$m)
  DeadFish<- Population*(1-TotalSurvive)
  DeadFish[,Patches$MPALocations==1]<- 0
  NumbersAtAgeCaught<- DeadFish*FishingProp
  WeightAtAgeCaught<- NumbersAtAgeCaught * WeightAtAge
  TotalCaught<- colSums(WeightAtAgeCaught)
  return(list(NatAge=NumbersAtAgeCaught,WatAge=WeightAtAgeCaught,Total=TotalCaught))
}

FindReferencePoint<- function(FtoOptim,Target,TargetValue)
{
  # FtoOptim<- 1.1*rep(FishToHalfBmsy$par,NumPatches)
  # Target<- 'Biomass'
  # TargetValue<- 0.5
   
  FTemp<- rep(exp(FtoOptim),NumPatches)
  It<-1
  if(lh$RecDevSTD>0)
  {
    It<- 10
  }
  FinalYield<- NULL
  FinalPopulation<- matrix(NA,nrow<- It,ncol=NumPatches)
  for (i in 1:It)
  {  
    TempPop<- GrowPopulation(UnfishedPopulation,FTemp,'EQ',0,NA)
	TempFinalYield<-  TempPop$Performance$Yields
    FinalYield[i]<- TempFinalYield[length(TempFinalYield)]
    FinalPopulation[i,]<- colSums(TempPop$FinalNumAtAge*WeightAtAge)  
  }
  
  FinalYield<- mean(FinalYield)
  FinalPopulation<- colMeans(FinalPopulation)
  if (Target=='FMSY')
  {
    return(-FinalYield)
  }
  
  if (Target=='BvBmsy')
  {
    
    if (sum(FinalPopulation)==0)
    {
      FinalPopulation<- colSums(UnfishedPopulation*WeightAtAge)*exp(FTemp)
    }
    Ratio<- FinalPopulation/lh$Bmsy
#     show(FinalPopulation)
#     show(FTemp)
    # show(sum((Ratio-TargetValue)^2))
    return(sum((Ratio-TargetValue)^2))
  }
}

MoveFleet<- function(BaseF,PercentClosed,FleetSpill,FleetModel)
{
  
  if (FleetModel==0 & FleetSpill==1)
  {
    FTemp<- BaseF*(1+PercentClosed)
    FTemp<- FTemp*as.numeric(Patches$MPALocations==0)
  }
  if (FleetModel==0 & FleetSpill==0)
  {
    FTemp<- BaseF*as.numeric(Patches$MPALocations==0)
  }
  
  return(FTemp)
}



FindOptimalMPASize<- function(MPASizes,FTemp,FleetSpill,StartPop) #Given a number of MPAs and a fishing pressure find optimal size of each MPA
{
#   Patches<- BasePatches
#   Patches$MPALocations<- c(1,0)
  NewPop<- AssignNTZ(MPASizes,StartPop)
  NewF<- MoveFleet(FTemp,MPASizes,FleetSpill,0)
  # NewPop<- ShiftMPABorders(StartPop)
  # show(NewF)
  TempYields<- GrowPopulation(NewPop,NewF,'EQ',0,"Eh")$Performance$Yields
  EQYields<- TempYields[length(TempYields)]
  # show(-EQYields)
  Patches<<- BasePatches
  return(-EQYields)	
}

LoanPayment<- function(principle,rate,nyears,ppyear)
{
#Calculate annual loan payment for a given principal, rate, and payment period, (and payments per year)
# r=rates[3]
# ppyear=4
# nyears=10
# ppyear=4
# principle=tankcost*test$minimum
	pr=rate/ppyear
	n=nyears*ppyear
	lp=(pr*(principle))/(1-(1+pr)^-n)
	return(lp)
}

MPAFunction<- function(OptVector,t,OptSize)
{
	# MPASize<- (OptVector[1]+OptVector[2]*t+(OptVector[3]*t)^2)/1000

	# MPASize<- OptVector[1]*(1+OptVector[2]*exp(-OptVector[3]*t))

	# MPASize<- OptVector[1]*(1+OptVector[2]*exp(-OptVector[3]*t))

	 # OptVector<- c(5,6)
	
    MPASize<-(OptSize/(1+exp(-log(19)*((t-OptVector[1])/(OptVector[2]-OptVector[1])))))
	
	 # 	MPASize<- pmin(MPASize,1)

	# MPASize<- pmax(MPASize,0)

	return(MPASize)
}

FindMPATrajectory<- function(OptVector,TimeFrame,FTemp,FleetSpill,StartPop,OptMode,BaseYields,OptSize)
{

#   OptVector<- TestOpt$par
#   FTemp= FScenarios[f]
#   TimeFrame=OptTime
#   FleetSpill=FleetSpill
#   StartPop=StartPop
#   OptSize=OptNTZSize$par
#   OptMode='Function'
#   BaseYields=BaseConditions$Yield[f]
#   
  Yields<- NULL

  PassPop<- StartPop
# 
#   Patches<- BasePatches
#   Patches$MPALocations<-  c(1,0)
	for (t in 1:TimeFrame)
	{

		if (OptMode=='Function')
		{
      MPASize<- MPAFunction(OptVector,t,OptSize)
		}
				
		NewPop<- AssignNTZ(MPASize,PassPop)
    NewF<- MoveFleet(FTemp,MPASize,FleetSpill,0)
		
		NewPop<- GrowPopulation(NewPop,NewF,1,0,'eh')

		Yields[t]<- NewPop$Performance$Yields

		PassPop<- NewPop$FinalNumAtAge
		
	}
	NPB<- Discount((Yields)-(BaseYields),Fleet$YieldDiscount,TimeFrame)$NPV

Patches<<- BasePatches

  return(-(NPB))

}

FindMaxInterestRate<- function(InterestRate,LoanTime,LoanAmount,Surplus)
{
	     
	     InterestRate<- exp(InterestRate)
	     
	      LoanPayments<- colSums(((1+Fleet$YieldDiscount)^-(0:(LoanTime-1))) %*% t(LoanPayment(LoanAmount,InterestRate,LoanTime,1)))
	      
	      LoanAdjustedBalance<- Surplus-LoanPayments
	      
	      return((LoanAdjustedBalance^2))
} 
