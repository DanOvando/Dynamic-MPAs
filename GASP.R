
####### G.A.S.P. (General age structured population model... sort of) ######

#Dan Ovando
#4.29/13
#Summary: This suite of functions run a general age structured fishery model. Will prepare a short summary shortly to use the set. 

GrowPopulation<- function(InitialPopulation,FishingPressure,Time,MakePlots,GroupFigName)
{
  
  RelativePatchSizes <- 1
  
  
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
  
  if (is.na(lh$B0[1]))
  {
    PopArray[1,1,]<- lh$R0
    
    for (a in 2:(lh$MaxAge-1))
    {
      PopArray[1,a,]<- PopArray[1,a-1,]*(1-lh$m)
    }
    
    PopArray[1,lh$MaxAge,]<- PopArray[1,lh$MaxAge-1,]*((1-lh$m)/(1-(1-lh$m)))
    
    lh$B0<<- colSums(PopArray[1,,]*WeightAtAge)
    
    lh$SSB0<<- colSums(PopArray[1,,]*WeightAtAge*MaturityAtAge)
    
    lh$SSB0_R0<<- lh$SSB0/lh$R0
    
    lh$SSN0<<- colSums(PopArray[1,,]*MaturityAtAge)
    
  }
  
  SSB<- colSums(PopArray[1,,]*MaturityAtAge*WeightAtAge)
  
  if (grepl('Shark',Species))
  {
    SSB<- colSums(PopArray[1,,]*MaturityAtAge)
  }
  
  while(t<TempTime) #Loop as long as you want
  {    
    
    
    RelocatedEggs<- MoveEggs(NumberOfEggs,'Common') #Move eggs
    
    RecruitsPerPatch<- Recruits(RelocatedEggs,SSB,lh$DDForm,'Random',c(lh$RecDevMean,lh$RecDevSTD),Species,lh) #Calculate recruits
    
    SelectivityAtAge<- FishingSelectivity(LengthAtAge,Fleet$s50,Fleet$s95,1000)	#Calculate selectivity curve
    
    FishingAtAge<- DistFishingAtAge(FishingPressure[1,],SelectivityAtAge) #Calculate fishing mortality at age
    
    FishingAtSpace<- DistFleet(FishingAtAge,SelectivityAtAge,PopArray[t,,]) #Distribute fishing fleet accordings to biomass
    
    Survival<- exp(-(FishingAtSpace+lh$m)) #calculate survival from each age class to the next
    
    FishingYields<- GoFish(FishingAtSpace,Survival,PopArray[t,,])     
    
    StoreFishedWeight[t,]<- FishingYields$Total 
    
    # PopArray[t,,]<- MoveAdults(PopArray[t,,],lh$Range,lh$MoveType)
    
    PlusGroup<- PopArray[t,lh$MaxAge,]*Survival[lh$MaxAge,] #The numbers in the plus group that survive
    
    PopArray[t+1,2:lh$MaxAge,]<- PopArray[t,1:(lh$MaxAge-1),] * Survival[1:(lh$MaxAge-1),] #Grow to the next age class
    
    PopArray[t+1,lh$MaxAge,]<- PlusGroup #Add in the survivors of the plus group
    
    PopArray[t+1,,]<- MoveAdults(PopArray[t+1,,],lh$Range,lh$MoveType)
    
    NumberOfEggs<- Eggs(PopArray[t+1,,],FecundityAtAge,MaturityAtAge) #Produce eggs by patch    
    
    SSB<- colSums(PopArray[t,,]*MaturityAtAge*WeightAtAge)
    
    if (grepl('Shark',Species))
    {
      SSB<- colSums(PopArray[t,,]*MaturityAtAge)
    }
    
    
    PopArray[t+1,1,]<- RecruitsPerPatch #Add recruits to population 
    
    #     PopArray[PopArray<1]<- 0
    
    WeightArray[t+1,,]<- PopArray[t+1,,] * WeightAtAge
    
    PopChange<- sum(PopTolerance< abs((colSums(WeightArray[t+1,,])-colSums(WeightArray[t,,])))) #Measure change in 
    
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
  }
  if (Time==1)
  {
    FinalNumAtAge<- PopArray
  }
  
  
  if (LookAtLengths==1)
  {
    
    if (Time>1)
    {
      FinalLengthFrequency<- LengthFrequency(PopArray[dim(PopArray)[1],,],LengthAtAge)
    }
    if (Time==1)
    {
      FinalLengthFrequency<- LengthFrequency(PopArray,LengthAtAge)
    }
  }
  
  if (lh$CarryingCapacity[1] !=-999)
  {
    Performance$FreqCollapse<- mean(TempTotalPop/sum(lh$CarryingCapacity)<= CollapseThreshold)
  }
  
  
  if (MakePlots==1)
  {   
    FormatFigure(paste(GroupFigName,'Total Numbers (N) Over Time.pdf'))
    matplot(rowSums(PopTrajectory),type='l',lwd=4,xlab='Time',ylab='Total Numbers',col=terrain.colors(2*NumPatches)[1:NumPatches],bty='n')
    #     legend('topright',legend=paste('Patch',1:NumPatches),col=terrain.colors(2*NumPatches)[1:NumPatches],lty=1:NumPatches,lwd=4,bty='n')
    dev.off()  
    
    if (lh$Bmsy[1]!=-999)
    {
      FormatFigure(paste(GroupFigName,'BvBmsy Over Time.pdf'))
      matplot(rowSums(WeightTrajectory)/sum(lh$Bmsy),type='l',lwd=4,xlab='Time',ylab='B/Bmsy',col=terrain.colors(2*NumPatches)[1:NumPatches],bty='n')
      #       legend('topright',legend=paste('Patch',1:NumPatches),col=terrain.colors(2*NumPatches)[1:NumPatches],lty=1:NumPatches,lwd=4,bty='n')
      dev.off()
      
    }    
    
    if (LookAtLengths==1)
    {
      FormatFigure(paste(GroupFigName,'Final Length Frequency Histogram.pdf'))
      
      FlatLengthFreq<- as.data.frame(FinalLengthFrequency)
      
      colnames(FlatLengthFreq)<- paste(1:NumPatches,sep='')
      
      FlatLengthFreq<- stack(as.data.frame(FlatLengthFreq))
      
      colnames(FlatLengthFreq)<- c('Length','Patch')
      
      PatchNumbers<- 1:NumPatches
      
      FlatLengthFreq$MPA[FlatLengthFreq$Patch %in% PatchNumbers[Patches$MPALocations==0]]<- 'Fished Area'
      
      FlatLengthFreq$MPA[FlatLengthFreq$Patch %in% PatchNumbers[Patches$MPALocations==1]]<- 'MPA'
      
      print(histogram(~Length | MPA ,data=FlatLengthFreq,xlab='Length (mm)',col=terrain.colors(2*NumPatches),type='count'))
      dev.off()
      
    }
    
    FormatFigure(paste(GroupFigName,'Fishing Yields.pdf'))
    
    FlatFishingYields<- as.data.frame(StoreFishedWeight)
    
    YearVector<- rep(as.numeric(rownames(StoreFishedWeight)),NumPatches)
    
    colnames(FlatFishingYields)<- paste(1:NumPatches,sep='')
    
    FlatFishingYields<- stack((FlatFishingYields))
    
    colnames(FlatFishingYields)<- c('Yield','Patch')
    
    FlatFishingYields$Year<- YearVector
    
    PatchNumbers<- 1:NumPatches
    
    FlatFishingYields$MPA[FlatFishingYields$Patch %in% PatchNumbers[Patches$MPALocations==0]]<- 'Fished Area'
    
    FlatFishingYields$MPA[FlatFishingYields$Patch %in% PatchNumbers[Patches$MPALocations==1]]<- 'MPA'
    
    FlatFishingYields<- ddply(FlatFishingYields,c('Year','MPA'),summarize,Yields=sum(Yield))
    
    print(xyplot(Yields ~Year  | MPA,data=FlatFishingYields,type='l',lwd=4))
    dev.off()
    
    FormatFigure(paste(GroupFigName,'Biomass Over Time.pdf'))
    
    FlatBiomass<- as.data.frame(WeightTrajectory)
    
    YearVector<- rep(as.numeric(rownames(FlatBiomass)),NumPatches)
    
    colnames(FlatBiomass)<- paste(1:NumPatches,sep='')
    
    FlatBiomass<- stack((FlatBiomass))
    
    colnames(FlatBiomass)<- c('Biomass','Patch')
    
    FlatBiomass$Year<- YearVector
    
    PatchNumbers<- 1:NumPatches
    
    FlatBiomass$MPA[FlatBiomass$Patch %in% PatchNumbers[Patches$MPALocations==0]]<- 'Fished Area'
    
    FlatBiomass$MPA[FlatBiomass$Patch %in% PatchNumbers[Patches$MPALocations==1]]<- 'MPA'
    
    FlatBiomass<- ddply(FlatBiomass,c('Year','MPA'),summarize,Biomass=sum(Biomass))
    
    print(xyplot(Biomass ~ Year | MPA,data=FlatBiomass,type='l',lwd=4))
    
    dev.off()
    
    FormatFigure(paste(GroupFigName,' Patch Biomass.pdf'))
    
    FlatFinalBiomass<- as.data.frame(WeightTrajectory[dim(WeightTrajectory)[1],])
    
    colnames(FlatFinalBiomass)<- 'Biomass'
    
    PatchNumbers<- 1:NumPatches
    
    FlatFinalBiomass$Patch<- PatchNumbers
    
    FlatFinalBiomass$MPA[FlatFinalBiomass$Patch %in% PatchNumbers[Patches$MPALocations==0]]<- 'Fished Area'
    
    FlatFinalBiomass$MPA[FlatFinalBiomass$Patch %in% PatchNumbers[Patches$MPALocations==1]]<- 'MPA'
    
    FlatFinalBiomass$Dummy<- 1
    
    pal <- colorRampPalette(c("blue", "red"))
    
    print (levelplot((Biomass+rnorm(NumPatches,0,1)) ~ Patch*Dummy ,ylab="",data=FlatFinalBiomass,col.regions=pal(2*NumPatches),
                     panel=function(x,...)
                     {
                       panel.levelplot(x,...)
                       if (sum(Patches$MPALocations)>0)
                       {
                         panel.abline(v=c(min(PatchNumbers[Patches$MPALocations==1],na.rm=T),max(PatchNumbers[Patches$MPALocations==1],na.rm=T)))
                       }
                     }))
    
    #     print (levelplot(Biomass ~ Patch*Dummy ,ylab="",data=FlatFinalBiomass,col.regions=pal(2*NumPatches),
    #                      panel=function(x,...)
    #                      {
    #                        panel.levelplot(x,...)
    #                        if (sum(Patches$MPALocations)>0)
    #                        {
    #                          panel.abline(v=c(min(PatchNumbers[Patches$MPALocations==1],na.rm=T),max(PatchNumbers[Patches$MPALocations==1],na.rm=T)))
    #                        }
    #                      }))
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
  #   Time<- dim(DataVec)[1]-1 #Get rid of stupid Timeframe thing, why did you do that
  #   show(dim(DataVec))
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

# LengthWiError2<- function(Lengths,NumDraws) #Produce a distribution 
# {
#   
#   HasNumber<- as.numeric(NumDraws>0)
#   NumDraws[NumDraws==0]<- 1
#   LenSD<- lh$VBSD*(1+lh$VBErrorSlope*Lengths/lh$Linf)
# 
#   LengthAtAgeDist<- ldply(lapply(1:length(Lengths),LengthDist,HasNumber=HasNumber,Lengths=Lengths,NumDraws=NumDraws,lh=lh,LenSD=LenSD))
#   
#   return(LengthAtAgeDist)
# }
# 
# LengthDist<- function(l,HasNumber,Lengths,NumDraws,lh,LenSD)
# {
#   LengthAtAgeDist<- HasNumber*(Lengths[l]* rlnorm(NumDraws,mean=lh$VBErrorMean,sd=LenSD[l]))
#   return(LengthAtAgeDist)
# }

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
  LengthDist<- LengthWiError(Lengths,100)
  
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
    LengthDist<- LengthWiError(Data,100)
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
    
    LengthDist<- LengthWiError(Data,100)
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
  
  #   NatAge<- FinalNumAtAge
  #   MoveType<- '2D'
  #   Distance<- 0.9
  #     
  
  if (Distance==0)
  {
    MovedAdults<- NatAge
  }
  
  if (MoveType=='2D')
  {
    
    MovedAdults<- t(lh$MovementArray %*% t(NatAge))
    
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


Recruits<- function(NofEggs,SSB,DDType,RecDevForm,RecDevValue,Species,lh) #Calculate recruitment
{
  #This function calculates recruitment incorporating density dependence of choice. So far only works with 'BH' (Beverton-Holt)
  #Eggs= number of eggs per patch
  #DDType= Density dependence type, so far only 'BH' will work
  #RecDevForm= the type of recruitment deviation, can either be random (from a distribution) or known (a fixed multiplier). 
  #RecDevValue= the values of the recruitment deviate. If form is random, this needs to be a vector, with the first value being the mean of the log normal distribution (probably leave at 0), and the second being the standard deviation 
  h<- lh$BH.Steepness
  if (grepl('Shark',DDType))
  {
    #     browser()
    #     Recs<- (lh$R0*SSB*5)/(lh$SSN0*5+(2.26-1)*SSB*5)
    Recs<- (0.8*lh$R0*h*SSB)/(0.2*lh$SSN0*(1-h)+(h-0.2)*SSB)
    
  }
  
  if (grepl('Babcock',DDType))
  {
    if (DDType=='Babcock1')
    {
      Recs<- (0.8*lh$R0*h*SSB)/(0.2*lh$SSB0*(1-h)+(h-0.2)*SSB)
    }
    
    if (DDType=='Babcock2')
    {
      Recs<- rep((1/NumPatches)*((0.8*sum(lh$R0)*h*sum(SSB))/(0.2*sum(lh$SSB0)*(1-h)+(h-0.2)*sum(SSB))),NumPatches)
    }
    if (DDType=='Babcock3')
    {
      Recs<- rep((1/NumPatches)*sum((0.8*lh$R0*h*SSB)/(0.2*lh$SSB0*(1-h)+(h-0.2)*SSB)),NumPatches)
    }
    if (DDType=='Babcock4')
    {
      Recs<- ((0.8*lh$R0*h)/(0.2*lh$SSB0*(1-h)+(h-0.2)*SSB))*(1/NumPatches)*sum(SSB)
      
    }
    if (DDType=='Babcock5')
    {
      Recs<- MoveEggs((0.8*lh$R0*h*SSB)/(0.2*lh$SSB0*(1-h)+(h-0.2)*SSB),'Common')
      
    }
    
  }
  
  if (Species=='Generic')
  {
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
      
      #     RelativePatchSizes<- Patches$PatchSizes/sum(Patches$PatchSizes)
      
      RelativePatchSizes<- 1
      
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
  }
  #   if (Species=='Yellowtail Snapper')
  #   {
  #     
  #     alpha<- ((4*lh$BH.Steepness*lh$SSB0)/lh$SSB0_R0)/(5*lh$BH.Steepness-1)
  #     
  #     beta<- (lh$SSB0*(1-lh$BH.Steepness))/(5*lh$BH.Steepness-1)
  #     
  #     Recs<- (alpha*SSB)/(beta+SSB)
  #     
  #   }
  #   
  
  
  return(Recs)
  
}



FishingSelectivity<- function(Lengths,s50,s95,NumDraws)
{
  # Lengths<- LengthAtAge
  # s50=35
  # s95=45
  # NumDraws=1000
  
    LengthDist<- LengthWiError(Lengths,100)
    
    GroupSelectivity<- 1/(1+exp(-log(19)*((LengthDist-s50)/(s95-s50))))	
    
    MeanSelectivity<- rowMeans(GroupSelectivity)	
    
  #   LengthDist<- LengthWiError(Lengths,100)
  
#   GroupSelectivity<- 1/(1+exp(-log(19)*((Lengths-s50)/(s95-s50))))  
#   
#   MeanSelectivity<- GroupSelectivity
  
  return(MeanSelectivity)	
}

DistFishingAtAge<- function(Effort,Selectivity)
{
  
  #   Effort<- 0.35
  #   Selectivity<- SelectivityAtAge
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

AssignNTZ<- function(NTZSize,Position)
{
  # AssignNTZ<- function(NTZSize,CurrentPopulation)
  # NTZSize<- 1
  Patches<<- BasePatches
  
  if (Position=='Edge')
  {
    NTZBorder<- min(NumPatches,max(0,round(NumPatches*NTZSize)))
    
    Patches$MPALocations[0:NTZBorder]<<- 1
  }
  if (Position=='Center')
  {
    
    CenterPatch<- (NumPatches/2)
    
    MPASize<- NTZSize*NumPatches
    
    Bounds<- round((CenterPatch+1-MPASize/2)):round(CenterPatch+MPASize/2)
    
    Patches$MPALocations[Bounds]<<- 1
    
  }
  #   # CurrentPopulation<- TempPop
  #   RelativePatchSize<- Patches$PatchSizes/sum(Patches$PatchSizes)
  #   OldPatches<- RelativePatchSize
  #   Patches$PatchSizes[Patches$MPALocations==1]<<- NTZSize/sum(Patches$MPALocations==1)
  #   Patches$PatchSizes[Patches$MPALocations==0]<<- (1-NTZSize)/sum(Patches$MPALocations==0)
  #   SizeChange<- (Patches$PatchSizes/OldPatches)-1
  #   SizeChange[is.nan(SizeChange)]<- 0
  #   ChangeSign<- sign(SizeChange)
  #   PopChange<- matrix(rep(ChangeSign,dim(CurrentPopulation)[1]),dim(CurrentPopulation),byrow=T) * (CurrentPopulation[,SizeChange<=0]*abs(SizeChange[SizeChange<=0]))
  #   NewPopulation<- CurrentPopulation+PopChange
  #   return(NewPopulation)
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
  #   NewPop<- AssignNTZ(MPASizes,StartPop)
  
  AssignNTZ(MPASizes,ReservePosition)
  #   show(Patches$MPALocation)
  NewF<- MoveFleet(FTemp,MPASizes,FleetSpill,0)
  # NewPop<- ShiftMPABorders(StartPop)
  #    show(NewF)
  TempYields<- GrowPopulation(StartPop,NewF,'EQ',0,"Eh")$Performance$Yields
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

# MPAFunction<- function(OptVector,t,OptSize,Mode,EvalTime,GrowMode)
# {
#   # MPASize<- (OptVector[1]+OptVector[2]*t+(OptVector[3]*t)^2)/1000
#   
#   # MPASize<- OptVector[1]*(1+OptVector[2]*exp(-OptVector[3]*t))
#   
#   # MPASize<- OptVector[1]*(1+OptVector[2]*exp(-OptVector[3]*t))
#   
#   # OptVector<- c(5,6)
#   
#   
#   if (Mode=='FreeLogistic')
#   {
#     
#     #     MPASize<- OptVector[3]*(1/(1+exp(-log(19)*((t-OptVector[1])/(OptVector[2]-OptVector[1])))))
#     
#     #     MPASize<- OptVector[3]*(1/(1+exp(-log(19)*((t-OptVector[1])/((OptVector[1]+OptVector[2])-OptVector[1])))))
#     
#     MPAFactor<- (1/(1+exp(-log(19)*((t-OptVector[1])/((OptVector[1]+OptVector[2])-OptVector[1])))))
#     
#     if (GrowMode=='Shrink')
#     {
#       MPAFactor<- MPAFactor+1
#     }
#     
#     MPASize<- min(OptSize*MPAFactor,1)
#     
#     #     MPASize<- min(OptSize[1]*t+OptSize[2],OptSize)a
#     
#     #     show(OptVector[2]>OptVector[1])
#     
#   }
#   if (Mode=='LockedLogistic')
#   {
#     
#     #     OptVector<- c(-2,EvalTime-1)
#     MPASize<- (1/(1+exp(-log(19)*((t-OptVector[1])/(OptVector[2]-OptVector[1])))))
#     
#   }
#   if (Mode=='Linear')
#   {
#     
#     Slope<- (OptSize-OptVector[1])/EvalTime
#     
#     MPASize<- Slope*t+OptVector[1]
#     
#   }
#   if (Mode=='Logit Grow')
#   {
#     MPASize<- pmax(0,pmin(1,(exp(OptVector[1]*t)/(OptVector[2]+exp(OptVector[1]*t)) - 1/(1+OptVector[2]))+OptVector[3]))
#   }
#   if (Mode=='Logit Shrink')
#   {
#     
#     MPASize<- 1-pmax(0,pmin(1,(exp(OptVector[1]*t)/(OptVector[2]+exp(OptVector[1]*t)) - 1/(1+OptVector[2]))+OptVector[3]))
#     
#   }
#   if (Mode=='LinearPlus')
#   {
#     #     Intercept=exp(OptVector[1])/(1+exp(OptVector[1]))
#     Intercept=OptVector[1]
#     
#     FlipYear<- max(1,min(round(EvalTime*OptVector[2]),EvalTime))
#     
#     Slope<- (OptSize-Intercept)/FlipYear
#     
#     MPASize<- Slope*t+Intercept
#     
#     MPASize<- pmax(0,pmin(1,MPASize))
#     
#     if (length(t)==1)
#     {
#       #       if(t>=OptVector[2]){MPASize<- OptSize}
#       if(t>=FlipYear){MPASize<- OptSize}
#       
#     }
#     if (length(t)>1)
#     {
#       MPASize[t>FlipYear]<- OptSize
#     }
#   }
#   
#   
#   # 	MPASize<- pmin(MPASize,1)
#   
#   # MPASize<- pmax(MPASize,0)
#   
#   return(MPASize)
# }

MPAFunction<- function(OptVector,t,OptSize,Mode,EvalTime,GrowMode)
{
  
  
  Intercept=OptVector[1]
  
  FlipYear<- max(1,min(round(EvalTime*OptVector[2]),EvalTime))
  
  Slope<- (OptSize-Intercept)/FlipYear
  
  MPASize<- Slope*t+Intercept
  
  MPASize<- pmax(0,pmin(1,MPASize))
  
  if (length(t)==1)
  {
    #       if(t>=OptVector[2]){MPASize<- OptSize}
    if(t>=FlipYear){MPASize<- OptSize}
    
  }
  if (length(t)>1)
  {
    MPASize[t>FlipYear]<- OptSize
  }
  
  
  return(MPASize)
}
FindMPATrajectory<- function(OptVector,TimeFrame,FTemp,FleetSpill,StartPop,OptMode,BaseYields,OptSize,Mode,EvalTime,Alpha,GrowMode)
{
  
  Yields<- rep(NA,TimeFrame)
  
  PassPop<- StartPop
  
  #   if (OptVector[1]!=OptVector[2])    
  #   {
#   OptVector<- exp(OptVector)
  for (t in 1:TimeFrame)
  {
    
    if (OptMode=='Function')
    {
      MPASize<- MPAFunction(OptVector,t,OptSize,Mode,EvalTime,GrowMode)
    }
    
    AssignNTZ(MPASize,ReservePosition)
    
    NewF<- MoveFleet(FTemp,MPASize,FleetSpill,0)
    
    NewPop<- GrowPopulation(PassPop,NewF,1,0,'eh')
    
    Yields[t]<- NewPop$Performance$Yields
    
    
    PassPop<- NewPop$FinalNumAtAge
    
  }
#         NPB<- Discount((Yields)-(BaseYields),Fleet$YieldDiscount,TimeFrame)$NPV  
  
  NPB<- sum((Yields - BaseYields)*(1+Fleet$YieldDiscount)^-(0:(TimeFrame-1)))
  
  #     show(Discount((Yields),Fleet$YieldDiscount,TimeFrame)$NPV)
  
  #     RelativeNPV<- Yields[length(Yields)]
  #     RelativeNPV<- Discount((Yields),Fleet$YieldDiscount,TimeFrame)$NPV
  #     RelativeNPV<- Discount((Yields),Fleet$YieldDiscount,TimeFrame)$NPV / Fleet$MSY_NPV
  
  #     RelativeNPV<- Discount((Yields)-(BaseYields),Fleet$YieldDiscount,TimeFrame)$NPV  
  
  #     RelativeNPV<- Discount((Yields)-(BaseYields),0.05,TimeFrame)$NPV  
  
  
  Objective<- NPB
  
#   if (OptVector[1]>1 | OptVector[2]> 1 )
#   {
#     Objective<- -Inf
#   }
  
  #     Objective<- Alpha*RelativeNPV + (1-Alpha)*(sum(colSums(PassPop*WeightAtAge))/sum(lh$CarryingCapacityWeight))
  
  #Change this to be alpha*NPB+ (1-alpha)*(NP)
  
  Patches<<- BasePatches
  
  #   } #Close if statement on starting objecive
  
  return(-(Objective))
  
}


FindMaxInterestRate<- function(InterestRate,LoanTime,LoanAmount,Surplus)
{
  
  InterestRate<- exp(InterestRate)
  
  LoanPayments<- colSums(((1+Fleet$YieldDiscount)^-(0:(LoanTime-1))) %*% t(LoanPayment(LoanAmount,InterestRate,LoanTime,1)))
  
  LoanAdjustedBalance<- Surplus-LoanPayments
  
  return((LoanAdjustedBalance^2))
} 

movArray<-function(SpaceC,sdy,Form)
{
  
  if (Form=='Wrap')
  {
    
    
    #     SpaceC<- NumPatches
    #      sdy<- lh$Range*NumPatches
    #     
    P<- SpaceC
    sigmaL<- sdy  
    ##################################
    #Movement
    ##################################
    
    #dispersal:common larval pool
    DL1=matrix(1/(P),nrow=P,ncol=P)
    
    #dispersal:gaussian movement
    #create movement probability matrix, start with distance matrix
    loc<-seq(-(P-1),2*P,1)
    area<-seq(1,P,1)
    
    #create distance matrix
    dist<-matrix(NA,nrow=P,ncol=P*3)
    for(i in 1:P)
    {
      for(j in 1:(P*3))
      {
        dist[i,j]<-area[i]-loc[j]
      }
    }
    #now create the movement matrix of probabilities of movement
    p.init<-round(exp(-((dist)^2)/(2*(sigmaL^2))),2)
    
    #now add matrices on ends to wrap movement and normalize so movement from any one area sums to one
    p.all<-matrix(NA,nrow=P,ncol=(3*P))
    for(i in 1:P)
    {
      for(j in 1:(3*P))
      {
        p.all[i,j]<-(p.init[i,j])/sum(p.init[i,])
      }
    }
    
    p1<-p.all[,1:P]
    p2<-p.all[,(2*P+1):(3*P)]
    parea<-p.all[,(P+1):(2*P)]
    SpaceIn<-p1+p2+parea
    
    FormatFigure('Movement Probabilities.pdf')
    contour(SpaceIn)
    dev.off()
    
  } #Close Form Loop
  if (Form=='Bounce')
  {
    
    SpaceR<- 1
    #      SpaceC<- NumPatches
    sdx<- 0.9
    #      sdy<- 0.9
    
    SpaceIn<-array(dim=c(SpaceR,SpaceC,SpaceC*SpaceR))
    # bivariate normal diffusion
    # list the coordinates 
    coords<-NULL
    for(x in 1:SpaceR)
    {
      temp<-cbind(rep(x,SpaceC),seq(1,SpaceC))
      coords<-rbind(coords,temp)
    }
    
    # population the probability array of moving 
    for(h in 1:nrow(coords))
      for(j in 1:SpaceR)
        for(k in 1:SpaceC)
          #         SpaceIn[j,k,h]<- exp(-(((coords[h,1]-j))^2/(2*sdx^2)+((coords[h,2]-k))^2/(2*sdy^2)))
          SpaceIn[j,k,h]<- exp(-(((coords[h,1]-j))^2/(2*sdx^2)+((coords[h,2]-k))^2/(2*sdy^2)))
    
    
    # normalize so that it sums to 1 (effectively wraps around)
    for(h in 1:nrow(coords))
      SpaceIn[,,h]<-SpaceIn[,,h]/sum(SpaceIn[,,h])
    
    dim(SpaceIn)<- c(NumPatches,NumPatches)
    
    FormatFigure('Movement Probabilities.pdf')
    contour(SpaceIn)
    dev.off()
  }
  
  return(SpaceIn)
}

DistFleet<- function(FishingAtAge,SelectivityAtAge,NumAtAge)
{
  ### Distribute fishing effort by proportional biomass in fishable patches
  
  #      NumAtAge<- PopArray[1,,]
  
  WeightsAtAge<- NumAtAge *WeightAtAge
  
  CommercialBiomass<- (WeightsAtAge * (SelectivityAtAge))
  
  ProportionalBiomass<- colSums(CommercialBiomass)* as.numeric(Patches$MPALocations==0)
  
  #   Patches$MPALocations[1]<- 1
  
  if (sum(Patches$MPALocations)>0)
  {
    
    FishingAtSpace<- t(t(FishingAtAge) *  (ProportionalBiomass/max(ProportionalBiomass)))
  }
  else {FishingAtSpace<- FishingAtAge}
  
  FishingAtSpace[is.na(FishingAtSpace)]<- 0
  return(FishingAtSpace)
}
