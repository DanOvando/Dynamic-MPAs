RunReserve<- function(r,RunMatrix,BasePatches,DefaultLifeHistory)
{
  
  lh<- DefaultLifeHistory
  
  Run<- RunMatrix[r,]
  
  for (d in 1:dim(Run)[2])
  {
    eval(parse(text=paste(colnames(Run)[d],'<- I(Run[d])',sep='')))
  }
  
  
  TotalStorage<- as.data.frame(matrix(NA,nrow= TimeToRun,ncol=16))
  
  colnames(TotalStorage)<- c('Species','m','f','Frate','FracNTZ','Year','Yield','Biomass','Numbers','SQYield','SQNumbers','SQBiomass','OptNTZYield','OptNTZNumbers','OptNTZBiomass','OptNTZ')
  
  StoreRun<- paste(BatchFolder,Species,sep='') #1 if you want to create a seeded folder to store results, 0 if you want it in the generic working folder
  
  Patches<- BasePatches
  
  LifeHistory<- LifeHistories[LifeHistories$CommName %in% Species,]
  
  for (l in 1:length(LifeVars))
  {
    
    Where1<- LifeColumns ==LifeVars[l]
    
    Where2<- names(lh)==LifeVars[l]
    
    lh[Where2]<- LifeHistory[Where1]
    
  }
  
  lh$LengthMa95<-  1.01* lh$LengthMa50
  
  lh$AgeMa95<-  1.01* lh$AgeMa50
  
  lh$MoveType<- '2D'
  
  lh$MovementArray<- movArray(NumPatches,((lh$Range)*NumPatches)/5,'Wrap')
  
  lh$Bmsy<- -999
  
  lh$LengthAtAge<- Length(1:lh$MaxAge) #Calculate length at age vector
  
  lh$WeightAtAge<- Weight(lh$LengthAtAge,lh$WeightForm) #Calculate weight at age vector
  
  lh$FecundityAtAge<- Fecundity(lh$WeightAtAge,'Weight') #Calculate Fecundity at age vector
  
  if (lh$NoFecundRelate==1)
  {
    lh$FecundityAtAge<- lh$WeightAtAge
  }
  
  lh$MaturityMode<- as.character(LifeHistory$MaturityMode)
  
  lh$MaturityAtAge<- Maturity(1:lh$MaxAge, lh$MaturityMode) # Calculate maturity at age vector
  
  if (Species=='Shark')
  {
    lh$R0<- 100*Patches$HabQuality
  }
  
  
  pdf(file=paste(FigureFolder,'Life History.pdf',sep=''),width=8,height=6)
  par(mfrow=c(2,2),mar=c(4,4,1,1),oma=c(4,6,2,6))
  plot(lh$LengthAtAge/10,ylab='Length(cm)',type='l',lwd=2,xlab=NA,xaxt='n')
  title(main=Species)
  plot(lh$WeightAtAge*2.2,ylab='Weight(lbs)',type='l',lwd=2,xlab=NA,xaxt='n')
  plot(lh$MaturityAtAge,ylab='Probability Mature',type='l',lwd=2,xlab='Age')
  plot(lh$FecundityAtAge,ylab='# of Eggs',type='l',lwd=2,xlab='Age')
  
  dev.off()
  
  show(Species)
  
  ####### SET UP INITIAL POPULATION ########
  #     Rprof(tmp <- tempfile(),line.profiling=T)
  
  EQPopulation<- GrowPopulation(1000,rep(0,NumPatches),'EQ',0,'EQ Run',Species=Species,lh=lh,Patches=Patches) #Run the population out to unfished equilibrium
  
  lh<- EQPopulation$lh
  #     Rprof() 
  #     summaryRprof(tmp)
  #     unlink(tmp)
  
  lh$CarryingCapacityWeight<- (colSums(EQPopulation$FinalNumAtAge*lh$WeightAtAge)) #Calculate carrying capacity in weight
  
  UnfishedPopulation<- EQPopulation$FinalNumAtAge #Unfished numbers at age
  
  ####### Calculate Reference Points ########
  
  Fmsy<- optimize(log(lh$m[1]),f=FindReferencePoint,Target='FMSY',TargetValue=NA,lower=-10,upper=4,Species=Species,lh=lh,UnfishedPopulation=UnfishedPopulation,Patches=Patches) #Find FMSY  
  
  Fmsy$par<- exp(Fmsy$minimum) 
  
  BmsyPopulation<- GrowPopulation(UnfishedPopulation,rep(Fmsy$par,NumPatches),'EQ',0,'Bmsy Run',Species=Species,lh=lh,Patches=Patches) #Bmsy Population
  
  lh$Nmsy<- (colSums(BmsyPopulation$FinalNumAtAge)) #Nmsy
  
  lh$Bmsy<- colSums(BmsyPopulation$FinalNumAtAge*lh$WeightAtAge) #Bmsy
  
  MsyFishing<- GrowPopulation(BmsyPopulation$FinalNumAtAge,rep(Fmsy$par,NumPatches),OptTime,0,'MSY Run',Species=Species,lh=lh,Patches=Patches) #Bmsy Population
  
  Fleet$MSY_NPV<- MsyFishing$Performance$DiscYields$NPV
  
  FatTarget<- optimize(log(2*Fmsy$par),f=FindReferencePoint,Target='BvBmsy',TargetValue=FLevel,lower=-10,upper=4,Species=Species,lh=lh,UnfishedPopulation=UnfishedPopulation,Patches=Patches) #Find F that results in target B/Bmsy
  
  FatTarget$par<- exp(FatTarget$minimum)
  
  BatTargetPopulation<- GrowPopulation(UnfishedPopulation, rep(FatTarget$par,NumPatches),'EQ',0,'BatTarget Run',Species=Species,lh=lh,Patches=Patches)
  
  ####### RUN MPA SIMULATIONS ########
  
  MPAs<- as.data.frame(matrix(NA,ncol=length(MPANames),nrow=OptTime+1))
  
  colnames(MPAs)<- MPANames
  
  #     TimeToRun<-dim(MPAs)[1]
  
  EvalTime<- OptTime #Time span to evaluate results on
  RunTime<- 'Custom'
  PropNames<- ReserveStrategy
  #          for (m in 1:dim(MPAs)[2])
  #          {
  #            PropNames[m]<-MPANames[m]      
  #          }
  
  FScenarios<- c(FatTarget$par) #load in fishing scenarios
  
  BaseConditions<- as.data.frame(matrix(NA,nrow= length(FScenarios),ncol=3))
  colnames(BaseConditions)<- c('Yield','Biomass','Numbers')
  
  OptimalConditions<- as.data.frame(matrix(NA,nrow= length(FScenarios),ncol=3))
  colnames(OptimalConditions)<- c('Yield','Biomass','Numbers')
  
  BasePop<- GrowPopulation(EQPopulation$FinalNumAtAge,FatTarget$par,'EQ',1,paste('FvFmsy is',round(FatTarget$par/Fmsy$par,2)),Species=Species,lh=lh,Patches=Patches)
  
  StartPop<- BasePop$FinalNumAtAge
  
  BaseConditions$Yield<- round(BasePop$Performance$Yields[length(BasePop$Performance$Yields)],2)
  
  BaseConditions$Biomass<- round(sum(lh$WeightAtAge %*% BasePop$FinalNumAtAge),2)
  
  BaseConditions$Numbers<- round(sum(BasePop$FinalNumAtAge),2)
  
  wtf<- seq(0,1,length.out=20)
  
  arg=ldply(lapply(wtf,FindOptimalMPASize,FTemp=FatTarget$par,StartPop=StartPop,FleetSpill=FleetSpill,lh=lh,BasePatches=BasePatches,Species=Species))
  
  mguess<- wtf[which(arg==min(arg))[1]]
  
  OptNTZSize<-   optim(mguess,f=FindOptimalMPASize,lower=0,upper=0.999,FTemp=FatTarget$par,StartPop=StartPop,FleetSpill=FleetSpill,BasePatches=BasePatches,Species=Species,lh=lh) #You need a better optimization here, gets really stuck with any kind of stochasticity
  
  OptNTZSize<-   optim(OptNTZSize$par,f=FindOptimalMPASize,lower=0,upper=0.999,FTemp=FatTarget$par,StartPop=StartPop,FleetSpill=FleetSpill,BasePatches=BasePatches,Species=Species,lh=lh) #You need a better optimization here, gets really stuck with any kind of stochasticity
  
  mcheck<- data.frame(wtf,-(arg))
  
  colnames(mcheck)<- c('MPASize','NPB')
  
  pdf(file=paste(FigureFolder,'MPA Check.pdf',sep=''))
  print(ggplot(data=mcheck,aes(MPASize,NPB))+geom_point()+geom_vline(xintercept=OptNTZSize$par))
  dev.off()
  
  Int=seq(0,1,length.out=5)   
  
  Flip=seq(.001,1,length.out=5)   
  
  ObjMat<- matrix(NA,nrow=length(Int)*length(Flip),ncol=3)
  
  l<- 0
  MPAs<- NULL
  
  if(ReserveStrategy=='OptNTZ')
  {
    
    for (i in 1:length(Int))
    {
      for (o in 1:length(Flip))
      {
        l=l+1
        show(l)
        
        ObjMat[l,3]<- -FindMPATrajectory(c(Int[i],Flip[o]),Mode='LinearPlus',EvalTime=OptTime,FTemp= FatTarget$par,TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,
                                         OptSize=OptNTZSize$par,Alpha=1,OptMode='Function',
                                         BaseYields=BaseConditions$Yield,GrowMode='Shrink',Species=Species,lh=lh,BasePatches=BasePatches)
        
        ObjMat[l,1:2]<- c(Int[i],Flip[o])
      }
    }
    
    ObjMat<- as.data.frame(ObjMat)
    colnames(ObjMat)<- c('Intercept','FlipYear','Obj')
    pdf(file=paste(FigureFolder,'GridSearch.pdf',sep=''))
    p <- ggplot(ObjMat, aes(x=Intercept,y=FlipYear))
    print(p + geom_tile(aes(fill=Obj)))
    dev.off()
    
    BestGuess<- ObjMat[ObjMat$Obj==max(ObjMat$Obj,na.rm=T) & is.na(ObjMat$Obj)==F,1:2]
    
    MPAs$OptNTZ<- c(0, MPAFunction(as.matrix(BestGuess),1:OptTime, OptNTZSize$par,'LinearPlus',OptTime))
    
  }
  
  OptNTZSize$par<- round(OptNTZSize$par,4)
  
  #   if (OptNTZSize$convergence>0)
  #   {
  #     warning(paste('Failed to converge on run'))
  #   }
  
  Patches<- AssignNTZ(OptNTZSize$par,ReservePosition,BasePatches)
  
  FTemp<- MoveFleet(FatTarget$par, OptNTZSize$par,FleetSpill,0)
  
  OptPop<- GrowPopulation(StartPop,FTemp,'EQ',0,paste('Opt MPA is',round(OptNTZSize$par,2),'when FvFmsy is',round(FatTarget$par/Fmsy$par,2)),Species=Species,lh=lh,Patches=Patches)
  
  OptimalConditions$Yield<-(OptPop$Performance$Yields[length(OptPop$Performance$Yields)])
  
  OptimalConditions$Biomass<-(sum(lh$WeightAtAge %*% OptPop$FinalNumAtAge))
  
  OptimalConditions$Numbers<-(sum(OptPop$FinalNumAtAge))
  
  MPAs$StatusQuo<- 0
  
  MPAs$EqNTZ<- c(0,rep(OptNTZSize$par,OptTime))
  
  MPAs$SNTZ<- c(0,seq(min(1,1.5*OptNTZSize$par),OptNTZSize$par,length.out=OptTime))
  
  MPAs$GNTZ<- c(0,seq(min(1,0.5*OptNTZSize$par),OptNTZSize$par,length.out=OptTime))
  
  #       MPAs$OptNTZ<- c(0, MPAFunction(OptMPAPath$par,1:OptTime, OptNTZSize$par,'LinearPlus',EvalTime))
  
  
  MPAs$CatchShareEqNTZ<- MPAs$EqNTZ
  
  
  MPA<- MPAs[[which(names(MPAs) %in% ReserveStrategy)]]
  
  #   pdf(file=paste(FigureFolder,'FvFmsy is',round(FatTarget$par/Fmsy$par,2),' MPAs.pdf'),family=Font,pointsize=12,width=6,height=4)
  #   par(mar=c(5.1,4.1,4.1,6.1),xpd=T)
  #   matplot((MPAs),type='l',col=rainbow(1.5*dim(MPAs)[2])[1:dim(MPAs)[2]],lty=1,xlab='Year',ylab='% NTZ',lwd=4,bty='n')
  #   legend('topright',inset=c(-.3,0),legend=MPANames,col=rainbow(1.5*dim(MPAs)[2])[1:dim(MPAs)[2]],lty=1,cex=.5)
  #   dev.off()
  
  TempPop<- StartPop
  
  Patches<- BasePatches
  
  for (y in 1:TimeToRun)
  {
    
    CurrentMPA<- MPA[y]
    
    if (y>length(MPA))
    {
      #             CurrentMPA<- MPAs[dim(MPAs)[2],m]
      CurrentMPA<-  OptNTZSize$par
      
    }
    
    Patches<- AssignNTZ(CurrentMPA,ReservePosition,BasePatches)
    
    FTemp<- FatTarget$par 
    
    FVec<- MoveFleet(FTemp,CurrentMPA,FleetSpill,0)
    
    if (ReserveStrategy=='CatchShareEqNTZ'& y>1)
    {
      #             FVec<- MoveFleet(FTemp,CurrentMPA,0,0)            
      FVec<- MoveFleet(Fmsy$par,CurrentMPA,0,0)            
      
    }
    
    PassPop<- GrowPopulation(TempPop,FVec,1,0,'TEST',Species=Species,lh=lh,Patches=Patches) #grow the new population
    
    TempPop<- PassPop$FinalNumAtAge
    
    TotalStorage[y,]<- data.frame(as.character(Species),as.character(ReserveStrategy),paste('f',100*FLevel,sep=''),FatTarget$par,CurrentMPA,y,PassPop$Performance$MeanYield,PassPop$Performance$MeanBiomass,PassPop$Performance$MeanNumbers, BaseConditions$Yield,BaseConditions$Numbers, 
                                  BaseConditions$Biomass, 
                                  OptimalConditions$Yield, OptimalConditions$Numbers, OptimalConditions$Biomass, OptNTZSize$par,stringsAsFactors=F)
  } #Close TimeToRun loop   
  return(TotalStorage)
} #Close function