rm(list=ls())
library(lattice)
require(stats)

####### Shrinking MPAs #########
#Wilson, Doughtery, Ovando et al something something

####### READ IN CONTROL FILE ########

setwd('/Users/danovando/Dropbox/Shrinking NTZ')

RunAnalysis<- 1

# SpeciesList<- c('YellowtailSnapper','RedGrouper') #c('Sheephead','Grouper','BlackPerch','OceanWhitefish','Opaleye','RedSeaUrchin','KelpBass', 'KelpRockfish','CaliforniaHalibut')

DataNames<- c('NPV of Yield','NPV of Biomass','Mean Yield','Mean Biomass','Mean Numbers','Yield Instability','Mean Changes in Yield','Percent Years with Profit Gains','Percent Years with Numbers Gains','Mean Percent Change in Yield','Mean Percent Change in Numbers','Percent Years With Numbers and Yield Gains','FiveyearYieldBalance','FiveyearBiomassBalance','FiveyearNPVBalance','TenyearYieldBalance','TenyearBiomassBalance','TenyearNPVBalance','YearsToYieldRecovery','YearsToBioRecovery','YearsToBalanceRecovery','TenYearNPSurplus','RequestedLoan','MaxInterestRate')


LongDataNames<- c('Species','Movement','Steepness','MPAScenario','FishingScenario','FvFmsy','OptNTZ','NPV of Yield','NPV of Biomass','Mean Yield','Mean Biomass','Mean Numbers','Yield Instability','Mean Changes in Yield','Percent Years with Profit Gains','Percent Years with Numbers Gains','Mean Percent Change in Yield','Mean Percent Change in Numbers','Percent Years With Numbers and Yield Gains','FiveyearYieldBalance','FiveyearBiomassBalance','FiveyearNPVBalance','TenyearYieldBalance','TenyearBiomassBalance','TenyearNPVBalance','YearsToYieldRecovery','YearsToBioRecovery','YearsToBalanceRecovery','TenYearNPSurplus','RequestedLoan','MaxInterestRate')

AllSpeciesStorage<- as.data.frame(matrix(NA,nrow=0,ncol=18))
colnames(AllSpeciesStorage)<- c('Species','m','f','Frate','FracNTZ','Year','Yield','Biomass','Numbers','SQYield','SQNumbers','SQBiomass','OptNTZYield','OptNTZNumbers','OptNTZBiomass','OptNTZ','YieldBalance','NPB')


AllSpeciesExperimentResults<- as.data.frame(matrix(NA,nrow=0,ncol=length(LongDataNames)))

colnames(AllSpeciesExperimentResults)<- c(LongDataNames)

BatchFolder<- 'Results/Scratch/'
dir.create(paste(BatchFolder,sep=''))

MPANames<- c('Status Quo','EqNTZ','Rotate','SNTZ','Basic','GNTZ','OptNTZ')



if (RunAnalysis==1)
{

for (s in 1:length(SpeciesList))
{
  Species<- SpeciesList[s] #species file to load
  
  StoreRun<- paste(BatchFolder,Species,sep='') #1 if you want to create a seeded folder to store results, 0 if you want it in the generic working folder
  source('SHRINKNTZ_CONTROLFILE.R')
  # Fleet$YieldDiscount<- 0
  
  lh$MoveType<- 'Simple-DO'
  # lh$MoveType<- 'None'
  show(SpeciesList[s])
  
  ####### SET UP INITIAL POPULATION ########
  
  EQPopulation<- GrowPopulation(100,rep(0,NumPatches),'EQ',1,'EQ Run')
  
  lh$CarryingCapacity<- (colSums(EQPopulation$FinalNumAtAge))
  
  lh$CarryingCapacityWeight<- (colSums(EQPopulation$FinalNumAtAge*WeightAtAge))
  
  UnfishedPopulation<- EQPopulation$FinalNumAtAge
  
  ####### Calculate Reference Points ########
  
  
  Fmsy<- optim(lh$m[1],FindReferencePoint,lower=0,upper=1,Target='FMSY',method="L-BFGS-B") #You need a better optimization here, gets really stuck with any kind of stochasticity
  
  BmsyPopulation<- GrowPopulation(UnfishedPopulation,rep(Fmsy$par,NumPatches),'EQ',1,'Bmsy Run')
  
  lh$Nmsy<- (colSums(BmsyPopulation$FinalNumAtAge))
  
  lh$Bmsy<- colSums(BmsyPopulation$FinalNumAtAge*WeightAtAge)
  
  FmsyPop<- GrowPopulation(UnfishedPopulation,Fmsy$par,'EQ',1,'FMSY RUN') #Fish down to MSY

   DIE<- GrowPopulation(UnfishedPopulation,4*Fmsy$par,100,1,'DIE')

  
  # OverFish<- GrowPopulation(UnfishedPopulation,10*rep(Fmsy$par,NumPatches),'EQ',1,'Overfished Run')
  
  ####### RUN MPA SIMULATIONS ########
  
  # Fmsy5xPop<- GrowPopulation(UnfishedPopulation,6*Fmsy$par,'EQ',1,'100x FMSY RUN') #Fish at X times MSY
  
  MPAs<- read.csv(paste(InputFolder,'SNTZ.csv',sep='')) #read in MPAs
  TimeToRun<-dim(MPAs)[2]
  # TimeToRun<- 600
  MPAs<- MPAs[1:7,]
  
  EvalTime<- 11 #Time span to evaluate results on
  RunTime<- 'Custom'
  PropNames<- NULL
  for (m in 1:dim(MPAs)[1])
  {
    PropNames[m]<- paste(round(MPAs[m,2:TimeToRun],2),collapse=':')
  }
  
  
  BasePatches<- Patches
  
  # FScenarios<- c(Fmsy$par,2*Fmsy$par,4*Fmsy$par) #load in fishing scenarios
  
  # FScenarios<- c(1.5*Fmsy$par,3*Fmsy$par,4.5*Fmsy$par) #load in fishing scenarios

  FScenarios<- c(1.5*Fmsy$par,3*Fmsy$par) #load in fishing scenarios
  
  
  ExperimentResults<-(array(NA,dim=c(dim(MPAs)[1],length(FScenarios),length(DataNames)))) #blank
  
  # dimnames(ExperimentResults)<- list(PropNames,c('Fmsy','4xFmsy','6xFmsy'),DataNames )
  
  dimnames(ExperimentResults)<- list(PropNames,c('1_5xFmsy','3xFmsy'),DataNames )
  
  
  FlatExperimentResults<- as.data.frame(matrix(NA,nrow=length(FScenarios)*dim(MPAs)[1],ncol=length(LongDataNames)))
  
  colnames(FlatExperimentResults)<- c(LongDataNames)
  
  BaseConditions<- as.data.frame(matrix(NA,nrow= length(FScenarios),ncol=3))
  colnames(BaseConditions)<- c('Yield','Biomass','Numbers')
  
  OptimalConditions<- as.data.frame(matrix(NA,nrow= length(FScenarios),ncol=3))
  colnames(OptimalConditions)<- c('Yield','Biomass','Numbers')
  
  if (RunTime=='Custom')
  {
    TimeToRun<- 75
  }
  
  
  BaseRelativeNumbers<- array(NA,dim=c(dim(MPAs)[1],length(FScenarios),TimeToRun))
  BaseRelativeYield<- array(NA,dim=c(dim(MPAs)[1],length(FScenarios),TimeToRun))
  BaseRelativeBiomass<- array(NA,dim=c(dim(MPAs)[1],length(FScenarios),TimeToRun))
  
  RawNumbers<- array(NA,dim=c(dim(MPAs)[1],length(FScenarios),TimeToRun))
  RawYield<- array(NA,dim=c(dim(MPAs)[1],length(FScenarios),TimeToRun))
  RawBiomass<- array(NA,dim=c(dim(MPAs)[1],length(FScenarios),TimeToRun))
  
  OptRelativeNumbers<- array(NA,dim=c(dim(MPAs)[1],length(FScenarios),TimeToRun))
  OptRelativeYield<- array(NA,dim=c(dim(MPAs)[1],length(FScenarios),TimeToRun))
  OptRelativeBiomass<- array(NA,dim=c(dim(MPAs)[1],length(FScenarios),TimeToRun))
  
  
  FleetSpill<- 1
  
  TotalStorage<- as.data.frame(matrix(NA,nrow<- length(FScenarios)*dim(MPAs)[1]*TimeToRun,ncol=18))
  
  colnames(TotalStorage)<- c('Species','m','f','Frate','FracNTZ','Year','Yield','Biomass','Numbers','SQYield','SQNumbers','SQBiomass','OptNTZYield','OptNTZNumbers','OptNTZBiomass','OptNTZ','YieldBalance','NPB')
  
  c<-0
  
  cc<- 0
  
  for (f in 1:length(FScenarios)) 
  {
    show(paste('F is ',f))
    Patches<- BasePatches
    BasePop<- GrowPopulation(EQPopulation$FinalNumAtAge,FScenarios[f],'EQ',1,paste('FvFmsy is',round(FScenarios[f]/Fmsy$par,2)))
    Patches$MPALocations<- c(1,0) #Create an MPA
    StartPop<- BasePop$FinalNumAtAge
    BaseConditions$Yield[f]<- round(BasePop$Performance$Yields[length(BasePop$Performance$Yields)],2)
    BaseConditions$Biomass[f]<- round(sum(WeightAtAge %*% BasePop$FinalNumAtAge),2)
    BaseConditions$Numbers[f]<- round(sum(BasePop$FinalNumAtAge),2)
    
    OptNTZSize<- 	optim(0.05,FindOptimalMPASize,lower=0,upper=0.999,FTemp=FScenarios[f],StartPop=StartPop,FleetSpill=FleetSpill,method="Brent") #You need a better optimization here, gets really stuck with any kind of stochasticity
    
   OptTime<- 10
   TestOpt<- optim(c(10,-1),FindMPATrajectory,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f])    

# # SizeVector<- -50:50
# # ResultMat<- matrix(nrow=length(SizeVector),ncol=length(SizeVector)) 

# MPAFunction(TestOpt$par,1:10, OptNTZSize$par)


# # for (ss in 1:length(SizeVector))
# {
	# for (sss in 1:length(SizeVector))
	# {
		
		# ResultMat[ss,sss]<- FindMPATrajectory(c(SizeVector[ss],(SizeVector+.0001)[sss]),OptTime, FScenarios[f],FleetSpill,StartPop,'Function',BaseConditions$Yield[f], OptNTZSize$par)
		
		# MPAFunction(c(SizeVector[ss],(SizeVector)[sss]),1:OptTime, OptNTZSize$par)
	# }
# }

# quartz()
# image(SizeVector,SizeVector,ResultMat)


     OptNTZSize$par<- round(OptNTZSize$par,4)
    if (OptNTZSize$convergence>0)
    {
      warning(paste('Failed to converge on run ',f))
    }
    
    Patches<- BasePatches  
    
    Patches$MPALocations<- c(1,0) #Create an MPA
    
    TempPop<- AssignNTZ(OptNTZSize$par,StartPop)
    
    FTemp<- MoveFleet(FScenarios[f], OptNTZSize$par,FleetSpill,0)
    
    OptPop<- GrowPopulation(TempPop,FTemp,'EQ',1,paste('Opt MPA is',round(OptNTZSize$par,2),'when FvFmsy is',round(FScenarios[f]/Fmsy$par,2)))
    
    OptimalConditions$Yield[f]<-(OptPop$Performance$Yields[length(OptPop$Performance$Yields)])
    
    OptimalConditions$Biomass[f]<-(sum(WeightAtAge %*% OptPop$FinalNumAtAge))
    
    OptimalConditions$Numbers[f]<-(sum(OptPop$FinalNumAtAge))
    
    
    MPAs[2,1]<- 0
    
    MPAs[2,2:dim(MPAs)[2]]<- OptNTZSize$par
    # MPAs[3,]<- c(0,1,1,OptNTZSize$par,OptNTZSize$par,OptNTZSize$par)
    
    MPAs[3,]<- c(0,0.25,0.25,0,0,0.25,0.25,0,0,0.25,0.25)
    
    Rotations<- rep(c(0,0,.25,.25),40)
    
    MPAs[4,]<- c(0,seq(min(1,1.5*OptNTZSize$par),OptNTZSize$par,length.out=10))
    
    MPAs[6,]<- c(0,seq(min(1,0.5*OptNTZSize$par),OptNTZSize$par,length.out=10))
    
    MPAs[7,]<- c(0, MPAFunction(TestOpt$par,1:10, OptNTZSize$par))
    
    
    pdf(file=paste(FigureFolder,'FvFmsy is',round(FScenarios[f]/Fmsy$par,2),' MPAs.pdf'),family=Font,pointsize=12,width=6,height=4)
    par(mar=c(5.1,4.1,4.1,6.1),xpd=T)
    matplot(t(MPAs[,2:7]),type='l',col=rainbow(1.5*dim(MPAs)[1])[1:dim(MPAs)[1]],lty=1,xlab='Year',ylab='% NTZ',lwd=4,bty='n')
    legend('topright',inset=c(-.3,0),legend=MPANames,col=rainbow(1.5*dim(MPAs)[1])[1:dim(MPAs)[1]],lty=1,cex=.5)
    dev.off()
    
    PropNames<- NULL
    for (mm in 1:dim(MPAs)[1])
    {
      PropNames[mm]<- paste(round(MPAs[mm,2:dim(MPAs)[2]],2),collapse=':')
    }
    
    for (m in 1:dim(MPAs)[1])
    {
      
      show(paste('M is ',m))
      ResultStorage<- as.data.frame(matrix(NA,nrow=dim(MPAs)[2],ncol=3))
      colnames(ResultStorage)<- c('Yield','Biomass','Numbers')
      TempPop<- StartPop
      Patches<- BasePatches
      Patches$MPALocations<- c(1,0)
      
      cc<- cc+1
      for (y in 1:TimeToRun)
      {
        c<- c+1
        
        CurrentMPA<- MPAs[m,y]
        if (y>dim(MPAs)[2])
        {
          CurrentMPA<- MPAs[m,dim(MPAs)[2]]
          
          if (m==3)
          {
            CurrentMPA<- Rotations[y-dim(MPAs)[2]]
          }
          
        }
        TempPop<- AssignNTZ(CurrentMPA,TempPop)
        FTemp<- FScenarios[f] 
        FVec<- MoveFleet(FTemp,CurrentMPA,FleetSpill,0)
        PassPop<- GrowPopulation(TempPop,FVec,1,0,NA) #grow the new population
        TempPop<- PassPop$FinalNumAtAge
        ResultStorage[y,1]<- (PassPop$Performance$Yields[length(PassPop$Performance$Yields)])
        ResultStorage[y,2]<- (sum(WeightAtAge %*% PassPop$FinalNumAtAge))
        ResultStorage[y,3]<- (sum(PassPop$FinalNumAtAge))
        
         # NPBTimeSeries<- Discount(ResultStorage$Yield[1:y]-BaseConditions$Yield[f],Fleet$YieldDiscount,length(ResultStorage$Yield))$NPV#Yield balance through 10

         NPBTimeSeries<- Discount(ResultStorage$Yield[1:y]-BaseConditions$Yield[f],Fleet$YieldDiscount,y)$NPV#Yield balance through 10
			
        
        TotalStorage[c,]<-c(SpeciesList[s],m,f,FScenarios[f],CurrentMPA,y,PassPop$Performance$MeanYield,PassPop$Performance$MeanBiomass,PassPop$Performance$MeanNumbers, BaseConditions$Yield[f],BaseConditions$Numbers[f], BaseConditions$Biomass[f], OptimalConditions$Yield[f], OptimalConditions$Numbers[f], OptimalConditions$Biomass[f], OptNTZSize$par, ResultStorage$Yield[y]-BaseConditions$Yield[f] ,NPBTimeSeries)

} #Close TimeToRun loop   
      #Calculate relative performance
      BaseRelativeYield[m,f,]<- ResultStorage$Yield/BaseConditions$Yield[f]-1
      BaseRelativeNumbers[m,f,]<- ResultStorage$Numbers/BaseConditions$Numbers[f]-1
      BaseRelativeBiomass[m,f,]<- ResultStorage$Biomass/BaseConditions$Biomass[f]-1
      
      OptRelativeYield[m,f,]<- ResultStorage$Yield/OptimalConditions$Yield[f]-1
      OptRelativeNumbers[m,f,]<- ResultStorage$Numbers/OptimalConditions$Numbers[f]-1
      OptRelativeBiomass[m,f,]<- ResultStorage$Biomass/OptimalConditions$Biomass[f]-1
      
      RawYield[m,f,]<- ResultStorage$Yield
      RawNumbers[m,f,]<- ResultStorage$Numbers
      RawBiomass[m,f,]<- ResultStorage$Biomass
      
      ExperimentResults[m,f,1]<- Discount(ResultStorage$Yield,Fleet$YieldDiscount,length(ResultStorage$Yield))$NPV #Discounted NPV
      ExperimentResults[m,f,2]<- Discount(ResultStorage$Biomass,Fleet$BiomassDiscount,length(ResultStorage$Yield))$NPV #Discounted Biomass
      ExperimentResults[m,f,3]<- mean(ResultStorage$Yield) #Average yield
      ExperimentResults[m,f,4]<- mean(ResultStorage$Biomass) #Average biomass
      ExperimentResults[m,f,5]<- min(ResultStorage$Numbers) #Average numbers
      YearToYearChange<-(ResultStorage$Yield[2:length( ResultStorage$Yield)]-ResultStorage$Yield[1:(length( ResultStorage$Yield)-1)])/ResultStorage$Yield[1:(length( ResultStorage$Yield)-1)]
      YearToYearChange[is.finite(YearToYearChange)==F]<- 0
      ExperimentResults[m,f,6]<- 	sd(YearToYearChange)	
      ExperimentResults[m,f,7]<- 	mean(YearToYearChange)	 
      ExperimentResults[m,f,8]<- 	100*mean(BaseRelativeYield[m,f,]>0)	 
      ExperimentResults[m,f,9]<- 	100*mean(BaseRelativeNumbers[m,f,]>0)	 
      ExperimentResults[m,f,10]<- 	100*mean(BaseRelativeYield[m,f,])	 
      ExperimentResults[m,f,11]<- 	100*mean(BaseRelativeNumbers[m,f,])	 
      ExperimentResults[m,f,12]<- 	100*mean(BaseRelativeNumbers[m,f,]>0 & BaseRelativeYield[m,f,]>0)	 
      
      ExperimentResults[m,f,13]<- ResultStorage$Yield[6]-BaseConditions$Yield[f] #Yield balance in year 5
      ExperimentResults[m,f,14]<- ResultStorage$Biomass[6]-BaseConditions$Biomass[f] #Biomass balance in year 5
      ExperimentResults[m,f,15]<- Discount(ResultStorage$Yield[2:6]-BaseConditions$Yield[f],Fleet$YieldDiscount,5)$NPV #Yield balance through year 5
      
      ExperimentResults[m,f,16]<- ResultStorage$Yield[11]-BaseConditions$Yield[f] #Yield balance in year 10
      
      
      ExperimentResults[m,f,17]<- ResultStorage$Biomass[11]-BaseConditions$Biomass[f] #Biomass balance in year 10
      
      ExperimentResults[m,f,18]<- Discount(ResultStorage$Yield[2:11]-BaseConditions$Yield[f],Fleet$YieldDiscount,10)$NPV # through 10
      
      ExperimentResults[m,f,19]<- which((ResultStorage$Yield-BaseConditions$Yield[f])>=1)[1]-1 #Years to positive yield balance
      
      # ExperimentResults[m,f,19]<- which((ResultStorage$Yield/BaseConditions$Yield[f])>=1)[1]-1 #Years to positive yield balance
      
      
      ExperimentResults[m,f,20]<- which((ResultStorage$Biomass-BaseConditions$Biomass[f])>=1)[1]-1 #Years to positive biomass balance
      
      
      ExperimentResults[m,f,22]<- Discount(pmax(0,ResultStorage$Yield[2:EvalTime]-BaseConditions$Yield[f]),Fleet$YieldDiscount,EvalTime-1)$NPV #Positive balance over selected time horizon for evaluations     
      
      
      FinalNPBPositive<- sign(Discount(ResultStorage$Yield-BaseConditions$Yield[f],Fleet$YieldDiscount,length(ResultStorage$Yield))$NPV)>0 
      
      ExperimentResults[m,f,21]<- which((Discount(ResultStorage$Yield-BaseConditions$Yield[f],Fleet$YieldDiscount,length(ResultStorage$Yield))$CumValues)>=1)[1]-1 #Years until positive net yield balance          
      
      if (FinalNPBPositive == F)
      {
        ExperimentResults[m,f,21]<- NA #If the final NPB is negative, set years to positive NPB as NA
      }
      
      
      YieldBalance<- ResultStorage$Yield-BaseConditions$Yield[f]
      
      YieldBalance<- YieldBalance[2:TimeToRun]
      
      NegativeYears<- length(YieldBalance[YieldBalance<0])
      
      NegativeYields<- Discount(YieldBalance[YieldBalance<0],Fleet$YieldDiscount,NegativeYears)$NPV #Discounted Requested Loan Amount. 
      
      ExperimentResults[m,f,23]<- -NegativeYields    #Store Discounted Requested Loan Amount
      
      FlatExperimentResults[cc,1]<- SpeciesList[s]
      
      FlatExperimentResults[cc,2:dim(FlatExperimentResults)[2]]<- c(lh$Range,lh$BH.Steepness,m,f,FScenarios[f]/Fmsy$par,OptNTZSize$par,as.numeric(ExperimentResults[m,f,]))
      
    } #Close loop over MPA proposals
    
    
    pdf(file=paste(FigureFolder,'Yield Trajectoty ',f,'.pdf'),family=Font,pointsize=12,width=7,height=6)
    # sizemat<- matrix(1:2,nrow=1,ncol=2)
    # sizelay<- layout(mat=sizemat)
    par(mar=c(5.1,4.1,4.1,8.1),xpd=T)
    yMax<- max(RawYield[,f,], BaseConditions$Yield[f], OptimalConditions$Yield[f])
    yMin<- min(RawYield[,f,], BaseConditions$Yield[f], OptimalConditions$Yield[f])
    yMin<- .9*yMin
    matplot(t(RawYield[,f,1:11]),type='l',col=topo.colors(100)[seq(70,1,length.out=dim(MPAs)[1])],lwd=6,bty='n',xlab='Time',ylab='Fishing Yields',main=paste('F is ',round(FScenarios[f],4)),ylim=c(yMin,yMax),lty=c(1:dim(MPAs)[1]))
    abline(h=c(OptimalConditions$Yield[f]),lty=2,xpd=F)
    legend('topright',inset=c(-.35,0),legend=c(MPANames,'Long Term Optimal EQ'),bty='y',col=c(topo.colors(100)[seq(70,1,length.out=dim(MPAs)[1])],1,1),lty=c(1:dim(MPAs)[1],2),cex=0.6)
    dev.off()
    
    pdf(file=paste(FigureFolder,'Biomass Trajectoty ',f,'.pdf'),family=Font,pointsize=12,width=7,height=6)
    par(mar=c(5.1,4.1,4.1,8.1),xpd=T)
    
    yMax<- max(RawBiomass[,f,], BaseConditions$Biomass[f], OptimalConditions$Biomass[f])
    yMin<- min(RawBiomass[,f,], BaseConditions$Biomass[f], OptimalConditions$Biomass[f])
    yMin<- .9*yMin
    matplot(t(RawBiomass[,f,1:11]),type='l',col=topo.colors(100)[seq(70,1,length.out=dim(MPAs)[1])],lwd=6,bty='n',xlab='Time',ylab='Biomass',main=paste('F is ',round(FScenarios[f],4)),ylim=c(yMin,yMax),lty=c(1:dim(MPAs)[1]))
    abline(h=c(OptimalConditions$Biomass[f]),lty=2,xpd=F)
    legend('topright',inset=c(-.35,0),legend=c(MPANames,'Long Term Optimal EQ'),bty='t',col=c(topo.colors(100)[seq(70,1,length.out=dim(MPAs)[1])],1),lty=c(1:dim(MPAs)[1],2),cex=0.6)
    dev.off()
    
    
    pdf(file=paste(FigureFolder,'Yield Balance Trajectoty ',f,'.pdf'),family=Font,pointsize=12,width=7,height=6)
    par(mar=c(5.1,4.1,4.1,8.1),xpd=T)
    
    matplot(t(RawYield[,f,1:11]-BaseConditions$Yield[f]),type='l',col=topo.colors(100)[seq(70,1,length.out=dim(MPAs)[1])],lwd=6,bty='n',xlab='Time',ylab='Surplus Fishing Yields',main=paste('F is ',round(FScenarios[f],4)),lty=1:dim(MPAs)[1])
    abline(h=c(OptimalConditions$Yield[f]-BaseConditions$Yield[f]),lty=2,xpd=F)
    legend('topright',inset=c(-.35,0),legend=c(MPANames,'Long Term Optimal EQ'),bty='t',col=c(topo.colors(100)[seq(70,1,length.out=dim(MPAs)[1])],1),lty=c(1:dim(MPAs)[1],2),cex=0.6)
    dev.off()
    
    TempSurplus<- RawYield[,f,]-BaseConditions$Yield[f]
    
    for (tt in 1:dim(TempSurplus)[1])
    {
      TempSurplus[tt,]<- cumsum(TempSurplus[tt,])
    }
    
    
    pdf(file=paste(FigureFolder,'NPB Trajectoty ',f,'.pdf'),family=Font,pointsize=12,width=7,height=6)
    par(mar=c(5.1,4.1,4.1,8.1),xpd=T)
    matplot(t(TempSurplus[,1:11]),type='l',col=topo.colors(100)[seq(70,1,length.out=dim(MPAs)[1])],lwd=6,bty='n',xlab='Time',ylab='NPB',main=paste('F is ',round(FScenarios[f],4)),lty=c(1:dim(MPAs)[1]))
    # abline(h=c(OptimalConditions$Yield[f]-BaseConditions$Yield[f]),lty=1:2)
    legend('topright',inset=c(-.35,0),legend=c(MPANames),bty='n',col=c(topo.colors(100)[seq(70,1,length.out=dim(MPAs)[1])]),lty=c(1:dim(MPAs)[1]),cex=.6)
    dev.off()  
    
    a=RawYield[,f,]-BaseConditions$Yield[f]

  } #Close loop over fishing scenarios
  
 
  #######Perform Loan Analysis ######
  
  LoanRates<- seq(.01,2,length.out=20)
  LoanTime<- 10
  
  PlotLayout<- matrix(1:length(FScenarios),nrow=length(FScenarios),ncol=1)
  
  Colors=terrain.colors(2*(dim(MPAs)[1]))[1:(dim(MPAs)[1])]
  
  
  for (d in 1:dim(FlatExperimentResults)[1])
  {
  	
  	
  	MaxInterestRate<- NA
  	
  	if (as.numeric(FlatExperimentResults$RequestedLoan[d])>0)
  	{
  		
  	
  	MaxInterestRate<- optim(.05,FindMaxInterestRate,LoanTime=LoanTime,LoanAmount=as.numeric(FlatExperimentResults$RequestedLoan[d]),Surplus=as.numeric(FlatExperimentResults$TenYearNPSurplus[d]),lower=-5,upper=10,method='Brent')$par
  	
  	}
 
  	FlatExperimentResults$MaxInterestRate[d]<- 100*exp(MaxInterestRate)
  	
  }
  
  
  
  ##### Store Things ####
  
  FlatExperimentResults$MaxInterestRate[FlatExperimentResults$MaxInterestRate< exp(-4)]<- 0

  FlatExperimentResults$MaxInterestRate[FlatExperimentResults$MaxInterestRate> 2]<- Inf

  FlatExperimentResults$MaxInterestRate<- round(FlatExperimentResults$MaxInterestRate,2)
  
   FlatExperimentResults$MPAScenario<- as.factor(FlatExperimentResults$MPAScenario)
  
  levels(FlatExperimentResults$MPAScenario)<- MPANames
  
  write.csv(file=paste(ResultFolder,'Total Results.csv'),TotalStorage)
  
  write.csv(file=paste(ResultFolder,'Experiment Results.csv'),FlatExperimentResults)
  
  TotalStorage$Species<- SpeciesList[s]
  FlatExperimentResults$Species<- SpeciesList[s]
  
  AllSpeciesStorage<- rbind(AllSpeciesStorage,TotalStorage)
  AllSpeciesExperimentResults<- rbind(AllSpeciesExperimentResults, FlatExperimentResults)
  
  
  
  ####### PLOT LOTS OF STUFF ########
  

  #### Make species specific 3d plots #####
  
  # pdf(file=paste(FigureFolder,'10 year NPV Balance Wireplot.pdf',sep=''),pointsize=14,width=12,height=8) 
  
  # print(wireframe( TenyearNPVBalance ~ MPAScenario * FvFmsy , data = FlatExperimentResults,scales=list(arrows=F) ,shade=T,zlab='NPB'))
  # dev.off()
  
  
  # pdf(file=paste(FigureFolder,'10 year Yield Balance Wireplot.pdf',sep=''),pointsize=14,width=12,height=8) 
  # print(wireframe( TenyearYieldBalance ~ MPAScenario * FvFmsy , data = FlatExperimentResults,scales=list(arrows=F) ,shade=T,zlab='Surplus'))
  # dev.off()
  
  # pdf(file=paste(FigureFolder,'10 year Biomass Balance Wireplot.pdf',sep=''),pointsize=14,width=12,height=8) 
  
  # print(wireframe( TenyearBiomassBalance ~ MPAScenario * FvFmsy , data = FlatExperimentResults,scales=list(arrows=F) ,shade=T,zlab='Biomass Surplus'))
  # dev.off()
  
  
  # pdf(file=paste(FigureFolder,'Year to Yield Recovery Cloud.pdf',sep=''),pointsize=14,width=12,height=8) 
  # print(cloud( YearsToYieldRecovery ~ MPAScenario * FvFmsy , data = FlatExperimentResults,scales=list(arrows=F),zlab='Year to Positive Yield' ))
  # dev.off()
  
  # pdf(file=paste(FigureFolder,'Years to Balance Recovery.pdf',sep=''),pointsize=14,width=12,height=8) 
  # print(cloud(YearsToBalanceRecovery ~ MPAScenario * FvFmsy , data = FlatExperimentResults,scales=list(arrows=F) ,shade=F,na.rm=T))
  # dev.off()
  
  
  
  write.csv(OptRelativeNumbers,file=paste(ResultFolder,'Numbers relative to optimum numbers.csv'))
  
  ResultNames<- dimnames(ExperimentResults)[3]
  ResultNames<- ResultNames[[1]]
  
  FNames<- dimnames(ExperimentResults)[2]
  FNames<- FNames[[1]]
  
  Comparisons<- c(1,2,3,4,3,5,4,6,3,6,4,10,3,11)
  TradeOffs<- matrix(Comparisons,ncol=2,nrow=length(Comparisons)/2,byrow=T)
  
  
  ### Make Paper Time Series Plots ###
  
  for (m in 1:dim(MPAs)[1])
  {
    pdf(file=paste(FigureFolder,'Scenario ',m,' Yields.pdf',sep=''),family=Font,pointsize=12,width=6,height=4)
    par(mar=c(5.1,4.1,4.1,7.1),xpd=T)
    
    matplot(t(RawYield[m,,]),type='l',lty=1,lwd=4,col=topo.colors(1.2*dim(MPAs)[1])[1:dim(MPAs)[1]],xlab='Year',ylab='Yields',bty='n',main=paste("Scenario ",m))
    legend('topright',inset=c(-.2,0),legend=round(FScenarios/Fmsy$par,2),title='F/Fmsy',col=topo.colors(1.2*dim(MPAs)[1])[1:dim(MPAs)[1]],lty=1,lwd=2,bty='n',cex=.6)
    dev.off()
    
    pdf(file=paste(FigureFolder,'Scenario ',m,' Biomass.pdf',sep=''),family=Font,pointsize=12,width=6,height=4)
    par(mar=c(5.1,4.1,4.1,7.1),xpd=T)
    
    matplot(t(RawBiomass[m,,]/sum(lh$Bmsy)),type='l',lty=1,lwd=4,col=topo.colors(1.2*dim(MPAs)[1])[1:dim(MPAs)[1]],xlab='Year',ylab='B/Bmsy',bty='n',main=paste("Scenario ",m))
    legend('topright',inset=c(-.2,0),legend=round(FScenarios/Fmsy$par,2),title='F/Fmsy',col=topo.colors(1.2*dim(MPAs)[1])[1:dim(MPAs)[1]],lty=1,lwd=2,bty='n',cex=.6)
    dev.off()
    
  }
  
  
  for (t in 1:dim(TradeOffs)[1])
  {
    
    TOff<- TradeOffs[t,]
    pdf(file=paste(FigureFolder,ResultNames[TOff[1]], ' vs ',ResultNames[TOff[2]],'.pdf',sep=''),family=Font,pointsize=12,width=7,height=6)
    par(mai=c(1,1,1,2),xpd=T)
    matplot(ExperimentResults[,,TOff[1]], ExperimentResults[,,TOff[2]],type='p',pch=19,xlab=ResultNames[TOff[1]],ylab=ResultNames[TOff[2]],bty='n',col=terrain.colors(100)[seq(70,1,length.out=length(FNames))])
    legend('topright',legend=FNames,pch=19,col=terrain.colors(100)[seq(70,1,length.out=length(FNames))],bty='n',inset=c(-.5,0))
    dev.off()
  }
  
} #Close species list loop


AllSpeciesExperimentResults$MPAScenario<- as.factor(AllSpeciesExperimentResults$MPAScenario)

levels(AllSpeciesExperimentResults$MPAScenario)<- MPANames

AllSpeciesExperimentResults<- AllSpeciesExperimentResults[AllSpeciesExperimentResults$MPAScenario!='Status Quo',]

TotalStorage$Year<- as.numeric(TotalStorage$Year)

TotalStorage$NPB<- as.numeric(TotalStorage$NPB)

TotalStorage$m<- as.factor(TotalStorage$m)

TotalStorage$f<- as.factor(TotalStorage$f)

levels(TotalStorage$m)<- MPANames

levels(TotalStorage$f)<- c('F/Fmsy=1.5','F/Fmsy=3')

AllSpeciesStorage$Year<- as.numeric(AllSpeciesStorage$Year)

AllSpeciesStorage$NPB<- as.numeric(AllSpeciesStorage$NPB)

AllSpeciesStorage $m<- as.factor(AllSpeciesStorage$m)

AllSpeciesStorage $f<- as.factor(AllSpeciesStorage$f)

levels(AllSpeciesStorage $m)<- MPANames

levels(AllSpeciesStorage $f)<- c('F/Fmsy=1.5','F/Fmsy=3')


SubScenarios<- c('EqNTZ','GNTZ','SNTZ','OptNTZ')

Where<-   TotalStorage$m %in% SubScenarios

PlotStorage<- AllSpeciesStorage[Where,]

PlotStorage$m<- as.character(PlotStorage$m)

PlotStorage$m<- as.factor(PlotStorage$m)

PlotStorage$FracNTZ<- as.numeric(PlotStorage$FracNTZ)

PlotStorage$PositiveYields<- 0

FNames<- levels(TotalStorage$f)

for (s in 1:length(SpeciesList))
{

for (m in 1:length(SubScenarios))
{
	
	for (ff in 1:length(FNames))
	{
		
		
		WherePositiveYields<- which(PlotStorage$Species==SpeciesList[s]& PlotStorage$m==SubScenarios[m] & PlotStorage$f==FNames[ff] & PlotStorage$YieldBalance>=0 & PlotStorage$Year>1)[1]
							
		PlotStorage$PositiveYields[WherePositiveYields]<- 1

	}
}

}
PlotStorage$PosYears<- PlotStorage$Year * PlotStorage$PositiveYields


PlotStorage$PosNPB<- PlotStorage$NPB * PlotStorage$PositiveYields


save.image(file=paste(ResultFolder,'Completed Workspace.rdata'))

} #Run Analysis Loop
if (RunAnalysis==0)
{
load(paste(BatchFolder,'Completed Workspace.rdata'))
}
# trellis.par.get()
Cols<-   scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65, h.start = 0, direction = 1)

# pdf(file=paste(FigureFolder,'NPB Trajectory.pdf',sep=''),col='rgb')
# xyplot(NPB ~ Year | f,group=m,data=PlotStorage,type='l',subset=(Year<=11),par.settings=simpleTheme(pch=20,lty=1:length(SubScenarios)),auto.key=list(space='right',points=F,lines=T),layout=c(1,2),index.cond=list(c(2,1)),lwd=4,drop.unused.levels=T,panel=function(x,y,...)
# {panel.xyplot(x,y,...)
# panel.abline(h=0,lty=2)
# # panel.points(PlotStorage$PosYears ,PlotStorage$PosNPB,type='p',cex=10,...)
# }) 


# # # + layer(panel.xyplot(PlotStorage$PosNPB ~ PlotStorage$PosYears | PlotStorage$f,pch=5,subscripts=T,layout=c(1,2),index.cond=list(c(2,1))))
# dev.off()

###Not as slick as xyplot, but gives you full control over what is plotted
for (s in 1:length(SpeciesList))
{

FigureFolder<- paste(BatchFolder,SpeciesList[s],'/Figures/',sep='')
ResultFolder<- paste(BatchFolder,SpeciesList[s],'/Results/',sep='')

plotyears=30  #no. years to plot
pdf(file=paste(FigureFolder,'NPB Trajectory.pdf',sep=''),col='rgb')
#subset data for each scenario (m)
m1=subset(PlotStorage,m=="EqNTZ" & Species== SpeciesList[s])
m2=subset(PlotStorage,m=="GNTZ" & Species == SpeciesList[s])
m3=subset(PlotStorage,m=="SNTZ"& Species == SpeciesList[s])

m4=subset(PlotStorage,m=="OptNTZ"& Species == SpeciesList[s])

# m4G=subset(PlotStorage,m=="OptPath"& Species == SpeciesList[2])

#get points where yield equals status quo
 y1=PlotStorage[PlotStorage$f=="F/Fmsy=1.5" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s],c(6,18) ]
 y2=PlotStorage[PlotStorage$f=="F/Fmsy=3" & PlotStorage$PositiveYields==1 & PlotStorage$Species== SpeciesList[s],c(6,18)]

par(mfrow=c(2,1),mar=c(0,0,0,0),oma=c(4,6,2,6))
 
plot(m1$Year[as.vector(m1$f)=="F/Fmsy=1.5"],m1$NPB[as.vector(m1$f)=="F/Fmsy=1.5"],
  type="l",ylab="",xlab="",ylim=c(min(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]]),1.25*max(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]])),
  xlim=c(1,plotyears),col="dodgerblue",lwd=4,las=1,xaxt="n" )
mtext("F/Fmsy = 1.5",line=-1)
text(1,max(PlotStorage$NPB)*0.99,"A.")
lines( m2$Year[as.vector(m2$f)=="F/Fmsy=1.5"],m2$NPB[as.vector(m2$f)=="F/Fmsy=1.5"],col=2,lwd=2,lty=3)
lines( m3$Year[as.vector(m3$f)=="F/Fmsy=1.5"],m3$NPB[as.vector(m3$f)=="F/Fmsy=1.5"],col="medium orchid",lwd=2,lty=4)
lines( m4$Year[as.vector(m4$f)=="F/Fmsy=1.5"],m4$NPB[as.vector(m4$f)=="F/Fmsy=1.5"],col="darkgreen",lwd=2,lty=5)
abline(h=0,lty=2)
points(y1[,1],y1[,2],col=c("dodgerblue","medium orchid",2,"darkgreen"),pch=16,cex=1.5)


plot(m1$Year[as.vector(m1$f)=="F/Fmsy=3"],m1$NPB[as.vector(m1$f)=="F/Fmsy=3"],
  type="l",ylab="",xlab="Year",ylim=c(min(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]]),1.25*max(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]])),
  xlim=c(1,plotyears),col="dodgerblue",lwd=4 ,las=1)
mtext("F/Fmsy = 3.0",line=-1)
text(1,max(PlotStorage$NPB)*0.99,"B.")
lines( m2$Year[as.vector(m2$f)=="F/Fmsy=3"],m2$NPB[as.vector(m2$f)=="F/Fmsy=3"],col=2,lwd=2,lty=3)
lines( m3$Year[as.vector(m3$f)=="F/Fmsy=3"],m3$NPB[as.vector(m3$f)=="F/Fmsy=3"],col="medium orchid",lwd=2,lty=4)
lines( m4$Year[as.vector(m4$f)=="F/Fmsy=3"],m4$NPB[as.vector(m4$f)=="F/Fmsy=3"],col="darkgreen",lwd=2,lty=5)
abline(h=0,lty=2)         
points(y2[,1],y2[,2],col=c("dodgerblue","medium orchid",2,"darkgreen"),pch=16,cex=1.5)

par(new=TRUE, mfrow=c(1,1), oma=c(0,0,0,0))

plot(1,
  xlim=c(0,1),
  ylim=c(0,1),
  xaxs="i",
  yaxs="i",
  type="n",
  axes=FALSE)


legend(0.82,0.9, legend= SubScenarios, lty=c(1,3,4,5),
col=c("dodgerblue",2,"medium orchid","darkgreen"),lwd=2,bty="n",cex=0.8)

 text(.5,.05,"Year",font=2)

 text(.03,.5,"Net Present Balance",srt=90,font=2)

dev.off()

pdf(file=paste(FigureFolder,'FracNTZ Trajectory.pdf',sep=''))
print(xyplot(FracNTZ ~ Year | f,group=m,data=PlotStorage,type='l',subset=(Year<=11 & PlotStorage$Species==SpeciesList[s]),auto.key=list(space='right',points=F,lines=T),layout=c(1,2),index.cond=list(c(2,1)),lwd=4,drop.unused.levels=T,panel=function(x,y,...)
{panel.xyplot(x,y,...)
panel.abline(h=0,lty=2)
},))
dev.off()

} #Close Species Loop

# pdf(file=paste(BatchFolder,'Years to Balance Recovery.pdf',sep=''))
# print(cloud(YearsToBalanceRecovery ~ Movement * FvFmsy  | MPAScenario,data = AllSpeciesExperimentResults,scales=list(arrows=F) ,shade=T,na.rm=T,pretty=T,ylab=c('FvFmsy'),xlab='Movement',zlab='',main='Years to Balance Recovery',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1),cex=2))
# dev.off()

# pdf(file=paste(BatchFolder,'Years to Balance Recovery Wireframe.pdf',sep=''))
# print(wireframe(YearsToBalanceRecovery ~ Movement * FvFmsy  | MPAScenario,data = AllSpeciesExperimentResults,scales=list(arrows=F) ,shade=T,na.rm=T,pretty=T,ylab=c('FvFmsy'),xlab='Movement',zlab='',main='Years to Balance Recovery',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1),cex=2))
# dev.off()

# pdf(file=paste(BatchFolder,'Years to Yield Recovery.pdf',sep=''))

# print(cloud(YearsToYieldRecovery ~ Movement * FvFmsy  | MPAScenario,data = AllSpeciesExperimentResults,scales=list(arrows=F) ,shade=T,na.rm=T,pretty=T,ylab=c('FvFmsy'),xlab='Movement',zlab='',main='Years to Yield Recovery',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1),cex=2))
# dev.off()

# pdf(file=paste(BatchFolder,'Years to Yield Recovery Wireframe.pdf',sep=''))

# print(wireframe(YearsToYieldRecovery ~ Movement * FvFmsy  | MPAScenario,data = AllSpeciesExperimentResults,scales=list(arrows=F) ,shade=T,na.rm=T,pretty=T,ylab=c('FvFmsy'),xlab='Movement',zlab='',main='Years to Yield Recovery',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1),cex=2))
# dev.off()


# pdf(file=paste(BatchFolder,'10 Year NPV Balance.pdf',sep=''))

# print(wireframe(TenyearNPVBalance ~ Movement * FvFmsy  | MPAScenario,data = AllSpeciesExperimentResults,scales=list(arrows=F) ,shade=T,na.rm=T,pretty=T,ylab=c('FvFmsy'),xlab='Movement',zlab='',main='10y NPV Balance',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1)))
# dev.off()

# pdf(file=paste(BatchFolder,'10 Year NPV Balance using Steepness.pdf',sep=''))

# print(wireframe(TenyearNPVBalance ~ Steepness * FvFmsy  | MPAScenario,data = AllSpeciesExperimentResults,scales=list(arrows=F) ,shade=T,na.rm=T,pretty=T,ylab=c('FvFmsy'),xlab='Steepness',zlab='',main='10y NPV Balance',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1)))
# dev.off()

# pdf(file=paste(BatchFolder,'10 Year Biomass Balance.pdf',sep=''))

# print(wireframe(TenyearBiomassBalance ~ Movement * FvFmsy  | MPAScenario,data = AllSpeciesExperimentResults,scales=list(arrows=F) ,shade=T,na.rm=T,pretty=T,ylab=c('FvFmsy'),xlab='Movement',zlab='',main='10y Biomass Balance',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1)))
# dev.off()

# pdf(file=paste(BatchFolder,'10 Year Yield Balance.pdf',sep=''))

# print(wireframe(TenyearYieldBalance ~ Movement * FvFmsy  | MPAScenario,data = AllSpeciesExperimentResults,scales=list(arrows=F) ,shade=T,na.rm=T,pretty=T,ylab=c('FvFmsy'),xlab='Movement',zlab='',main='10y Yield Balance',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1)))
# dev.off()


# Colors=colorRampPalette(c('green','red'))(8)

# pdf(file=paste(BatchFolder,'Histogram of Years to Balance Recovery.pdf',sep=''))
# print(histogram(~YearsToBalanceRecovery | MPAScenario,data = AllSpeciesExperimentResults,scales=list(arrows=F) ,shade=T,na.rm=T,pretty=F,xlab=c('Years to Positive Net Present Yield Balance'),ylab='Frequency',zlab='',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1),type='count',col=Colors))
# dev.off()


# Colors=colorRampPalette(c('green','red'))(8)

# pdf(file=paste(BatchFolder,'Histogram of Ten Year NPB.pdf',sep=''))
# print(histogram(~ TenyearNPVBalance | MPAScenario,data = AllSpeciesExperimentResults,scales=list(arrows=F) ,shade=T,na.rm=T,pretty=F,xlab=c('10 Year NPB'),ylab='Frequency',zlab='',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1),type='count',col=Colors))
# dev.off()


# BoxData<- AllSpeciesExperimentResults[AllSpeciesExperimentResults$MPAScenario!='Status Quo',]

# BoxData$MPAScenario<- as.character(BoxData$MPAScenario)

# levels(BoxData$MPAScenario)<- MPANames[2:length(MPANames)]


# OrderCats<- order(as.numeric(by(BoxData$YearsToBalanceRecovery,BoxData$MPAScenario,median,na.rm=T)))

# BoxData$MPAScenario<- ordered(BoxData$MPAScenario,levels=levels(BoxData$MPAScenario)[OrderCats])

# pdf(file=paste(BatchFolder,'Boxplot of Years to Balance Recovery.pdf',sep=''))
# par(bty='n')
# boxplot(YearsToBalanceRecovery ~MPAScenario,data=BoxData,notch=T,varwidth=T,outline=F,bty='n',ylab='Years to Positive NPB',col='grey')
# dev.off()


# OrderCats<- order(as.numeric(by(BoxData$YearsToYieldRecovery,BoxData$MPAScenario,median,na.rm=T)))

# BoxData$MPAScenario<- ordered(BoxData$MPAScenario,levels=levels(BoxData$MPAScenario)[OrderCats])


# pdf(file=paste(BatchFolder,'Boxplot of Years to Yield Recovery.pdf',sep=''))
# par(bty='n')
# boxplot(YearsToYieldRecovery ~MPAScenario,data=BoxData,notch=T,varwidth=T,outline=F,bty='n',ylab='Years to Positive Yield Balance',col='grey')
# dev.off()


# OrderCats<- order(as.numeric(by(BoxData$TenyearNPVBalance,BoxData$MPAScenario,median,na.rm=T)))

# BoxData$MPAScenario<- ordered(BoxData$MPAScenario,levels=levels(BoxData$MPAScenario)[OrderCats])


# pdf(file=paste(BatchFolder,'Boxplot of 10 year NPB.pdf',sep=''))
# par(bty='n')
# boxplot(TenyearNPVBalance ~MPAScenario,data=BoxData,notch=T,varwidth=T,outline=F,bty='n',ylab='NPB in Year 10',col='grey')
# dev.off()


# pdf(file=paste(BatchFolder,'Histogram of Years to Balance Recovery no 1FvFmsy.pdf',sep=''))
# print(histogram(~YearsToBalanceRecovery | MPAScenario,data = AllSpeciesExperimentResults[AllSpeciesExperimentResults$FvFmsy!=1,],scales=list(arrows=F) ,shade=T,na.rm=T,pretty=F,xlab=c('Years to Positive Net Present Yield Balance'),ylab='Frequency',zlab='',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1),type='count',col=Colors))
# dev.off()

# pdf(file=paste(BatchFolder,'Histogram of Years to Yield Recovery.pdf',sep=''))
# print(histogram(~YearsToYieldRecovery | MPAScenario,data = AllSpeciesExperimentResults,scales=list(arrows=F) ,shade=T,na.rm=T,pretty=F,xlab=c('Years to Positive Yield Balance'),ylab='Frequency',zlab='',strip = strip.custom(strip.levels = c(TRUE,TRUE),style=1),type='count',col=Colors))
# dev.off()


write.csv(file=paste(BatchFolder,'All Species Experiment Results.csv',sep=''),AllSpeciesExperimentResults)

AllSpeciesStorage$m<- as.factor(AllSpeciesStorage$m)
levels(AllSpeciesStorage$m)<- MPANames

write.csv(file=paste(BatchFolder,'All Species All Results.csv',sep=''), AllSpeciesStorage)

save.image(file=paste(BatchFolder,'Completed Workspace.rdata'))



# ShortMPANames<- MPANames[MPANames!='Status Quo']

# AllSpeciesStorage<- AllSpeciesStorage[AllSpeciesStorage$m!='Status Quo',]

# MeanYearMatrix<- as.data.frame(matrix(NA,nrow=TimeToRun,ncol=length(ShortMPANames)))
# colnames(MeanYearMatrix)<- ShortMPANames

# TopYearMatrix<- as.data.frame(matrix(NA,nrow=TimeToRun,ncol=length(ShortMPANames)))
# colnames(TopYearMatrix)<- ShortMPANames

# BotYearMatrix<- as.data.frame(matrix(NA,nrow=TimeToRun,ncol=length(ShortMPANames)))
# colnames(BotYearMatrix)<- ShortMPANames



# for (m in 1:length(ShortMPANames))
# {
  
  # TempData<- AllSpeciesStorage[AllSpeciesStorage$m==ShortMPANames[m],]
  # for (y in 1:TimeToRun)
  # {
    # Where<- TempData$Year==y
    
    # SummaryData<- boxplot(as.numeric(TempData$NPB[Where]),outline=F,plot=F)
    
    # MeanYearMatrix[y,m]<- mean(as.numeric(TempData$NPB[Where]),na.rm=T)
    
    # # MeanYearMatrix[y,m]<- SummaryData$stats[3]
    
    # BotYearMatrix[y,m]<- SummaryData$stats[1]
    
    # TopYearMatrix[y,m]<- SummaryData$stats[5]
    
    
  # }
# }

# # MeanYearMatrix<- cbind(1:TimeToRun,MeanYearMatrix)
# # colnames(MeanYearMatrix)[1]<- 'Year'
# # Test<- reshape(MeanYearMatrix,ids=colnames(MeanYearMatrix), direction='long',varying='Year')

# RunWindow<- 15

# MeanTimeStack<- stack(MeanYearMatrix[1:RunWindow,])

# BotTimeStack<- stack(BotYearMatrix[1: RunWindow,])

# TopTimeStack<- stack(TopYearMatrix[1: RunWindow,])

# TotalStack<- cbind(MeanTimeStack,BotTimeStack[,1],TopTimeStack[,1])

# colnames(TotalStack)<- c('values','MPAStrategy','Bottom','Top')

# MinVal<- signif(ceiling(min(TotalStack$Bottom,na.rm=T)),1)

# MaxVal<- signif(ceiling(max(TotalStack$Top,na.rm=T)),1)


# pdf(file=paste(BatchFolder,'Mean Time Trend of NPB2.pdf',sep=''),width=10,height=6)
# l<- matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T)
# l2<- layout(mat=l,widths=c(3,3,3))

# for (r in 1:length(ShortMPANames)) {
  
  # Where<- TotalStack$MPAStrategy==ShortMPANames[r]
  # plot(1: RunWindow,TotalStack$values[Where],bty='n',type='l',lwd=4,col='dodgerblue',ylim=c(MinVal,MaxVal))
  # lines(1:RunWindow,TotalStack$Bottom[Where],lty=2)
  # lines(1:RunWindow,TotalStack$Top[Where],lty=2)
  
  # title(main=ShortMPANames[r])
  
# }
# dev.off()



# pdf(file=paste(BatchFolder,'Mean Time Trend of NPB.pdf',sep=''),width=10,height=6)
# xyplot(values ~ rep((1:15),length(levels(TotalStack$MPAStrategy))) | MPAStrategy ,data=TotalStack,pretty=T,xlab='Time',ylab='NPB',type='l',lwd=4)
# dev.off()


# for (m in 1:length(ShortMPANames))
# {
  
  # TempData<- AllSpeciesStorage[AllSpeciesStorage$m==ShortMPANames[m],]
  # for (y in 1:TimeToRun)
  # {
    # Where<- TempData$Year==y
    
    # SummaryData<- boxplot(as.numeric(TempData$YieldBalanc[Where]),outline=F,plot=F)
    
    # MeanYearMatrix[y,m]<- mean(as.numeric(TempData$YieldBalance[Where]),na.rm=T)
    
    # BotYearMatrix[y,m]<- SummaryData$stats[1]
    
    # TopYearMatrix[y,m]<- SummaryData$stats[5]
    
    
  # }
# }

# # MeanYearMatrix<- cbind(1:TimeToRun,MeanYearMatrix)
# # colnames(MeanYearMatrix)[1]<- 'Year'
# # Test<- reshape(MeanYearMatrix,ids=colnames(MeanYearMatrix), direction='long',varying='Year')

# MeanTimeStack<- stack(MeanYearMatrix[1:15,])

# BotTimeStack<- stack(BotYearMatrix[1:15,])

# TopTimeStack<- stack(TopYearMatrix[1:15,])

# TotalStack<- cbind(MeanTimeStack,BotTimeStack[,1],TopTimeStack[,1])

# colnames(TotalStack)<- c('values','MPAStrategy','Bottom','Top')

# pdf(file=paste(BatchFolder,'Mean Time Trend of Yield Balance.pdf',sep=''),width=10,height=6)
# xyplot(values ~ rep((1:15),length(levels(TotalStack$MPAStrategy))) | MPAStrategy ,data=TotalStack,pretty=T,xlab='Time',ylab='Yield Balance',type='l',lwd=4)
# dev.off()

