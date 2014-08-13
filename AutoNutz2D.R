rm(list=ls())
library(lattice)
require(stats)
library(proftools)
library(plyr)
####### Shrinking MPAs #########
#Wilson, Doughtery, Ovando et al something something

####### READ IN CONTROL FILE ########

# setwd('/Users/danovando/Dropbox/Shrinking NTZ')

RunAnalysis<- 0

DataNames<- c('NPV of Yield','NPV of Biomass','Mean Yield','Mean Biomass','Mean Numbers','Yield Instability','Mean Changes in Yield','Percent Years with Profit Gains','Percent Years with Numbers Gains','Mean Percent Change in Yield','Mean Percent Change in Numbers','Percent Years With Numbers and Yield Gains','FiveyearYieldBalance','FiveyearBiomassBalance','FiveyearNPVBalance','TenyearYieldBalance','TenyearBiomassBalance','TenyearNPVBalance','YearsToYieldRecovery','YearsToBioRecovery','YearsToBalanceRecovery','TenYearNPSurplus','RequestedLoan','MaxInterestRate')

LongDataNames<- c('Species','Movement','Steepness','MPAScenario','FishingScenario','FvFmsy','OptNTZ','NPV of Yield','NPV of Biomass','Mean Yield','Mean Biomass','Mean Numbers','Yield Instability','Mean Changes in Yield','Percent Years with Profit Gains','Percent Years with Numbers Gains','Mean Percent Change in Yield','Mean Percent Change in Numbers','Percent Years With Numbers and Yield Gains','FiveyearYieldBalance','FiveyearBiomassBalance','FiveyearNPVBalance','TenyearYieldBalance','TenyearBiomassBalance','TenyearNPVBalance','YearsToYieldRecovery','YearsToBioRecovery','YearsToBalanceRecovery','TenYearNPSurplus','RequestedLoan','MaxInterestRate')

AllSpeciesStorage<- as.data.frame(matrix(NA,nrow=0,ncol=18))

colnames(AllSpeciesStorage)<- c('Species','m','f','Frate','FracNTZ','Year','Yield','Biomass','Numbers','SQYield','SQNumbers','SQBiomass','OptNTZYield','OptNTZNumbers','OptNTZBiomass','OptNTZ','YieldBalance','NPB')

AllSpeciesExperimentResults<- as.data.frame(matrix(NA,nrow=0,ncol=length(LongDataNames)))

colnames(AllSpeciesExperimentResults)<- c(LongDataNames)

BatchFolder<- 'Results/August_12_14_Wrapped/'
dir.create(paste(BatchFolder,sep=''))

# MPANames<- c('Status Quo','EqNTZ','Rotate','SNTZ','Basic','GNTZ','OptNTZ')

MPANames<- c('StatusQuo','EqNTZ','SNTZ','GNTZ','OptNTZ')

LifeHistories<- read.csv('Inputs/Life History Inputs.csv')

# LifeHistories<- LifeHistories[7,]

LifeColumns<- colnames(LifeHistories)

LifeVars<- c('Range','MaxAge','m','k','Linf','t0','AgeMa50','LengthMa50','MaturityMode','BH.Steepness','wa','wb','WeightForm','fa','fb')

LifeHistories$SciName<- as.character(levels(LifeHistories$SciName))[LifeHistories$SciName]

LifeHistories$CommName<- as.character(levels(LifeHistories$CommName))[LifeHistories$CommName]

LifeHistories$WeightForm<- as.character(levels(LifeHistories$WeightForm))[LifeHistories$WeightForm]

SpeciesList<- LifeHistories$CommName

if (RunAnalysis==1)
{
  
  for (s in 1:length(SpeciesList))
  {
    Species<- SpeciesList[s] #species file to load
    
    StoreRun<- paste(BatchFolder,Species,sep='') #1 if you want to create a seeded folder to store results, 0 if you want it in the generic working folder
    
    source('SHRINKNTZ_CONTROLFILE2D.R')
    
    BasePatches<- Patches
    
    LifeHistory<- LifeHistories[s,]
    
    for (l in 1:length(LifeVars))
    {
      
      Where1<- LifeColumns ==LifeVars[l]
      
      Where2<- names(lh)==LifeVars[l]
      
      lh[Where2]<- LifeHistory[Where1]
      
    }
    
    lh$LengthMa95<-  1.01* lh$LengthMa50
    
    lh$AgeMa95<-  1.01* lh$AgeMa50
    
    lh$MoveType<- '2D'
    
    lh$MovementArray<- movArray(NumPatches,(lh$Range/2)*NumPatches,'Wrap')
    
    lh$Bmsy<- -999
    
    LengthAtAge<- Length(1:lh$MaxAge) #Calculate length at age vector
    WeightAtAge<- Weight(LengthAtAge,lh$WeightForm) #Calculate weight at age vector
    FecundityAtAge<- Fecundity(WeightAtAge,'Weight') #Calculate Fecundity at age vector
    
    MaturityMode<- as.character(LifeHistory$MaturityMode)
    
    MaturityAtAge<- Maturity(1:lh$MaxAge, MaturityMode) # Calculate maturity at age vector
    
    lh$R0<- lh$B0*(1-exp(-lh$m[1])) #Recruitment must equal deaths at EQ
    
    SpawningPotential<- 10 *sum(MaturityAtAge * FecundityAtAge *exp(-lh$m*(1:lh$MaxAge)))
    
    if (Species=='Blacknose shark')
    {
      SpawningPotential<- 10*SpawningPotential
    }
    
    # lh$R0<- lh$B0*(1-exp(-lh$m[1])) #Recruitment must equal deaths at EQ
    
    lh$R0<- SpawningPotential
    
    pdf(file=paste(FigureFolder,'Life History.pdf',sep=''),width=8,height=6)
    par(mfrow=c(2,2),mar=c(4,4,1,1),oma=c(4,6,2,6))
    plot(LengthAtAge/10,ylab='Length(cm)',type='l',lwd=2,xlab=NA,xaxt='n')
    title(main=Species)
    plot(WeightAtAge*2.2,ylab='Weight(lbs)',type='l',lwd=2,xlab=NA,xaxt='n')
    plot(MaturityAtAge,ylab='Probability Mature',type='l',lwd=2,xlab='Age')
    plot(FecundityAtAge,ylab='# of Eggs',type='l',lwd=2,xlab='Age')
    
    dev.off()
    
    show(SpeciesList[s])
    
    ####### SET UP INITIAL POPULATION ########
    
    EQPopulation<- GrowPopulation(1000,rep(0,NumPatches),'EQ',1,'EQ Run') #Run the population out to unfished equilibrium
    
    lh$CarryingCapacityWeight<- (colSums(EQPopulation$FinalNumAtAge*WeightAtAge)) #Calculate carrying capacity in weight
    
    UnfishedPopulation<- EQPopulation$FinalNumAtAge #Unfished numbers at age
    
    ####### Calculate Reference Points ########
    
    Fmsy<- optimize(log(lh$m[1]),f=FindReferencePoint,Target='FMSY',TargetValue=NA,lower=-10,upper=4) #Find FMSY  
    
    #     SearchSeq<- seq(from=0,to=2,length.out=25)
    #     
    #     ARG<- NULL
    #     for (ff in 1:length(SearchSeq))
    #     {
    #       
    #       ARG[ff]<- FindReferencePoint(log(SearchSeq[ff]),'FMSY',NA)
    #       
    #     }
    #     
    #     pdf(file=paste(FigureFolder,'Production Check.pdf',sep=''))   #Check production model
    #     plot((SearchSeq),-ARG,type='l',ylim=c(min(-ARG),1.1*-Fmsy$objective))
    #     points(exp(Fmsy$minimum) , -Fmsy$objective)
    #     dev.off()
    
    Fmsy$par<- exp(Fmsy$minimum) 
    
    BmsyPopulation<- GrowPopulation(UnfishedPopulation,rep(Fmsy$par,NumPatches),'EQ',0,'Bmsy Run') #Bmsy Population
    
    lh$Nmsy<- (colSums(BmsyPopulation$FinalNumAtAge)) #Nmsy
    
    lh$Bmsy<- colSums(BmsyPopulation$FinalNumAtAge*WeightAtAge) #Bmsy
    
    # DIE<- GrowPopulation(UnfishedPopulation,-log(exp(-Fmsy$par)/1.5),100,1,'DIE')
    
    F50<- optimize(log(2*Fmsy$par),f=FindReferencePoint,Target='BvBmsy',TargetValue=0.5,lower=-10,upper=2) #Find F that results in target B/Bmsy
    
    F50$par<- exp(F50$minimum)
    
    #     Rprof()
    F25<- optimize(log(2*Fmsy$par),f=FindReferencePoint,Target='BvBmsy',TargetValue=0.25,lower=-10,upper=4) #Find F that results in target B/Bmsy
    
    F25$par<- exp(F25$minimum)
    
    B50Population<- GrowPopulation(UnfishedPopulation, rep(F50$par,NumPatches),'EQ',1,'B50 Run')
    #      Rprof()
    
    B25Population<- GrowPopulation(UnfishedPopulation, rep(F25$par,NumPatches),'EQ',1,'B25 Run')
    #      Rprof(NULL)
    #      RProfData<- readProfileData('Rprof.out')
    #      flatProfile(RProfData,byTotal=TRUE)
    ####### RUN MPA SIMULATIONS ########
    
    MPAs<- as.data.frame(matrix(NA,ncol=length(MPANames),nrow=OptTime+1))
    
    colnames(MPAs)<- MPANames
    
    TimeToRun<-dim(MPAs)[1]
    
    EvalTime<- OptTime+1 #Time span to evaluate results on
    RunTime<- 'Custom'
    PropNames<- NULL
    for (m in 1:dim(MPAs)[2])
    {
      PropNames[m]<-MPANames[m]
      #       PropNames[m]<- paste(round(MPAs[m,2:TimeToRun],2),collapse=':')
      
    }
    
    FScenarios<- c(F50$par, F25$par) #load in fishing scenarios
    
    ExperimentResults<-(array(NA,dim=c(dim(MPAs)[2],length(FScenarios),length(DataNames)))) #blank
    
    dimnames(ExperimentResults)<- list(PropNames,c('F50','F25'),DataNames )
    
    
    FlatExperimentResults<- as.data.frame(matrix(NA,nrow=length(FScenarios)*dim(MPAs)[2],ncol=length(LongDataNames)))
    
    colnames(FlatExperimentResults)<- c(LongDataNames)
    
    BaseConditions<- as.data.frame(matrix(NA,nrow= length(FScenarios),ncol=3))
    colnames(BaseConditions)<- c('Yield','Biomass','Numbers')
    
    OptimalConditions<- as.data.frame(matrix(NA,nrow= length(FScenarios),ncol=3))
    colnames(OptimalConditions)<- c('Yield','Biomass','Numbers')
    
    if (RunTime=='Custom')
    {
      TimeToRun<- 75
    }
    
    
    BaseRelativeNumbers<- array(NA,dim=c(dim(MPAs)[2],length(FScenarios),TimeToRun))
    BaseRelativeYield<- array(NA,dim=c(dim(MPAs)[2],length(FScenarios),TimeToRun))
    BaseRelativeBiomass<- array(NA,dim=c(dim(MPAs)[2],length(FScenarios),TimeToRun))
    
    RawNumbers<- array(NA,dim=c(dim(MPAs)[2],length(FScenarios),TimeToRun))
    RawYield<- array(NA,dim=c(dim(MPAs)[2],length(FScenarios),TimeToRun))
    RawBiomass<- array(NA,dim=c(dim(MPAs)[2],length(FScenarios),TimeToRun))
    
    OptRelativeNumbers<- array(NA,dim=c(dim(MPAs)[2],length(FScenarios),TimeToRun))
    OptRelativeYield<- array(NA,dim=c(dim(MPAs)[2],length(FScenarios),TimeToRun))
    OptRelativeBiomass<- array(NA,dim=c(dim(MPAs)[2],length(FScenarios),TimeToRun))
    
    
    FleetSpill<- 1
    
    TotalStorage<- as.data.frame(matrix(NA,nrow<- length(FScenarios)*dim(MPAs)[2]*TimeToRun,ncol=18))
    
    colnames(TotalStorage)<- c('Species','m','f','Frate','FracNTZ','Year','Yield','Biomass','Numbers','SQYield','SQNumbers','SQBiomass','OptNTZYield','OptNTZNumbers','OptNTZBiomass','OptNTZ','YieldBalance','NPB')
    
    c<-0
    
    cc<- 0
    
    for (f in 1:length(FScenarios)) 
    {
      show(paste('F is ',f))
      Patches<- BasePatches
      BasePop<- GrowPopulation(EQPopulation$FinalNumAtAge,FScenarios[f],'EQ',1,paste('FvFmsy is',round(FScenarios[f]/Fmsy$par,2)))
      
      
      StartPop<- BasePop$FinalNumAtAge
      BaseConditions$Yield[f]<- round(BasePop$Performance$Yields[length(BasePop$Performance$Yields)],2)
      BaseConditions$Biomass[f]<- round(sum(WeightAtAge %*% BasePop$FinalNumAtAge),2)
      BaseConditions$Numbers[f]<- round(sum(BasePop$FinalNumAtAge),2)
      
      OptNTZSize<- 	optim(0.2,FindOptimalMPASize,lower=0,upper=0.999,FTemp=FScenarios[f],StartPop=StartPop,FleetSpill=FleetSpill,method="Brent") #You need a better optimization here, gets really stuck with any kind of stochasticity
      
      Patches<- BasePatches
      
      OptMPAPath1<- (optim(c(OptTime,-OptTime),FindMPATrajectory,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f]))    
      
      OptMPAPath2<- (optim(c(-1,5),FindMPATrajectory,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f],method='BFGS'))    
      
      WhichBest<- which(c(OptMPAPath1$value,OptMPAPath2$value)==min(c(OptMPAPath1$value,OptMPAPath2$value)))
      
      if (WhichBest==1)
      {
        OptMPAPath<- (optim(OptMPAPath1$par,FindMPATrajectory,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f],method='BFGS'))    
        
      }
      if (WhichBest==2)
      {
        OptMPAPath<- (optim(OptMPAPath2$par,FindMPATrajectory,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f],method='BFGS'))    
        
      }
      rev(MPAFunction(c(8,8.00001),1:10, 1))
      
      MPAFunction(OptMPAPath1$par,1:10, OptNTZSize$par)
      
      FindMPATrajectory(OptMPAPath$par,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f])
      
      FindMPATrajectory(OptMPAPath$par,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f])
      
      
      #       Patches$MPALocations<- c(1,0) #Create an MPA
      
      
      #       ii=seq(-20,20,length.out=10)
      #       
      #       xx<- rev(ii)*1.01
      #       
      #       MPAMat<- (matrix(NA,nrow=10,ncol=10))
      #       
      #       for (i in 1:10)
      #       {
      #         for (x in 1:10)
      #         {
      #           MPAMat[i,x]<- -FindMPATrajectory(c(ii[i],xx[x]),FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f])
      #         }
      #       }
      #         colnames(MPAMat)<- xx
      #       rownames(MPAMat)<- ii
      
      
      # FindMPATrajectory(c(5,-1.000),FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f])
      
      OptNTZSize$par<- round(OptNTZSize$par,4)
      if (OptNTZSize$convergence>0)
      {
        warning(paste('Failed to converge on run ',f))
      }
      
      Patches<- BasePatches  
      
      #       Patches$MPALocations<- c(1,0) #Create an MPA
      
      AssignNTZ(OptNTZSize$par,ReservePosition)
      
      FTemp<- MoveFleet(FScenarios[f], OptNTZSize$par,FleetSpill,0)
      
      OptPop<- GrowPopulation(StartPop,FTemp,'EQ',1,paste('Opt MPA is',round(OptNTZSize$par,2),'when FvFmsy is',round(FScenarios[f]/Fmsy$par,2)))
      
      OptimalConditions$Yield[f]<-(OptPop$Performance$Yields[length(OptPop$Performance$Yields)])
      
      OptimalConditions$Biomass[f]<-(sum(WeightAtAge %*% OptPop$FinalNumAtAge))
      
      OptimalConditions$Numbers[f]<-(sum(OptPop$FinalNumAtAge))
      
      
      MPAs$StatusQuo<- 0
      
      MPAs$EqNTZ<- c(0,rep(OptNTZSize$par,OptTime))
      
      MPAs$SNTZ<- c(0,seq(min(1,1.5*OptNTZSize$par),OptNTZSize$par,length.out=OptTime))
      
      MPAs$GNTZ<- c(0,seq(min(1,0.5*OptNTZSize$par),OptNTZSize$par,length.out=OptTime))
      
      MPAs$OptNTZ<- c(0, MPAFunction(OptMPAPath$par,1:OptTime, OptNTZSize$par))
      
      
      pdf(file=paste(FigureFolder,'FvFmsy is',round(FScenarios[f]/Fmsy$par,2),' MPAs.pdf'),family=Font,pointsize=12,width=6,height=4)
      par(mar=c(5.1,4.1,4.1,6.1),xpd=T)
      matplot((MPAs),type='l',col=rainbow(1.5*dim(MPAs)[2])[1:dim(MPAs)[2]],lty=1,xlab='Year',ylab='% NTZ',lwd=4,bty='n')
      legend('topright',inset=c(-.3,0),legend=MPANames,col=rainbow(1.5*dim(MPAs)[2])[1:dim(MPAs)[2]],lty=1,cex=.5)
      dev.off()
      
      #       PropNames<- NULL
      #       for (mm in 1:dim(MPAs)[1])
      #       {
      #         PropNames[mm]<- paste(round(MPAs[mm,2:dim(MPAs)[2]],2),collapse=':')
      #       }
      #       
      for (m in 1:dim(MPAs)[2])
      {
        
        show(paste('M is ',m))
        ResultStorage<- as.data.frame(matrix(NA,nrow=dim(MPAs)[1],ncol=3))
        colnames(ResultStorage)<- c('Yield','Biomass','Numbers')
        TempPop<- StartPop
        Patches<- BasePatches
        
        cc<- cc+1
        for (y in 1:TimeToRun)
        {
          c<- c+1
          
          CurrentMPA<- MPAs[y,m]
          if (y>dim(MPAs)[1])
          {
            CurrentMPA<- MPAs[dim(MPAs)[2],m]
            
          }
          
          AssignNTZ(CurrentMPA,ReservePosition)
          
          FTemp<- FScenarios[f] 
          
          FVec<- MoveFleet(FTemp,CurrentMPA,FleetSpill,0)
          
          PassPop<- GrowPopulation(TempPop,FVec,1,0,'TEST') #grow the new population
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
        
        ExperimentResults[m,f,19]<- which((ResultStorage$Yield[2:TimeToRun]-BaseConditions$Yield[f])>=0)[1] #Years to positive yield balance
        
        # ExperimentResults[m,f,19]<- which((ResultStorage$Yield/BaseConditions$Yield[f])>=1)[1]-1 #Years to positive yield balance
        
        
        ExperimentResults[m,f,20]<- which((ResultStorage$Biomass[2:TimeToRun]-BaseConditions$Biomass[f])>=0)[1] #Years to positive biomass balance
        
        
        ExperimentResults[m,f,22]<- Discount(pmax(0,ResultStorage$Yield[2:EvalTime]-BaseConditions$Yield[f]),Fleet$YieldDiscount,EvalTime-1)$NPV # Surplus over the selected time horizon available to pay loans   
        
        
        FinalNPBPositive<- sign(Discount(ResultStorage$Yield-BaseConditions$Yield[f],Fleet$YieldDiscount,length(ResultStorage$Yield))$NPV)>0 
        
        ExperimentResults[m,f,21]<- which((Discount(ResultStorage$Yield[2:TimeToRun]-BaseConditions$Yield[f],Fleet$YieldDiscount,length(ResultStorage$Yield[2:TimeToRun]))$CumValues)>=0)[1] #Years until positive net yield balance          
        
        if (FinalNPBPositive == F)
        {
          ExperimentResults[m,f,21]<- NA #If the final NPB is negative, set years to positive NPB as NA
        }
        
        YieldBalance<- ResultStorage$Yield-BaseConditions$Yield[f]
        
        YieldBalance<- YieldBalance[2:EvalTime]
        
        NegativeYears<- length(YieldBalance[YieldBalance<0])
        
        RequestedLoan<- Discount(pmin(0,YieldBalance),Fleet$YieldDiscount,EvalTime-1)$NPV #Discounted Requested Loan Amount. 
        
        ExperimentResults[m,f,23]<- -RequestedLoan    #Store Discounted Requested Loan Amount
        
        FlatExperimentResults[cc,1]<- SpeciesList[s]
        
        FlatExperimentResults[cc,2:dim(FlatExperimentResults)[2]]<- c(lh$Range,lh$BH.Steepness,m,f,FScenarios[f]/Fmsy$par,OptNTZSize$par,as.numeric(ExperimentResults[m,f,]))
        # 
        #         FlatExperimentResults[cc,2:dim(FlatExperimentResults)[2]]<- data.frame(lh$Range,lh$BH.Steepness,m,f,FScenarios[f]/Fmsy$par,OptNTZSize$par,as.numeric(ExperimentResults[m,f,]))
        #         
        #         colnames(FlatExperimentResults)<- LongDataNames
        #         
      } #Close loop over MPA proposals
      
    } #Close loop over fishing scenarios
    
    
    #######Perform Loan Analysis ######
    
    LoanRates<- seq(.01,2,length.out=20)
    LoanTime<- 10
    
    PlotLayout<- matrix(1:length(FScenarios),nrow=length(FScenarios),ncol=1)
    
    Colors=terrain.colors(2*(dim(MPAs)[2]))[1:(dim(MPAs)[2])]
    
    for (d in 1:dim(FlatExperimentResults)[1])
    {
      
      MaxInterestRate<- NA
      
      if (as.numeric(FlatExperimentResults$RequestedLoan[d])>0 & as.numeric(FlatExperimentResults$RequestedLoan[d])<=as.numeric(FlatExperimentResults$TenYearNPSurplus[d]) )
      {
        
        MaxInterestRate<- optim(-4,FindMaxInterestRate,LoanTime=LoanTime,LoanAmount=as.numeric(FlatExperimentResults$RequestedLoan[d]),Surplus=as.numeric(FlatExperimentResults$TenYearNPSurplus[d]),lower=-10,upper=10,method='Brent')$par
        
      }
      FlatExperimentResults$MaxInterestRate[d]<- 100*exp(MaxInterestRate)
    }
    
    
    ##### Store Things ####
    
    # FlatExperimentResults$MaxInterestRate[FlatExperimentResults$MaxInterestRate< exp(-4)]<- 0
    
    # FlatExperimentResults$MaxInterestRate[FlatExperimentResults$MaxInterestRate> 2]<- Inf
    
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
    
    write.csv(OptRelativeNumbers,file=paste(ResultFolder,'Numbers relative to optimum numbers.csv'))
    
    ResultNames<- dimnames(ExperimentResults)[3]
    ResultNames<- ResultNames[[1]]
    
    FNames<- dimnames(ExperimentResults)[2]
    FNames<- FNames[[1]]
    
  } #Close species list loop
  
  save.image(file=paste(ResultFolder,'Completed Workspace.rdata'))
  
  AllSpeciesExperimentResults$MPAScenario<- as.factor(AllSpeciesExperimentResults$MPAScenario)
  
  levels(AllSpeciesExperimentResults$MPAScenario)<- MPANames
  
  AllSpeciesExperimentResults<- AllSpeciesExperimentResults[AllSpeciesExperimentResults$MPAScenario!='StatusQuo',]
  
  TotalStorage$Year<- as.numeric(TotalStorage$Year)
  
  TotalStorage$NPB<- as.numeric(TotalStorage$NPB)
  
  TotalStorage$m<- as.factor(TotalStorage$m)
  
  TotalStorage$f<- as.factor(TotalStorage$f)
  
  levels(TotalStorage$m)<- MPANames
  
  levels(TotalStorage$f)<- c('F50','F25')
  
  AllSpeciesStorage$Year<- as.numeric(AllSpeciesStorage$Year)
  
  AllSpeciesStorage$NPB<- as.numeric(AllSpeciesStorage$NPB)
  
  AllSpeciesStorage $m<- as.factor(AllSpeciesStorage$m)
  
  AllSpeciesStorage $f<- as.factor(AllSpeciesStorage$f)
  
  levels(AllSpeciesStorage $m)<- MPANames
  
  levels(AllSpeciesStorage $f)<- c('F50','F25')
  
  SubScenarios<- c('EqNTZ','GNTZ','SNTZ','OptNTZ')
  
  Where<-   TotalStorage$m %in% SubScenarios
  
  PlotStorage<- AllSpeciesStorage[Where,]
  
  PlotStorage$m<- as.character(PlotStorage$m)
  
  PlotStorage$m<- as.factor(PlotStorage$m)
  
  PlotStorage$FracNTZ<- as.numeric(PlotStorage$FracNTZ)
  
  PlotStorage$Yield<- as.numeric(PlotStorage$Yield)
  
  PlotStorage$SQYield<- as.numeric(PlotStorage$SQYield)

  
  
  PlotStorage$PositiveYields<- 0
  
  FNames<- levels(TotalStorage$f)
  
  
} #Run Analysis Loop
if (RunAnalysis==0)
{
  load(paste(BatchFolder,'Completed Workspace.rdata'))
}

PlotStorage$YieldBalance<- as.numeric(PlotStorage$YieldBalance)

for (s in 1:length(SpeciesList))
{
  
  for (m in 1:length(SubScenarios))
  {
    
    for (ff in 1:length(FNames))
    {
      
      WhereIsIt1<- (AllSpeciesExperimentResults$Species==SpeciesList[s]& AllSpeciesExperimentResults$MPAScenario==SubScenarios[m] & AllSpeciesExperimentResults$FishingScenario==ff)
      
      WhereIsIt2<- (PlotStorage$Species==SpeciesList[s]& PlotStorage$m==SubScenarios[m] & PlotStorage$f==FNames[ff])
      
      WherePositiveYields<- which(PlotStorage$Species==SpeciesList[s]& PlotStorage$m==SubScenarios[m] & PlotStorage$f==FNames[ff] & PlotStorage$YieldBalance>=0 & PlotStorage$Year>1)[1]
      
      PlotStorage$PositiveYields[WherePositiveYields]<- 1
      
      WhereYields<- which(PlotStorage$Species==SpeciesList[s]& PlotStorage$m==SubScenarios[m] & PlotStorage$f==FNames[ff] & PlotStorage$Year>1 & PlotStorage$Year<=EvalTime)
      
      NegativeYields<-  Discount(pmin(0,PlotStorage$YieldBalance[WhereYields]),Fleet$YieldDiscount,EvalTime-1)$NPV
      
      StatusQuoYields<-  Discount(PlotStorage$SQYield[WhereYields]*as.numeric(PlotStorage$YieldBalance[WhereYields]<=0),Fleet$YieldDiscount,EvalTime-1)$NPV
      
      PriceIncreaseNeeded<- 100*(StatusQuoYields/(StatusQuoYields+NegativeYields)-1) # %increase in prices needed to match status quo profits
      
      YieldBalance<- PlotStorage$YieldBalance[WhereYields]
                  
      RequestedLoan<- Discount(pmin(0,YieldBalance),Fleet$YieldDiscount,EvalTime-1)$NPV #Discounted Requested Loan Amount. 

      AvailableSurplus<- Discount(pmax(0,YieldBalance),Fleet$YieldDiscount,EvalTime-1)$NPV #Discounted Requested Loan Amount. 
            
      if (RequestedLoan!=0)
       {
         
        MaxInterestRate<- optim(-4,FindMaxInterestRate,LoanTime=LoanTime,LoanAmount=-RequestedLoan,Surplus=AvailableSurplus,lower=-10,upper=10,method='Brent')$par
        
        PlotStorage$MaxInterestRate[WhereIsIt2]<- 100*exp(MaxInterestRate)
        
        AllSpeciesExperimentResults$MaxInterestRate[WhereIsIt1]<- round(100*exp(MaxInterestRate),2)
        
       }

      PlotStorage$PriceIncreaseNeeded[WhereIsIt2]<- PriceIncreaseNeeded
      
      AllSpeciesExperimentResults$PriceIncreaseNeeded[WhereIsIt1]<- PriceIncreaseNeeded


    }
  }
  
}


PlotStorage$PosYears<- PlotStorage$Year * PlotStorage$PositiveYields

PlotStorage$PosNPB<- PlotStorage$NPB * PlotStorage$PositiveYields

Cols<-   scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65, h.start = 0, direction = 1)

plotyears=11  #no. years to plot

###Not as slick as xyplot, but gives you full control over what is plotted
for (s in 1:length(SpeciesList))
{
  
  MaxNPB<- max(PlotStorage$NPB[PlotStorage$Year<=plotyears])
  
  FigureFolder<- paste(BatchFolder,SpeciesList[s],'/Figures/',sep='')
  ResultFolder<- paste(BatchFolder,SpeciesList[s],'/Results/',sep='')
  
  pdf(file=paste(FigureFolder,'NPB Trajectory.pdf',sep=''),col='rgb')
  #subset data for each scenario (m)
  m1=subset(PlotStorage,m=="EqNTZ" & Species== SpeciesList[s])
  m2=subset(PlotStorage,m=="GNTZ" & Species == SpeciesList[s])
  m3=subset(PlotStorage,m=="SNTZ"& Species == SpeciesList[s])
  m4=subset(PlotStorage,m=="OptNTZ"& Species == SpeciesList[s])
  
  # m4G=subset(PlotStorage,m=="OptPath"& Species == SpeciesList[2])
  
  y1_FirstYields<-(as.numeric(PlotStorage$YieldBalance[PlotStorage$f=="F50" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s]]))
  
  y2_FirstYields<-(as.numeric(PlotStorage$YieldBalance[PlotStorage$f=="F25" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s]]))
  
  
  MaxYields<- max(y1_FirstYields, y2_FirstYields)
  
  y1_FirstYields<- y1_FirstYields/MaxYields
  y2_FirstYields<- y2_FirstYields/MaxYields
  
  #get points where yield equals status quo
  y1=PlotStorage[PlotStorage$f=="F50" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s],c(6,18) ]
  y2=PlotStorage[PlotStorage$f=="F25" & PlotStorage$PositiveYields==1 & PlotStorage$Species== SpeciesList[s],c(6,18)]
  
  par(mfrow=c(2,1),mar=c(0,0,0,0),oma=c(4,6,2,6))
  
  plot(m1$Year[as.vector(m1$f)=="F50"],m1$NPB[as.vector(m1$f)=="F50"],
       type="l",ylab="",xlab="",ylim=c(min(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears]),1.25*max(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears])),
       xlim=c(1,plotyears),col="dodgerblue",lwd=4,las=1,xaxt="n" )
  mtext("Moderate Overfishing",line=-1)
  text(1,max(PlotStorage$NPB)*0.99,"A.")
  lines( m2$Year[as.vector(m2$f)=="F50"],m2$NPB[as.vector(m2$f)=="F50"],col=2,lwd=2,lty=3)
  lines( m3$Year[as.vector(m3$f)=="F50"],m3$NPB[as.vector(m3$f)=="F50"],col="medium orchid",lwd=2,lty=4)
  lines( m4$Year[as.vector(m4$f)=="F50"],m4$NPB[as.vector(m4$f)=="F50"],col="darkgreen",lwd=2,lty=5)
  abline(h=0,lty=2)
  points(y1[,1],y1[,2],col=c("dodgerblue","medium orchid",2,"darkgreen"),pch=16,cex=2)
  
  plot(m1$Year[as.vector(m1$f)=="F25"],m1$NPB[as.vector(m1$f)=="F25"],
       type="l",ylab="",xlab="Year",ylim=c(min(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears ]),1.25*max(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears])),
       xlim=c(1,plotyears),col="dodgerblue",lwd=4 ,las=1)
  mtext("Heavy Overfishing",line=-1)
  text(1,max(PlotStorage$NPB)*0.99,"B.")
  lines( m2$Year[as.vector(m2$f)=="F25"],m2$NPB[as.vector(m2$f)=="F25"],col=2,lwd=2,lty=3)
  lines( m3$Year[as.vector(m3$f)=="F25"],m3$NPB[as.vector(m3$f)=="F25"],col="medium orchid",lwd=2,lty=4)
  lines( m4$Year[as.vector(m4$f)=="F25"],m4$NPB[as.vector(m4$f)=="F25"],col="darkgreen",lwd=2,lty=5)
  abline(h=0,lty=2)         
  points(y2[,1],y2[,2],col=c("dodgerblue","medium orchid",2,"darkgreen"),pch=16, cex=2)
  
  par(new=TRUE, mfrow=c(1,1), oma=c(0,0,0,0))
  
  plot(1,
       xlim=c(0,1),
       ylim=c(0,1),
       xaxs="i",
       yaxs="i",
       type="n",
       axes=FALSE)
  
  MaxCircle<- round((MaxYields),0) 
  
  Difference<- 10-MaxCircle%%10
  
  MaxCircle <- MaxCircle + Difference
  
  CircleSize<- ((seq(from=MaxCircle/10,to=MaxCircle,length.out=5)))
  
  legend(0.85,0.9, legend= SubScenarios, lty=c(1,3,4,5),
         col=c("dodgerblue",2,"medium orchid","darkgreen"),lwd=2,bty="n",cex=0.8)
  
  # legend(0.85,0.75, legend= CircleSize,pch=16,pt.cex=(CircleSize),col=c("black"),title='Yield Surplus',bty='n',y.intersp=1.5)
  
  
  text(.5,.05,"Year",font=2)
  
  text(.03,.5,"Net Present Balance",srt=90,font=2)
  
  dev.off()
  
  pdf(file=paste(FigureFolder,'FracNTZ Trajectory.pdf',sep=''))
  print(xyplot(FracNTZ ~ Year | f,group=m,data=PlotStorage,type='l',subset=(Year<=11 & PlotStorage$Species==SpeciesList[s]),auto.key=list(space='right',points=F,lines=T),layout=c(1,2),index.cond=list(c(2,1)),lwd=4,drop.unused.levels=T,panel=function(x,y,...)
  {panel.xyplot(x,y,...)
   panel.abline(h=0,lty=2)
  },))
  dev.off()
  
  pdf(file=paste(FigureFolder,'Heavy Fishing NPB Trajectory.pdf',sep=''),col='rgb')
  #subset data for each scenario (m)
  m1=subset(PlotStorage,m=="EqNTZ" & Species== SpeciesList[s])
  m2=subset(PlotStorage,m=="GNTZ" & Species == SpeciesList[s])
  m3=subset(PlotStorage,m=="SNTZ"& Species == SpeciesList[s])
  m4=subset(PlotStorage,m=="OptNTZ"& Species == SpeciesList[s])
  
  # m4G=subset(PlotStorage,m=="OptPath"& Species == SpeciesList[2])
  
  y1_FirstYields<-(as.numeric(PlotStorage$YieldBalance[PlotStorage$f=="F50" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s]]))
  
  y2_FirstYields<-(as.numeric(PlotStorage$YieldBalance[PlotStorage$f=="F25" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s]]))
  
  
  MaxYields<- max(y1_FirstYields, y2_FirstYields)
  
  y1_FirstYields<- y1_FirstYields/MaxYields
  y2_FirstYields<- y2_FirstYields/MaxYields
  
  #get points where yield equals status quo
  y1=PlotStorage[PlotStorage$f=="F50" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s],c(6,18) ]
  y2=PlotStorage[PlotStorage$f=="F25" & PlotStorage$PositiveYields==1 & PlotStorage$Species== SpeciesList[s],c(6,18)]
  
  par(mar=c(0,0,0,0),oma=c(4,6,2,6))
  
  plot(m1$Year[as.vector(m1$f)=="F25"],m1$NPB[as.vector(m1$f)=="F25"],
       type="l",ylab="",xlab="Year",ylim=c(min(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears ]),1.25*max(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears])),
       xlim=c(1,plotyears),col="dodgerblue",lwd=4 ,las=1)
  mtext("Heavy Overfishing",line=-1)
  text(1,max(PlotStorage$NPB)*0.99,"B.")
  lines( m2$Year[as.vector(m2$f)=="F25"],m2$NPB[as.vector(m2$f)=="F25"],col=2,lwd=2,lty=3)
  lines( m3$Year[as.vector(m3$f)=="F25"],m3$NPB[as.vector(m3$f)=="F25"],col="medium orchid",lwd=2,lty=4)
  lines( m4$Year[as.vector(m4$f)=="F25"],m4$NPB[as.vector(m4$f)=="F25"],col="darkgreen",lwd=2,lty=5)
  abline(h=0,lty=2)         
  points(y2[,1],y2[,2],col=c("dodgerblue","medium orchid",2,"darkgreen"),pch=16, cex=2)
  
  par(new=TRUE, mfrow=c(1,1), oma=c(0,0,0,0))
  
  plot(1,
       xlim=c(0,1),
       ylim=c(0,1),
       xaxs="i",
       yaxs="i",
       type="n",
       axes=FALSE)
  
  MaxCircle<- round((MaxYields),0) 
  
  Difference<- 10-MaxCircle%%10
  
  MaxCircle <- MaxCircle + Difference
  
  CircleSize<- ((seq(from=MaxCircle/10,to=MaxCircle,length.out=5)))
  
  legend(0.85,0.9, legend= SubScenarios, lty=c(1,3,4,5),
         col=c("dodgerblue",2,"medium orchid","darkgreen"),lwd=2,bty="n",cex=0.8)
  
  # legend(0.85,0.75, legend= CircleSize,pch=16,pt.cex=(CircleSize),col=c("black"),title='Yield Surplus',bty='n',y.intersp=1.5)
  
  
  text(.5,.05,"Year",font=2)
  
  text(.03,.5,"Net Present Balance",srt=90,font=2)
  
  dev.off()
  
  pdf(file=paste(FigureFolder,'Heavy Overfishing FracNTZ Trajectory.pdf',sep=''))
  print(xyplot(FracNTZ ~ Year,group=m,data=PlotStorage,type='l',subset=(Year<=11 & PlotStorage$Species==SpeciesList[s] & f=='F25'),auto.key=list(space='right',points=F,lines=T),lwd=4,drop.unused.levels=T,panel=function(x,y,...)
  {panel.xyplot(x,y,...)
   panel.abline(h=0,lty=2)
  },))
  dev.off()
  
} #Close Species Loop


PaperTable<- ddply(AllSpeciesExperimentResults,c('Species','FishingScenario','MPAScenario'),summarize,
                   PositiveNPBYear=YearsToBalanceRecovery,PositivePBYear=YearsToBalanceRecovery,MaxInterestRate=MaxInterestRate,PriceBoostNeeded=PriceIncreaseNeeded)

PaperTable$FishingScenario<- as.factor(PaperTable$FishingScenario)

levels(PaperTable$FishingScenario)<- c('Moderate Overfishing','Heavy Overfishing')

write.csv(file=paste(BatchFolder,'Summary Table for Paper.csv',sep=''),PaperTable)

write.csv(file=paste(BatchFolder,'All Species Experiment Results.csv',sep=''),AllSpeciesExperimentResults)

AllSpeciesStorage$m<- as.factor(AllSpeciesStorage$m)
levels(AllSpeciesStorage$m)<- MPANames

write.csv(file=paste(BatchFolder,'All Species All Results.csv',sep=''), AllSpeciesStorage)

write.csv(file=paste(BatchFolder,'All Species Paper Results.csv',sep=''), PlotStorage)


save.image(file=paste(BatchFolder,'Completed Workspace.rdata'))

