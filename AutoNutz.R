rm(list=ls())
library(lattice)
require(stats)
library(proftools)
#Is this working?
####### Shrinking MPAs #########
#Wilson, Doughtery, Ovando et al something something

####### READ IN CONTROL FILE ########

# setwd('/Users/danovando/Dropbox/Shrinking NTZ')

RunAnalysis<- 1

DataNames<- c('NPV of Yield','NPV of Biomass','Mean Yield','Mean Biomass','Mean Numbers','Yield Instability','Mean Changes in Yield','Percent Years with Profit Gains','Percent Years with Numbers Gains','Mean Percent Change in Yield','Mean Percent Change in Numbers','Percent Years With Numbers and Yield Gains','FiveyearYieldBalance','FiveyearBiomassBalance','FiveyearNPVBalance','TenyearYieldBalance','TenyearBiomassBalance','TenyearNPVBalance','YearsToYieldRecovery','YearsToBioRecovery','YearsToBalanceRecovery','TenYearNPSurplus','RequestedLoan','MaxInterestRate')

LongDataNames<- c('Species','Movement','Steepness','MPAScenario','FishingScenario','FvFmsy','OptNTZ','NPV of Yield','NPV of Biomass','Mean Yield','Mean Biomass','Mean Numbers','Yield Instability','Mean Changes in Yield','Percent Years with Profit Gains','Percent Years with Numbers Gains','Mean Percent Change in Yield','Mean Percent Change in Numbers','Percent Years With Numbers and Yield Gains','FiveyearYieldBalance','FiveyearBiomassBalance','FiveyearNPVBalance','TenyearYieldBalance','TenyearBiomassBalance','TenyearNPVBalance','YearsToYieldRecovery','YearsToBioRecovery','YearsToBalanceRecovery','TenYearNPSurplus','RequestedLoan','MaxInterestRate')

AllSpeciesStorage<- as.data.frame(matrix(NA,nrow=0,ncol=18))

colnames(AllSpeciesStorage)<- c('Species','m','f','Frate','FracNTZ','Year','Yield','Biomass','Numbers','SQYield','SQNumbers','SQBiomass','OptNTZYield','OptNTZNumbers','OptNTZBiomass','OptNTZ','YieldBalance','NPB')

AllSpeciesExperimentResults<- as.data.frame(matrix(NA,nrow=0,ncol=length(LongDataNames)))

colnames(AllSpeciesExperimentResults)<- c(LongDataNames)

BatchFolder<- 'Results/TEST/'
dir.create(paste(BatchFolder,sep=''))

MPANames<- c('Status Quo','EqNTZ','Rotate','SNTZ','Basic','GNTZ','OptNTZ')

LifeHistories<- read.csv('Inputs/Life History Inputs.csv')

# LifeHistories<- LifeHistories[7,]

LifeColumns<- colnames(LifeHistories)

LifeVars<- c('Range','MaxAge','m','k','Linf','t0','AgeMa50','LengthMa50','MaturityMode','BH.Steepness','wa','wb','WeightForm','fa','fb')

LifeHistories$SciName<- as.character(levels(LifeHistories$SciName))[LifeHistories$SciName]

LifeHistories$CommName<- as.character(levels(LifeHistories$CommName))[LifeHistories$CommName]

LifeHistories$WeightForm<- as.character(levels(LifeHistories$WeightForm))[LifeHistories$WeightForm]

# LifeHistories<- LifeHistories[1,]

SpeciesList<- LifeHistories$CommName

if (RunAnalysis==1)
{
  
  for (s in 1:length(SpeciesList))
  {
    Species<- SpeciesList[s] #species file to load
    
    StoreRun<- paste(BatchFolder,Species,sep='') #1 if you want to create a seeded folder to store results, 0 if you want it in the generic working folder
    
    source('SHRINKNTZ_CONTROLFILE.R')
    # Fleet$YieldDiscount<- 0
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
    
    lh$MoveType<- 'Simple-DO'
    
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
    
    # lh$MoveType<- 'None'
    show(SpeciesList[s])
    
    ####### SET UP INITIAL POPULATION ########
    
    EQPopulation<- GrowPopulation(1000,rep(0,NumPatches),'EQ',1,'EQ Run') #Run the population out to unfished equilibrium
    
    lh$CarryingCapacity<- (colSums(EQPopulation$FinalNumAtAge)) #Calculate carrying capacity in numbers
    
    lh$CarryingCapacityWeight<- (colSums(EQPopulation$FinalNumAtAge*WeightAtAge)) #Calculate carrying capacity in weight
    
    UnfishedPopulation<- EQPopulation$FinalNumAtAge #Unfished numbers at age
    
    ####### Calculate Reference Points ########
    
    Fmsy<- optimize(log(lh$m[1]),f=FindReferencePoint,Target='FMSY',TargetValue=NA,lower=-10,upper=4) #Find FMSY  
    
    SearchSeq<- seq(from=0,to=5,length.out=25)
    
    ARG<- NULL
    for (ff in 1:length(SearchSeq))
    {
      
      ARG[ff]<- FindReferencePoint(log(SearchSeq[ff]),'FMSY',NA)
      
    }
    
    pdf(file=paste(FigureFolder,'Production Check.pdf',sep=''))   #Check production model
    plot((SearchSeq),-ARG,type='l',ylim=c(min(-ARG),1.1*-Fmsy$objective))
    points(exp(Fmsy$minimum) , -Fmsy$objective)
    dev.off()
    
    Fmsy$par<- exp(Fmsy$minimum) 
    
    BmsyPopulation<- GrowPopulation(UnfishedPopulation,rep(Fmsy$par,NumPatches),'EQ',0,'Bmsy Run') #Bmsy Population
    
    lh$Nmsy<- (colSums(BmsyPopulation$FinalNumAtAge)) #Nmsy
    
    lh$Bmsy<- colSums(BmsyPopulation$FinalNumAtAge*WeightAtAge) #Bmsy
    
    BmsyPopulation<- GrowPopulation(UnfishedPopulation,rep(Fmsy$par,NumPatches),'EQ',1,'Bmsy Run') #Rerun Bmsy Population
    
    # DIE<- GrowPopulation(UnfishedPopulation,-log(exp(-Fmsy$par)/1.5),100,1,'DIE')
    
    F50<- optimize(log(2*Fmsy$par),f=FindReferencePoint,Target='BvBmsy',TargetValue=0.5,lower=-10,upper=4) #Find F that results in target B/Bmsy
    
    F50$par<- exp(F50$minimum)
    
    
    #     Rprof()
    F25<- optimize(log(2*Fmsy$par),f=FindReferencePoint,Target='BvBmsy',TargetValue=0.25,lower=-10,upper=4) #Find F that results in target B/Bmsy
    
    
    F25$par<- exp(F25$minimum)
    
    B50Population<- GrowPopulation(UnfishedPopulation, rep(F50$par,NumPatches),'EQ',1,'B50 Run')
    # Rprof()
    
    B25Population<- GrowPopulation(UnfishedPopulation, rep(F25$par,NumPatches),'EQ',1,'B25 Run')
    # Rprof(NULL)
    # RProfData<- readProfileData('Rprof.out')
    # flatProfile(RProfData,byTotal=TRUE)
    ####### RUN MPA SIMULATIONS ########
    
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
    
    FScenarios<- c(F50$par, F25$par) #load in fishing scenarios
    
    ExperimentResults<-(array(NA,dim=c(dim(MPAs)[1],length(FScenarios),length(DataNames)))) #blank
    
    dimnames(ExperimentResults)<- list(PropNames,c('F50','F25'),DataNames )
    
    
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
      
      Patches<- BasePatches
      
      Patches$MPALocations<- c(1,0) #Create an MPA
      
      OptMPAPath1<- (optim(c(10,-10),FindMPATrajectory,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f],method='BFGS'))    
      
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
      MPAFunction(c(10,-1),1:10, OptNTZSize$par)
      
      MPAFunction(OptMPAPath2$par,1:10, OptNTZSize$par)
      
      FindMPATrajectory(OptMPAPath$par,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f])
      
      FindMPATrajectory(OptMPAPath$par,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,OptMode='Function',BaseYields=BaseConditions$Yield[f])
      
      
      
      
      Patches$MPALocations<- c(1,0) #Create an MPA
      
      
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
      
      MPAs[7,]<- c(0, MPAFunction(OptMPAPath$par,1:10, OptNTZSize$par))
      
      
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
        
        YieldBalance<- YieldBalance[2:TimeToRun]
        
        NegativeYears<- length(YieldBalance[YieldBalance<0])
        
        NegativeYields<- Discount(YieldBalance[YieldBalance<0],Fleet$YieldDiscount,NegativeYears)$NPV #Discounted Requested Loan Amount. 
        
        ExperimentResults[m,f,23]<- -NegativeYields    #Store Discounted Requested Loan Amount
        
        FlatExperimentResults[cc,1]<- SpeciesList[s]
        
        FlatExperimentResults[cc,2:dim(FlatExperimentResults)[2]]<- c(lh$Range,lh$BH.Steepness,m,f,FScenarios[f]/Fmsy$par,OptNTZSize$par,as.numeric(ExperimentResults[m,f,]))
# 
#         FlatExperimentResults[cc,2:dim(FlatExperimentResults)[2]]<- data.frame(lh$Range,lh$BH.Steepness,m,f,FScenarios[f]/Fmsy$par,OptNTZSize$par,as.numeric(ExperimentResults[m,f,]))
#         
#         colnames(FlatExperimentResults)<- LongDataNames
#         
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
    
    Comparisons<- c(1,2,3,4,3,5,4,6,3,6,4,10,3,11)
    TradeOffs<- matrix(Comparisons,ncol=2,nrow=length(Comparisons)/2,byrow=T)
    
    
    ### Make Time Series Plots ###
    
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
  
save.image(file=paste(ResultFolder,'Completed Workspace.rdata'))




  AllSpeciesExperimentResults$MPAScenario<- as.factor(AllSpeciesExperimentResults$MPAScenario)
  
  levels(AllSpeciesExperimentResults$MPAScenario)<- MPANames
  
  AllSpeciesExperimentResults<- AllSpeciesExperimentResults[AllSpeciesExperimentResults$MPAScenario!='Status Quo',]
  
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
        
        WhereNegativeYields<- which(PlotStorage$Species==SpeciesList[s]& PlotStorage$m==SubScenarios[m] & PlotStorage$f==FNames[ff] & PlotStorage$YieldBalance<0 & PlotStorage$Year>1)
        
        Yields<-  Discount(PlotStorage$Yield[WhereNegativeYields],Fleet$YieldDiscount,length(WhereNegativeYields))$NPV
        
        StatusQuoYields<-  Discount(PlotStorage$SQYield[WhereNegativeYields],Fleet$YieldDiscount,length(WhereNegativeYields))$NPV
        
        PriceIncreaseNeeded<- 100*(StatusQuoYields/Yields-1) # %increase in prices needed to match status quo profits
        
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
  
} #Close Species Loop


sum(is.na(AllSpeciesExperimentResults$MaxInterestRate))/dim(AllSpeciesExperimentResults)[1]

write.csv(file=paste(BatchFolder,'All Species Experiment Results.csv',sep=''),AllSpeciesExperimentResults)

AllSpeciesStorage$m<- as.factor(AllSpeciesStorage$m)
levels(AllSpeciesStorage$m)<- MPANames

write.csv(file=paste(BatchFolder,'All Species All Results.csv',sep=''), AllSpeciesStorage)

write.csv(file=paste(BatchFolder,'All Species Paper Results.csv',sep=''), PlotStorage)


save.image(file=paste(BatchFolder,'Completed Workspace.rdata'))

