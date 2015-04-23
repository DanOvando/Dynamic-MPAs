rm(list=ls())
library(lattice)
require(stats)
library(proftools)
library(plyr)
library(grid)
library(gridExtra)
library(optimx)
library(ggplot2)
library(gridExtra)
library(plyr)
library(dplyr)
library(dfoptim)
library(tidyr)
library(broom)

####### Shrinking MPAs #########
#Wilson, Doughtery, Ovando et al something something

####### READ IN CONTROL FILE ########

PopTolerance<- .1 #the cutoff for identifying point where population stops changing
InitialPopulation<- 1000 #Seed for initial population 
CollapseThreshold<- 0.1
LookAtLengths<- 0
ReservePosition<- 'Center'
OptTime<- 65 #Time Horizon to optimize over
Alpha<- 0.5

TimeToRun<- OptTime+1

# setwd('/Users/danovando/Dropbox/Shrinking NTZ')
BatchFolder<- 'Results/3.0/'

RunAnalysis<- 1

DataNames<- c('NPV of Yield','NPV of Biomass','Mean Yield','Mean Biomass','Mean Numbers','Yield Instability','Mean Changes in Yield','Percent Years with Profit Gains','Percent Years with Numbers Gains','Mean Percent Change in Yield','Mean Percent Change in Numbers','Percent Years With Numbers and Yield Gains','FiveyearYieldBalance','FiveyearBiomassBalance','FiveyearNPVBalance','TenyearYieldBalance','TenyearBiomassBalance','TenyearNPVBalance','YearsToYieldRecovery','YearsToBioRecovery','YearsToBalanceRecovery','TenYearNPSurplus','RequestedLoan','MaxInterestRate')

LongDataNames<- c('Species','Movement','Steepness','MPAScenario','FishingScenario','FvFmsy','OptNTZ','NPV of Yield','NPV of Biomass','Mean Yield','Mean Biomass','Mean Numbers','Yield Instability','Mean Changes in Yield','Percent Years with Profit Gains','Percent Years with Numbers Gains','Mean Percent Change in Yield','Mean Percent Change in Numbers','Percent Years With Numbers and Yield Gains','FiveyearYieldBalance','FiveyearBiomassBalance','FiveyearNPVBalance','TenyearYieldBalance','TenyearBiomassBalance','TenyearNPVBalance','YearsToYieldRecovery','YearsToBioRecovery','YearsToBalanceRecovery','TenYearNPSurplus','RequestedLoan','MaxInterestRate')

dir.create(paste(BatchFolder,sep=''))

# MPANames<- c('Status Quo','EqNTZ','Rotate','SNTZ','Basic','GNTZ','OptNTZ')

MPANames<- c('StatusQuo','EqNTZ','SNTZ','GNTZ','OptNTZ','CatchShareEqNTZ')

LifeHistories<- read.csv('Inputs/Life History Inputs.csv',stringsAsFactors=F)

# LifeHistories<- LifeHistories[7,]

LifeColumns<- colnames(LifeHistories)

LifeVars<- c('Range','MaxAge','m','k','Linf','t0','AgeMa50','LengthMa50','MaturityMode','BH.Steepness','wa','wb','WeightForm','fa','fb','DDForm')

SpeciesList<- LifeHistories$CommName

SystemBmsyStorage<- as.data.frame(matrix(NA,nrow=length(SpeciesList),ncol=2))

colnames(SystemBmsyStorage)<- c('Species','Bmsy')

SystemBmsyStorage$Species<- as.character(SystemBmsyStorage$Species)

SpeciesList<- SpeciesList[1:2]

NumFs<- 1

TotalStorage<- as.data.frame(matrix(NA,nrow= NumFs*length(MPANames)*TimeToRun*length(SpeciesList),ncol=16))

colnames(TotalStorage)<- c('Species','m','f','Frate','FracNTZ','Year','Yield','Biomass','Numbers','SQYield','SQNumbers','SQBiomass','OptNTZYield','OptNTZNumbers','OptNTZBiomass','OptNTZ')

c<- 0

if (RunAnalysis==1)
{
  
  for (s in 1:length(SpeciesList))
  {
    Species<- SpeciesList[s] #species file to load
    
    StoreRun<- paste(BatchFolder,Species,sep='') #1 if you want to create a seeded folder to store results, 0 if you want it in the generic working folder
    
    source('BetterNutzControlfile.R')
    
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
    
    lh$MovementArray<- movArray(NumPatches,((lh$Range)*NumPatches)/5,'Wrap')
    
    lh$Bmsy<- -999
    
    LengthAtAge<- Length(1:lh$MaxAge) #Calculate length at age vector
    
    WeightAtAge<- Weight(LengthAtAge,lh$WeightForm) #Calculate weight at age vector
    
    FecundityAtAge<- Fecundity(WeightAtAge,'Weight') #Calculate Fecundity at age vector
    
    MaturityMode<- as.character(LifeHistory$MaturityMode)
    
    MaturityAtAge<- Maturity(1:lh$MaxAge, MaturityMode) # Calculate maturity at age vector
    
    if (grepl('Shark',Species))
    {
      lh$R0<- 100*Patches$HabQuality
    }
    
    
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
    #     Rprof(tmp <- tempfile(),line.profiling=T)
    
    EQPopulation<- GrowPopulation(1000,rep(0,NumPatches),'EQ',0,'EQ Run') #Run the population out to unfished equilibrium
    
    #     Rprof() 
    #     summaryRprof(tmp)
    #     unlink(tmp)
    
    lh$CarryingCapacityWeight<- (colSums(EQPopulation$FinalNumAtAge*WeightAtAge)) #Calculate carrying capacity in weight
    
    UnfishedPopulation<- EQPopulation$FinalNumAtAge #Unfished numbers at age
    
    ####### Calculate Reference Points ########
    
    Fmsy<- optimize(log(lh$m[1]),f=FindReferencePoint,Target='FMSY',TargetValue=NA,lower=-10,upper=4) #Find FMSY  
    
    Fmsy$par<- exp(Fmsy$minimum) 
    
    BmsyPopulation<- GrowPopulation(UnfishedPopulation,rep(Fmsy$par,NumPatches),'EQ',0,'Bmsy Run') #Bmsy Population
    
    lh$Nmsy<- (colSums(BmsyPopulation$FinalNumAtAge)) #Nmsy
    
    lh$Bmsy<- colSums(BmsyPopulation$FinalNumAtAge*WeightAtAge) #Bmsy
    
    SystemBmsyStorage[s,]<- data.frame(Species,sum(lh$Bmsy),stringsAsFactors=F)
    
    MsyFishing<- GrowPopulation(BmsyPopulation$FinalNumAtAge,rep(Fmsy$par,NumPatches),OptTime,0,'MSY Run') #Bmsy Population
    
    Fleet$MSY_NPV<- MsyFishing$Performance$DiscYields$NPV
    
    F25<- optimize(log(2*Fmsy$par),f=FindReferencePoint,Target='BvBmsy',TargetValue=0.25,lower=-10,upper=4) #Find F that results in target B/Bmsy
    
    F25$par<- exp(F25$minimum)
    
    B25Population<- GrowPopulation(UnfishedPopulation, rep(F25$par,NumPatches),'EQ',0,'B25 Run')
    
    ####### RUN MPA SIMULATIONS ########
    
    MPAs<- as.data.frame(matrix(NA,ncol=length(MPANames),nrow=OptTime+1))
    
    colnames(MPAs)<- MPANames
    
#     TimeToRun<-dim(MPAs)[1]
    
    EvalTime<- OptTime #Time span to evaluate results on
    RunTime<- 'Custom'
    PropNames<- NULL
    for (m in 1:dim(MPAs)[2])
    {
      PropNames[m]<-MPANames[m]      
    }
    
    FScenarios<- c(F25$par) #load in fishing scenarios
    
    BaseConditions<- as.data.frame(matrix(NA,nrow= length(FScenarios),ncol=3))
    colnames(BaseConditions)<- c('Yield','Biomass','Numbers')
    
    OptimalConditions<- as.data.frame(matrix(NA,nrow= length(FScenarios),ncol=3))
    colnames(OptimalConditions)<- c('Yield','Biomass','Numbers')
        
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
      
      wtf<- seq(0,1,length.out=20)
      
      arg=ldply(lapply(wtf,FindOptimalMPASize,FTemp=FScenarios[f],StartPop=StartPop,FleetSpill=FleetSpill))
      
      mguess<- wtf[which(arg==min(arg))[1]]
      
      OptNTZSize<- 	optim(mguess,f=FindOptimalMPASize,lower=0,upper=0.999,FTemp=FScenarios[f],StartPop=StartPop,FleetSpill=FleetSpill) #You need a better optimization here, gets really stuck with any kind of stochasticity
      
      OptNTZSize<-   optim(OptNTZSize$par,f=FindOptimalMPASize,lower=0,upper=0.999,FTemp=FScenarios[f],StartPop=StartPop,FleetSpill=FleetSpill) #You need a better optimization here, gets really stuck with any kind of stochasticity
      
      mcheck<- data.frame(wtf,-(arg))
      
      colnames(mcheck)<- c('MPASize','NPB')
      
      pdf(file=paste(FigureFolder,'MPA Check.pdf',sep=''))
      print(ggplot(data=mcheck,aes(MPASize,NPB))+geom_point()+geom_vline(xintercept=OptNTZSize$par))
      dev.off()
      
      Patches<- BasePatches
      
      Int=seq(0,1,length.out=20)   
      
      Flip=seq(.001,1,length.out=20)   
      
      ObjMat<- matrix(NA,nrow=length(Int)*length(Flip),ncol=3)
      
      l<- 0
      for (i in 1:length(Int))
      {
        for (o in 1:length(Flip))
        {
          l=l+1
          show(l)
          
          ObjMat[l,3]<- -FindMPATrajectory(c(Int[i],Flip[o]),Mode='LinearPlus',EvalTime=OptTime,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,Alpha=1,OptMode='Function',BaseYields=BaseConditions$Yield[f],GrowMode='Shrink')
          
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
      
      #       OptMPAPath<- (nlminb(log(jitter(as.numeric(BestGuess))),FindMPATrajectory,lower=log(c(0,.001)),upper=log(c(1,0.6)),
      #                            control=list(trace=T,step.min=.01,step.max=.1),
      #                            Mode='LinearPlus',EvalTime=OptTime,FTemp= FScenarios[f],
      #                            TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,
      #                            Alpha=1,GrowMode='Grow',OptMode='Function',BaseYields=BaseConditions$Yield[f]))    
      #       
      #       hchange=function(par,lower,upper,dist,round=TRUE,...)
      #       { D=length(par) # dimension
      #         show(par)
      #         step=dist(D,...) # slight step
      #         if(round) step=round(step)
      #         par1=par+step
      #         # return par1 within [lower,upper]:
      #         return(ifelse(par1<lower,lower,ifelse(par1>upper,upper,par1)))
      #       }
      #       
      #       
      #       bchange=function(par) # binary change
      #       { 
      #         D=length(par)
      #         hchange(par,lower=c(0,.001),upper=c(1,0.6),rnorm,mean=0,sd=1)
      #       }
      #         
      #       C=list(maxit=10000,temp=1000,trace=TRUE,REPORT=1)
      #       
      #       s=(optim(BestGuess,f=FindMPATrajectory,gr=bchange,method="SANN",control=C,Mode='LinearPlus',EvalTime=OptTime,FTemp= FScenarios[f],
      #               TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,
      #               Alpha=1,GrowMode='Grow',OptMode='Function',BaseYields=BaseConditions$Yield[f]))
      #       
      #       s=(optim(BestGuess,f=FindMPATrajectory,method='SANN',control=C,Mode='LinearPlus',EvalTime=OptTime,FTemp= FScenarios[f],
      #                TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,
      #                Alpha=1,GrowMode='Grow',OptMode='Function',BaseYields=BaseConditions$Yield[f]))
      #       
      #       
      #       
      #       TestGrow<- (nlminb(c(0.3,0.14,0.01),FindMPATrajectory,lower=c(0,0,0),upper=c(8,200,1),Mode='Logit Grow',EvalTime=OptTime,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,Alpha=1,GrowMode='Grow',OptMode='Function',BaseYields=BaseConditions$Yield[f]))    
      #       
      #       TestShrink1<- (nlminb(c(.2,.8,.5),FindMPATrajectory,lower=c(0,0,0),upper=c(8,10,1),Mode='Logit Shrink',EvalTime=OptTime,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,Alpha=1,GrowMode='Grow',OptMode='Function',BaseYields=BaseConditions$Yield[f]))    
      #       
      #       TestShrink2<- (nlminb(c(.1,.1,.1),FindMPATrajectory,lower=c(0,0,0),upper=c(8,10,1),Mode='Logit Shrink',EvalTime=OptTime,FTemp= FScenarios[f],TimeFrame=OptTime,FleetSpill=FleetSpill,StartPop=StartPop,OptSize=OptNTZSize$par,Alpha=1,GrowMode='Grow',OptMode='Function',BaseYields=BaseConditions$Yield[f]))    
      #       
      
      #       
      #       quartz()
      #       plot(MPAFunction(TestGrow$par,0:OptTime, OptNTZSize$par,'Logit Grow',EvalTime),ylim=c(0,1)) 
      #       lines(MPAFunction(TestShrink$par,0:OptTime, OptNTZSize$par,'Logit Shrink',EvalTime)) 
      #       
      
      OptNTZSize$par<- round(OptNTZSize$par,4)
      if (OptNTZSize$convergence>0)
      {
        warning(paste('Failed to converge on run ',f))
      }
      
      Patches<- BasePatches  
      
      AssignNTZ(OptNTZSize$par,ReservePosition)
      
      FTemp<- MoveFleet(FScenarios[f], OptNTZSize$par,FleetSpill,0)
      
      OptPop<- GrowPopulation(StartPop,FTemp,'EQ',0,paste('Opt MPA is',round(OptNTZSize$par,2),'when FvFmsy is',round(FScenarios[f]/Fmsy$par,2)))
      
      OptimalConditions$Yield[f]<-(OptPop$Performance$Yields[length(OptPop$Performance$Yields)])
      
      OptimalConditions$Biomass[f]<-(sum(WeightAtAge %*% OptPop$FinalNumAtAge))
      
      OptimalConditions$Numbers[f]<-(sum(OptPop$FinalNumAtAge))
      
      MPAs$StatusQuo<- 0
      
      MPAs$EqNTZ<- c(0,rep(OptNTZSize$par,OptTime))
      
      MPAs$SNTZ<- c(0,seq(min(1,1.5*OptNTZSize$par),OptNTZSize$par,length.out=OptTime))
      
      MPAs$GNTZ<- c(0,seq(min(1,0.5*OptNTZSize$par),OptNTZSize$par,length.out=OptTime))
      
      #       MPAs$OptNTZ<- c(0, MPAFunction(OptMPAPath$par,1:OptTime, OptNTZSize$par,'LinearPlus',EvalTime))
      
      MPAs$OptNTZ<- c(0, MPAFunction(as.matrix(BestGuess),1:OptTime, OptNTZSize$par,'LinearPlus',EvalTime))
      
      MPAs$CatchShareEqNTZ<- MPAs$EqNTZ
      
      pdf(file=paste(FigureFolder,'FvFmsy is',round(FScenarios[f]/Fmsy$par,2),' MPAs.pdf'),family=Font,pointsize=12,width=6,height=4)
      par(mar=c(5.1,4.1,4.1,6.1),xpd=T)
      matplot((MPAs),type='l',col=rainbow(1.5*dim(MPAs)[2])[1:dim(MPAs)[2]],lty=1,xlab='Year',ylab='% NTZ',lwd=4,bty='n')
      legend('topright',inset=c(-.3,0),legend=MPANames,col=rainbow(1.5*dim(MPAs)[2])[1:dim(MPAs)[2]],lty=1,cex=.5)
      dev.off()
      
      for (m in 1:dim(MPAs)[2])
      {
        
        show(paste('M is ',m))
        
        TempPop<- StartPop
        Patches<- BasePatches
        
        cc<- cc+1
        for (y in 1:TimeToRun)
        {
          c<- c+1
          
          CurrentMPA<- MPAs[y,m]
          if (y>dim(MPAs)[1])
          {
            #             CurrentMPA<- MPAs[dim(MPAs)[2],m]
            CurrentMPA<-  OptNTZSize$par
            
          }
          
          AssignNTZ(CurrentMPA,ReservePosition)
          
          FTemp<- FScenarios[f] 
          
          FVec<- MoveFleet(FTemp,CurrentMPA,FleetSpill,0)
          
          if (MPANames[m]=='CatchShareEqNTZ')
          {
            #             FVec<- MoveFleet(FTemp,CurrentMPA,0,0)            
            FVec<- MoveFleet(Fmsy$par,CurrentMPA,0,0)            
            
          }
          
          PassPop<- GrowPopulation(TempPop,FVec,1,0,'TEST') #grow the new population
          
          TempPop<- PassPop$FinalNumAtAge
          
          #           NPBTimeSeries<- Discount(ResultStorage$Yield[1:y]-BaseConditions$Yield[f],Fleet$YieldDiscount,y)$NPV#Yield balance through 10
          
          TotalStorage[c,]<-data.frame(SpeciesList[s],MPANames[m],f,FScenarios[f],CurrentMPA,y,PassPop$Performance$MeanYield,PassPop$Performance$MeanBiomass,PassPop$Performance$MeanNumbers, BaseConditions$Yield[f],BaseConditions$Numbers[f], BaseConditions$Biomass[f], 
                                       OptimalConditions$Yield[f], OptimalConditions$Numbers[f], OptimalConditions$Biomass[f], OptNTZSize$par,stringsAsFactors=F)
          
        } #Close TimeToRun loop   
        
      } #Close loop over MPA proposals
      
    } #Close loop over fishing scenarios
    
    
  } #Close species list loop
  
  browser()
  TotalStorage$YieldBalance<- TotalStorage$Yield-TotalStorage$SQYield
  
  TotalStorage$ScenId<- with(TotalStorage,paste(Species,m,f,sep='-'))
  
  TotalStorage<- ddply(TotalStorage,c('ScenId'),mutate,
                       PresentYield=Yield*(1+Fleet$YieldDisc)^-(Year-1),PresentBalance=(Yield-SQYield)*(1+Fleet$YieldDisc)^-(Year-1),
                       NPY=cumsum(PresentYield),NPB=cumsum(PresentBalance),RequestedLoan = sum(PresentBalance[YieldBalance<0]))
  quartz()
ggplot(TotalStorage,aes(Year,NPB,color=m))+geom_line()+facet_wrap(~Species)
  
  write.csv(file=paste(ResultFolder,'Total Results.csv'),TotalStorage)
  
  save.image(file=paste(ResultFolder,'Completed Workspace.rdata'))
  
  AllSpeciesExperimentResults$MPAScenario<- as.factor(AllSpeciesExperimentResults$MPAScenario)
  
  levels(AllSpeciesExperimentResults$MPAScenario)<- MPANames
  
  AllSpeciesExperimentResults<- AllSpeciesExperimentResults[AllSpeciesExperimentResults$MPAScenario!='StatusQuo',]
  
  TotalStorage$Year<- as.numeric(TotalStorage$Year)
  
  TotalStorage$NPB<- as.numeric(TotalStorage$NPB)
  
  TotalStorage$m<- as.factor(TotalStorage$m)
  
  TotalStorage$f<- as.factor(TotalStorage$f)
  
  levels(TotalStorage$m)<- MPANames
  
  #   levels(TotalStorage$f)<- c('F50','F25')
  
  levels(TotalStorage$f)<- c('F25')
  
  AllSpeciesStorage$Year<- as.numeric(AllSpeciesStorage$Year)
  
  AllSpeciesStorage$NPB<- as.numeric(AllSpeciesStorage$NPB)
  
  AllSpeciesStorage$m<- as.factor(AllSpeciesStorage$m)
  
  AllSpeciesStorage$f<- as.factor(AllSpeciesStorage$f)
  
  levels(AllSpeciesStorage$m)<- MPANames
  
  levels(AllSpeciesStorage$f)<- c('F25')
  #   levels(AllSpeciesStorage$f)<- c('F50','F25')
  
  SubScenarios<- c('EqNTZ','GNTZ','SNTZ','OptNTZ','CatchShareEqNTZ')
  
  Where<-   AllSpeciesStorage$m %in% SubScenarios
  
  PlotStorage<- AllSpeciesStorage[Where,]
  
  PlotStorage$m<- as.character(PlotStorage$m)
  
  PlotStorage$m<- as.factor(PlotStorage$m)
  
  PlotStorage$FracNTZ<- as.numeric(PlotStorage$FracNTZ)
  
  PlotStorage$Yield<- as.numeric(PlotStorage$Yield)
  
  PlotStorage$SQYield<- as.numeric(PlotStorage$SQYield)
  
  PlotStorage$PositiveYields<- 0
  
  FNames<- levels(TotalStorage$f)
  
  #   PlotStorage$YieldBalance<- as.numeric(PlotStorage$YieldBalance)
  
  # levels(PlotStorage$m)<- c('Equilibrium','Grow','Optimize','Shrink')
  PlotStorage$OldM<- PlotStorage$m
  
  levels(PlotStorage$m)<- c('Catch Share','Long Term Optimal','Grow','Short Term Optimal','Shrink')
  
  
} #Run Analysis If
if (RunAnalysis==0)
{
  load(paste(BatchFolder,'Completed Workspace.rdata'))
}

# quartz()
# ggplot(data=PlotStorage,aes(x=Year,y=NPB,color=m))+geom_line()+facet_wrap(~Species)

Font<- "Helvetica" #Font type for figures


# levels(PlotStorage$m)<- c('Long Term Optimal','Grow','Shrink','Short Term Optimal')

# SubScenarios<- c('Long Term Optimal','Short Term Optimal','Grow','Shrink')

SubScenarios<- c('Long Term Optimal','Short Term Optimal','Grow','Shrink','Catch Share')

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
      
      LoanTime<- 20
      
      RequestedLoan<- Discount(pmin(0,YieldBalance[1:LoanTime]),Fleet$YieldDiscount,LoanTime)$NPV #Discounted Requested Loan Amount. 
      
      AvailableSurplus<- Discount(pmax(0,YieldBalance[1:LoanTime]),Fleet$YieldDiscount,LoanTime)$NPV #Discounted Requested Loan Amount. 
      
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

PlotStorage<- PlotStorage[,colnames(PlotStorage)!='Bmsy']

PlotStorage<- join(PlotStorage,SystemBmsyStorage,by='Species')

PlotStorage$PosYears<- PlotStorage$Year * PlotStorage$PositiveYields

PlotStorage$PosNPB<- PlotStorage$NPB * PlotStorage$PositiveYields

TimeToNPB<- ddply(PlotStorage[PlotStorage$f=='F25',],c('Species','m'),summarize,BenefitYear=which(NPB>0)[1],PercentYieldIncrease=round(100*(Yield[max(Year)]/Yield[1]-1),0),TotalCatch=sum(Yield),Fish=sum(Biomass))

#theme_get
FontColor<- 'black'

EQTimeToNPB<- subset(TimeToNPB,m=='Equilibrium')

SpeciesOrder<- order(EQTimeToNPB$PercentYieldIncrease)

KeynoteTheme<- theme(legend.position='top',plot.background=element_rect(color=NA),rect=element_rect(fill='transparent',color=NA),text=element_text(size=22,family=Font,color=FontColor),
                     axis.text=element_text(color=FontColor),axis.title.y=element_text(size=25,hjust=0.5,angle=0),axis.text.y=element_text(size=30),axis.text.x=element_text(angle=35, vjust=0.9,hjust=0.9,color=FontColor,size=22),
                     legend.text=element_text(size=14,color='black'),legend.title=element_text(size=16,color='black'),legend.background=element_rect(fill="gray90"))

PaperTheme<- theme(legend.position='top',text=element_text(size=22,family=Font,color=FontColor),
                   axis.text=element_text(color=FontColor),axis.title.y=element_text(size=25,hjust=0.5,angle=0),axis.text.y=element_text(size=30),axis.text.x=element_text(angle=35, vjust=0.9,hjust=0.9,color=FontColor,size=22),
                   legend.text=element_text(size=14,color='black'),legend.title=element_text(size=16,color='black'))

TimeToNPB$Species<- as.factor(TimeToNPB$Species)

TimeToNPB$Species <- reorder(TimeToNPB$Species, 1/TimeToNPB$BenefitYear)


pdf(file=paste(BatchFolder,'Time to Positive NPB Barchart.pdf',sep=''),width=12,height=10)

# TimePlot<- ggplot(data=subset(TimeToNPB,m=='Equilibrium'), aes(x = factor(Species), y = BenefitYear,fill=round(PercentYieldIncrease,2)))+KeynoteTheme + geom_bar(stat = "identity")
TimePlot<- ggplot(data=subset(TimeToNPB,m=='Long Term Optimal'), aes(x = factor(Species), y = BenefitYear)) +PaperTheme + geom_bar(fill='dodgerblue3',stat = "identity")

print(TimePlot+xlab('')
      +ylab('Years to Net Benefit')+geom_hline(aes(yintercept=sum(BenefitYear*TotalCatch)/sum(TotalCatch),size=3),color='red2'))
dev.off()



KeynoteTheme<- theme(legend.position='top',plot.background=element_rect(color=NA),rect=element_rect(fill='transparent',color=NA),text=element_text(size=22,family=Font,color=FontColor),
                     axis.text=element_text(color=FontColor),axis.title.y=element_text(size=25),axis.text.y=element_text(size=30),axis.title.x=element_text(size=25),axis.text.x=element_text(angle=35, vjust=0.9,hjust=0.9,color=FontColor,size=25),
                     legend.text=element_text(size=20,color='black'),legend.title=element_text(size=16,color='black'),legend.background=element_rect(fill="gray90"),legend.key.size = unit(1.5, "cm"))

KeynoteTheme2<- theme(legend.position='top',plot.background=element_rect(color=NA),rect=element_rect(fill='transparent',color=NA),text=element_text(size=22,family=Font,color=FontColor),
                      axis.text=element_text(color=FontColor),axis.title.y=element_text(size=25,angle=0),axis.text.y=element_text(size=30),axis.title.x=element_text(size=25),axis.text.x=element_text(color=FontColor,size=25),
                      legend.text=element_text(size=20,color='black'),legend.title=element_text(size=25,color='black'),legend.background=element_rect(fill="gray90"),legend.key.size = unit(3, "cm"))

PaperTheme2<- theme(legend.position='top',text=element_text(size=22,family=Font,color=FontColor),
                    axis.text=element_text(color=FontColor),axis.title.y=element_text(size=25,angle=0),axis.text.y=element_text(size=30),axis.title.x=element_text(size=25),axis.text.x=element_text(color=FontColor,size=25),
                    legend.text=element_text(size=20,color='black'),legend.title=element_text(size=25,color='black'),legend.key.size = unit(3, "cm"))



pdf(file=paste(BatchFolder,'Time to Positive NPB for EQ and Opt Barchart.pdf',sep=''),width=12,height=10)
TimePlot<- ggplot(data=subset(TimeToNPB,m=='Long Term Optimal' | m=='Short Term Optimal'), aes(x = Species, y = BenefitYear,fill=m)) + geom_bar(stat = "identity",position='dodge')

(TimePlot+PaperTheme+xlab('')
 +ylab('Years to Net Benefit')+scale_fill_manual(values=c('green4','deepskyblue3'),name=""))
dev.off()


PlotStorage$MaxInterestRate[PlotStorage$MaxInterestRate>25]<- 25

PlotStorage$LoanType[PlotStorage$MaxInterestRate<5]<- 'Philanthropy'

PlotStorage$LoanType[PlotStorage$MaxInterestRate>=5 & PlotStorage$MaxInterestRate<15]<- 'Social Investment'

PlotStorage$LoanType[PlotStorage$MaxInterestRate>=15]<- 'Commercial'

pdf(file=paste(BatchFolder,'Max Loan Payable Barchart.pdf',sep=''),width=12,height=10)

(ggplot(data=subset(PlotStorage,Year==max(Year) & f=='F25' & (m=='Long Term Optimal' | m=='Short Term Optimal')),aes(x=m,y=MaxInterestRate,fill=LoanType))
 +geom_bar(stat='identity')+facet_wrap(~Species)+PaperTheme+xlab('')+ylab('Max % Interest Rate')
 +theme(strip.text=element_text(color='black'),axis.text.x=element_text(size=20))+scale_fill_manual(values=c('goldenrod1','green4','deepskyblue3'),name=""))
dev.off()

pdf(file=paste(BatchFolder,'Required Price Increase Barchart.pdf',sep=''),width=12,height=10)

(ggplot(data=subset(PlotStorage,Year==max(Year) & f=='F25' & (m=='Long Term Optimal' | m=='Short Term Optimal')),aes(x=m,y=PriceIncreaseNeeded,fill=Species))
 +geom_bar(stat='identity')+facet_wrap(~Species)+PaperTheme+xlab('')+ylab('% Price Increase Needed')
 +theme(strip.text=element_text(color='black'))+scale_fill_discrete(name="")+theme(legend.position="none"))
dev.off()

PlotStorage$BadYears<- PlotStorage$YieldBalance<0

Breaks<- c(-50000,0,50000)

pdf(file=paste(BatchFolder,'Single Yield Trajectory.pdf',sep=''),width=16,height=10)
a=(ggplot(data=subset(PlotStorage, f=='F25' & Year<=20 & m=='Short Term Optimal' & Species=='Yellowtail Snapper'),
          aes(x=Year,y=YieldBalance,fill=NPB))+xlab('Years with Reserve')+ylab("Change in Catch")
   +geom_bar(stat='identity',color='black') +
     scale_fill_gradient2(low='red',mid='white',high='black',name='Total Benefits',breaks=Breaks,labels=c('-50,000','0','50,000'))+PaperTheme2
   +geom_hline(aes(yintercept=0),size=2))
a+theme(text=element_text(size=20),legend.position='right')
dev.off()


pdf(file=paste(BatchFolder,'Single Biomass Trajectory.pdf',sep=''),width=16,height=10)
a=(ggplot(data=subset(PlotStorage, f=='F25' & Year<=15 & m=='Long Term Optimal' & Species=='Yellowtail Snapper'),
          aes(x=Year,y=Biomass))+xlab('Years with Reserve')+ylab("Fish")
   +geom_bar(stat='identity',color='black',fill='steelblue2') +
     PaperTheme2)
a+theme(text=element_text(size=20))
dev.off()


pdf(file=paste(BatchFolder,'Single Yield Trajectory Line.pdf',sep=''),width=16,height=10)
a=(ggplot(data=subset(PlotStorage, f=='F25' & Year<=15 & m=='Long Term Optimal' & Species=='Yellowtail Snapper'),
          aes(x=Year,y=YieldBalance,color=NPB))+xlab('Years with Reserve')+ylab("Change in Catch")
   +geom_line(size=4) +PaperTheme2+scale_color_gradient(low='red',high='black',name='Total Benefit')
   +geom_hline(aes(yintercept=0),size=2))
a+theme(text=element_text(size=20),legend.position='right')
dev.off()

pdf(file=paste(BatchFolder,'Single NPB Trajectory Line.pdf',sep=''),width=16,height=10)
a=(ggplot(data=subset(PlotStorage, f=='F25' & Year<=15 & m=='Long Term Optimal' & Species=='Yellowtail Snapper'),
          aes(x=Year,y=NPB,color=NPB))+xlab('Years with Reserve')+ylab("Change in Catch")
   +geom_line(size=4) +PaperTheme2+scale_color_gradient(low='red',high='black',name='Total Benefit')
   +geom_hline(aes(yintercept=0),size=2))
a+theme(text=element_text(size=20),legend.position='right')
dev.off()



pdf(file=paste(BatchFolder,'NPB Trajectory.pdf',sep=''),width=12,height=10)
a=(ggplot(data=subset(PlotStorage, f=='F25' & Year<=65),aes(x=Year,y=NPB,color=m))
   +geom_line(size=2)+facet_wrap(~Species,scales='free')+geom_hline(aes(yintercept=0))+scale_color_discrete(name=""))
a+theme(text=element_text(size=20))
dev.off()


pdf(file=paste(BatchFolder,'Yield Trajectory.pdf',sep=''),width=12,height=10)
a=(ggplot(data=subset(PlotStorage, f=='F25' & Year<=25),aes(x=Year,y=Yield,color=m))
   +geom_line(size=2)+facet_wrap(~Species,scales='free')+geom_hline(aes(yintercept=Yield[1]))+scale_color_discrete(name=""))
a+theme(text=element_text(size=20))
dev.off()

pdf(file=paste(BatchFolder,'Reserve Trajectory.pdf',sep=''),width=12,height=10)
a=(ggplot(data=subset(PlotStorage, f=='F25' & Year<=75),aes(x=Year,y=100*FracNTZ,color=m))
   +geom_line(size=2)+facet_wrap(~Species,scales='free')+scale_color_discrete(name="")+ylab('% in Reserve'))
a+theme(text=element_text(size=20))+PaperTheme
dev.off()



MPASize<- ddply(subset(PlotStorage,f=='F25'),c('Species'),summarize,OptMPA= 100*unique(OptNTZ))

MPASize$X<- 50

MPASize$Y<- 50

ggplot(data=MPASize,aes(x=X,y=Y,size=OptMPA,color=Species,alpha=0.1))+scale_size_continuous(range=c(2*min(MPASize$OptMPA),2*max(MPASize$OptMPA)),guide=FALSE)+geom_point()

TimeToNPB$BenefitYear[is.na(TimeToNPB$BenefitYear)]<- TimeToRun

# pdf(file=paste(BatchFolder,'Time to Positive NPB.pdf',sep=''))
# barchart(~BenefitYear | Species,group=m,data=TimeToNPB,auto.key=T,xlab='Years to Positive NPB')
# dev.off()


Cols<-   scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65, h.start = 0, direction = 1)

plotyears=max(PlotStorage$Year)  #no. years to plot


# pdf(file=paste(BatchFolder,'Aggregate NPB Trajectory.pdf',sep=''))
# print(xyplot(NPB ~ Year | Species,group=m,data=PlotStorage,scales=list(y='free'),ylab='Net Profit Balance ($)',
#              subset=f=='F25' & Year<= plotyears,auto.key=T,type='l',lwd=2,panel=function(x,y,...)
#              {
#                panel.xyplot(x,y,...)
#                panel.abline(h=0,lwd=2,lty=2)
#              }))
# dev.off()

# pdf(file=paste(BatchFolder,'Aggregate PB Trajectory.pdf',sep=''))
# 
# print(xyplot(YieldBalance ~ Year | Species,group=m,data=PlotStorage,scales=list(y='free'),ylab='Profit Balance ($)',
#              subset=f=='F25' & Year<= plotyears,auto.key=T,type='l',lwd=2,panel=function(x,y,...)
#              {
#                panel.xyplot(x,y,...)
#                panel.abline(h=0,lwd=2,lty=2)
#              }))
# 
# dev.off()

# pdf(file=paste(BatchFolder,'Aggregate BvBmsy Trajectory.pdf',sep=''))
# 
# print(xyplot(as.numeric(Biomass)/(Bmsy)  ~ Year | Species,group=m,data=PlotStorage,ylab='B/Bmsy',
#              subset=f=='F25' & Year<= plotyears,auto.key=T,type='l',lwd=2,panel=function(x,y,...)
#              {
#                panel.xyplot(x,y,...)
#                panel.abline(h=1,lwd=2,lty=2)
#              }))
# 
# dev.off()

# pdf(file=paste(BatchFolder,'Tradeoff Analysis.pdf',sep=''))
# 
# print(xyplot((Biomass)/(Bmsy)  ~ NPB | Species,group=m,data=PlotStorage,ylab='B/Bmsy',
#              subset=f=='F25' & Year== max(Year,na.rm=T),cex=2,cex.axis=0.75,pch=19,auto.key=T,lwd=2,scales=list(x='free'),
#              panel=function(x,y,...)
#              {
#                panel.xyplot(x,y,...)
#                panel.abline(v=0,lty=2)
#                panel.abline(h=1,lty=2)
#                
#              }))
# 
# dev.off()

# pdf(file=paste(BatchFolder,'Reserve Trajectory Analysis.pdf',sep=''))
# 
# print(xyplot(100*FracNTZ  ~ Year | Species,group=m,data=PlotStorage,ylab='% Reserve',
#              subset=f=='F25' & Year<= plotyears,pch=19,auto.key=T,type='l',lwd=2,
#              panel=function(x,y,...)
#              {
#                panel.xyplot(x,y,...)
#                #                panel.abline(v=0,lty=2)
#                #                panel.abline(h=1,lty=2)
#                
#              }))
# 
# dev.off()

PlotStorage$MaxInterestRate[PlotStorage$MaxInterestRate>100]<- 100

PlotStorage$MaxInterestRate[is.na(PlotStorage$MaxInterestRate)]<- 100

scales=list(labels=c("",paste(seq(0,80,by=20),'%',sep=''),'>=100%'))

# pdf(file=paste(BatchFolder,'Aggregate Financial Analysis.pdf',sep=''))
# 
# 
# IntPlot<- (barchart(~MaxInterestRate | Species,scales=list(cex=0.7,labels=c("",paste(seq(0,80,by=20),'%',sep=''),'>100%'))
#                     ,group=m,main='A',data=PlotStorage,auto.key=T,subset=Year==1 & f=='F25',par.strip.text=list(cex=.75),xlab=list(label='Maximum % Interest Rate',fontsize=10),
#                     
# ))

# PlotStorage$PriceIncreaseNeeded[is.na(PlotStorage$PriceIncreaseNeeded)]<- 0
# 
# a<- ddply(PlotStorage,c('f','m','Species'),summarize,BenefitTime=which(NPB>0)[1])
# 
# PricePlot<- (barchart(~PriceIncreaseNeeded | Species,group=m,main='B',data=PlotStorage,scales=list(cex=0.7,labels=c("",paste((seq(0,round(max(PlotStorage$PriceIncreaseNeeded,na.rm=T),-1),length.out=6)),'%',sep='')))
#                       ,subset=Year==1 & f=='F25',par.strip.text=list(cex=.75),auto.key=F,xlab=list(label='% Price Increase Needed',fontsize=10)))
# 
# grid.arrange(IntPlot,PricePlot, nrow=2)
# 
# dev.off()
# 
# pdf(file=paste(BatchFolder,'Aggregate EQ Reserve.pdf',sep=''))
# 
# # PlotStorage$OptNTZ<- 100*as.numeric(PlotStorage$OptNTZ)
# 
# print(barchart(~(100*OptNTZ), group= Species,data=AllSpeciesExperimentResults,scales=list(cex=0.7,labels=c("",paste(seq(0,100,by=20),'%',sep='')))
#                ,subset= FishingScenario==2 & MPAScenario=='EqNTZ',auto.key=T,xlab='Optimal Equilibrium Reserve %'))
# 
# dev.off()

SubScenarios<- c('Equilibrium','Grow','Shrink','Optimize')

for (s in 1:length(SpeciesList))
{
  
  
  MaxNPB<- max(PlotStorage$NPB[PlotStorage$Year<=plotyears])
  
  FigureFolder<- paste(BatchFolder,SpeciesList[s],'/Figures/',sep='')
  ResultFolder<- paste(BatchFolder,SpeciesList[s],'/Results/',sep='')
  
  #   pdf(file=paste(FigureFolder,'NPB Trajectory.pdf',sep=''),col='rgb')
  #   #subset data for each scenario (m)
  #   m1=subset(PlotStorage,m=="Equilibrium" & Species== SpeciesList[s])
  #   m2=subset(PlotStorage,m=="Grow" & Species == SpeciesList[s])
  #   m3=subset(PlotStorage,m=="Shrink"& Species == SpeciesList[s])
  #   m4=subset(PlotStorage,m=="Optimize"& Species == SpeciesList[s])
  #   
  #   # m4G=subset(PlotStorage,m=="OptPath"& Species == SpeciesList[2])
  #   
  # #   y1_FirstYields<-(as.numeric(PlotStorage$YieldBalance[PlotStorage$f=="F50" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s]]))
  #   
  #   y2_FirstYields<-(as.numeric(PlotStorage$YieldBalance[PlotStorage$f=="F25" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s]]))
  #   
  #   
  # #   MaxYields<- max(y1_FirstYields, y2_FirstYields)
  # 
  # MaxYields<- max(y2_FirstYields)
  # 
  # 
  # #   y1_FirstYields<- y1_FirstYields/MaxYields
  #   y2_FirstYields<- y2_FirstYields/MaxYields
  #   
  #   #get points where yield equals status quo
  # #   y1=PlotStorage[PlotStorage$f=="F50" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s],c(6,18) ]
  #   y2=PlotStorage[PlotStorage$f=="F25" & PlotStorage$PositiveYields==1 & PlotStorage$Species== SpeciesList[s],c(6,18)]
  #   
  #   par(mfrow=c(2,1),mar=c(0,0,0,0),oma=c(4,6,2,6))
  #   
  #   plot(m1$Year[as.vector(m1$f)=="F50"],m1$NPB[as.vector(m1$f)=="F50"],
  #        type="l",ylab="",xlab="",ylim=c(min(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears]),1.25*max(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears])),
  #        xlim=c(1,plotyears),col="dodgerblue",lwd=4,las=1,xaxt="n" )
  #   mtext("Moderate Overfishing",line=-1)
  #   text(1,max(PlotStorage$NPB)*0.99,"A.")
  #   lines( m2$Year[as.vector(m2$f)=="F50"],m2$NPB[as.vector(m2$f)=="F50"],col=2,lwd=2,lty=3)
  #   lines( m3$Year[as.vector(m3$f)=="F50"],m3$NPB[as.vector(m3$f)=="F50"],col="medium orchid",lwd=2,lty=4)
  #   lines( m4$Year[as.vector(m4$f)=="F50"],m4$NPB[as.vector(m4$f)=="F50"],col="darkgreen",lwd=2,lty=5)
  #   abline(h=0,lty=2)
  #   points(y1[,1],y1[,2],col=c("dodgerblue","medium orchid",2,"darkgreen"),pch=16,cex=2)
  #   
  #   plot(m1$Year[as.vector(m1$f)=="F25"],m1$NPB[as.vector(m1$f)=="F25"],
  #        type="l",ylab="",xlab="Year",ylim=c(min(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears ]),1.25*max(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears])),
  #        xlim=c(1,plotyears),col="dodgerblue",lwd=4 ,las=1)
  #   mtext("Heavy Overfishing",line=-1)
  #   text(1,max(PlotStorage$NPB)*0.99,"B.")
  #   lines( m2$Year[as.vector(m2$f)=="F25"],m2$NPB[as.vector(m2$f)=="F25"],col=2,lwd=2,lty=3)
  #   lines( m3$Year[as.vector(m3$f)=="F25"],m3$NPB[as.vector(m3$f)=="F25"],col="medium orchid",lwd=2,lty=4)
  #   lines( m4$Year[as.vector(m4$f)=="F25"],m4$NPB[as.vector(m4$f)=="F25"],col="darkgreen",lwd=2,lty=5)
  #   abline(h=0,lty=2)         
  #   points(y2[,1],y2[,2],col=c("dodgerblue","medium orchid",2,"darkgreen"),pch=16, cex=2)
  #   
  #   par(new=TRUE, mfrow=c(1,1), oma=c(0,0,0,0))
  #   
  #   plot(1,
  #        xlim=c(0,1),
  #        ylim=c(0,1),
  #        xaxs="i",
  #        yaxs="i",
  #        type="n",
  #        axes=FALSE)
  #   
  #   MaxCircle<- round((MaxYields),0) 
  #   
  #   Difference<- 10-MaxCircle%%10
  #   
  #   MaxCircle <- MaxCircle + Difference
  #   
  #   CircleSize<- ((seq(from=MaxCircle/10,to=MaxCircle,length.out=5)))
  #   
  #   legend(0.85,0.9, legend= SubScenarios, lty=c(1,3,4,5),
  #          col=c("dodgerblue",2,"medium orchid","darkgreen"),lwd=2,bty="n",cex=0.8)
  #   
  #   # legend(0.85,0.75, legend= CircleSize,pch=16,pt.cex=(CircleSize),col=c("black"),title='Yield Surplus',bty='n',y.intersp=1.5)
  #   
  #   
  #   text(.5,.05,"Year",font=2)
  #   
  #   text(.03,.5,"Net Present Balance ($)",srt=90,font=2)
  #   
  #   dev.off()
  
  pdf(file=paste(FigureFolder,'FracNTZ Trajectory.pdf',sep=''))
  print(xyplot(FracNTZ ~ Year | f,group=m,data=PlotStorage,type='l',subset=(Year<=11 & PlotStorage$Species==SpeciesList[s]),auto.key=list(space='right',points=F,lines=T),layout=c(1,2),index.cond=list(c(2,1)),lwd=4,drop.unused.levels=T,panel=function(x,y,...)
  {panel.xyplot(x,y,...)
   panel.abline(h=0,lty=2)
  },))
  dev.off()
  
  #   pdf(file=paste(FigureFolder,'Heavy Fishing NPB Trajectory.pdf',sep=''),col='rgb')
  #   #subset data for each scenario (m)
  #   m1=subset(PlotStorage,m=="Equilibrium" & Species== SpeciesList[s])
  #   m2=subset(PlotStorage,m=="Grow" & Species == SpeciesList[s])
  #   m3=subset(PlotStorage,m=="Shrink"& Species == SpeciesList[s])
  #   m4=subset(PlotStorage,m=="Optimize"& Species == SpeciesList[s])
  #   
  #   # m4G=subset(PlotStorage,m=="OptPath"& Species == SpeciesList[2])
  #   
  #   y1_FirstYields<-(as.numeric(PlotStorage$YieldBalance[PlotStorage$f=="F50" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s]]))
  #   
  #   y2_FirstYields<-(as.numeric(PlotStorage$YieldBalance[PlotStorage$f=="F25" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s]]))
  #   
  #   
  #   MaxYields<- max(y1_FirstYields, y2_FirstYields)
  #   
  #   y1_FirstYields<- y1_FirstYields/MaxYields
  #   y2_FirstYields<- y2_FirstYields/MaxYields
  #   
  #   #get points where yield equals status quo
  #   y1=PlotStorage[PlotStorage$f=="F50" & PlotStorage$PositiveYields==1 &  PlotStorage$Species== SpeciesList[s],c(6,18) ]
  #   y2=PlotStorage[PlotStorage$f=="F25" & PlotStorage$PositiveYields==1 & PlotStorage$Species== SpeciesList[s],c(6,18)]
  #   
  #   par(mar=c(0,0,0,0),oma=c(4,6,2,6))
  #   
  #   plot(m1$Year[as.vector(m1$f)=="F25"],m1$NPB[as.vector(m1$f)=="F25"],
  #        type="l",ylab="",xlab="Year",ylim=c(min(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears ]),1.25*max(PlotStorage$NPB[PlotStorage$Species== SpeciesList[s]& PlotStorage$Year<= plotyears])),
  #        xlim=c(1,plotyears),col="dodgerblue",lwd=4 ,las=1)
  #   #   mtext("Heavy Overfishing",line=-1)
  #   text(1,max(PlotStorage$NPB)*0.99,"B.")
  #   lines( m2$Year[as.vector(m2$f)=="F25"],m2$NPB[as.vector(m2$f)=="F25"],col=2,lwd=2,lty=3)
  #   lines( m3$Year[as.vector(m3$f)=="F25"],m3$NPB[as.vector(m3$f)=="F25"],col="medium orchid",lwd=2,lty=4)
  #   lines( m4$Year[as.vector(m4$f)=="F25"],m4$NPB[as.vector(m4$f)=="F25"],col="darkgreen",lwd=2,lty=5)
  #   abline(h=0,lty=2)         
  #   points(y2[,1],y2[,2],col=c("dodgerblue","medium orchid",2,"darkgreen"),pch=16, cex=2)
  #   
  #   par(new=TRUE, mfrow=c(1,1), oma=c(0,0,0,0))
  #   
  #   plot(1,
  #        xlim=c(0,1),
  #        ylim=c(0,1),
  #        xaxs="i",
  #        yaxs="i",
  #        type="n",
  #        axes=FALSE)
  #   
  #   MaxCircle<- round((MaxYields),0) 
  #   
  #   Difference<- 10-MaxCircle%%10
  #   
  #   MaxCircle <- MaxCircle + Difference
  #   
  #   CircleSize<- ((seq(from=MaxCircle/10,to=MaxCircle,length.out=5)))
  #   
  #   legend(0.85,0.9, legend= SubScenarios, lty=c(1,3,4,5),
  #          col=c("dodgerblue",2,"medium orchid","darkgreen"),lwd=2,bty="n",cex=0.7)
  #   
  #   # legend(0.85,0.75, legend= CircleSize,pch=16,pt.cex=(CircleSize),col=c("black"),title='Yield Surplus',bty='n',y.intersp=1.5)
  #   
  #   
  #   text(.5,.05,"Year",font=2)
  #   
  #   text(.03,.5,"Net Present Balance ($)",srt=90,font=2)
  #   
  #   dev.off()
  #   
  #   pdf(file=paste(FigureFolder,'Heavy Overfishing FracNTZ Trajectory.pdf',sep=''))
  #   print(xyplot(100*FracNTZ ~ Year,group=m,data=PlotStorage,ylab=' % in Reserve',type='l',subset=(Year<=11 & PlotStorage$Species==SpeciesList[s] & f=='F25'),auto.key=list(space='right',points=F,lines=T),lwd=4,drop.unused.levels=T,panel=function(x,y,...)
  #   {panel.xyplot(x,y,...)
  #    panel.abline(h=0,lty=2)
  #   },))
  #   dev.off()
  
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

