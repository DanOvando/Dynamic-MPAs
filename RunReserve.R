RunReserve<- function(Run,BasePatches)
{
      

  for (d in 1:dim(Run)[2])
  {
    eval(parse(text=paste(colnames(Run)[d],'<- I(Run[d])',sep='')))
  }
  
  
       StoreRun<- paste(BatchFolder,Species,sep='') #1 if you want to create a seeded folder to store results, 0 if you want it in the generic working folder
         
         source('BetterNutzControlfile.R')
         
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
         
         LengthAtAge<- Length(1:lh$MaxAge) #Calculate length at age vector
         
         WeightAtAge<- Weight(LengthAtAge,lh$WeightForm) #Calculate weight at age vector
         
         FecundityAtAge<- Fecundity(WeightAtAge,'Weight') #Calculate Fecundity at age vector
         
         MaturityMode<- as.character(LifeHistory$MaturityMode)
         
         MaturityAtAge<- Maturity(1:lh$MaxAge, MaturityMode) # Calculate maturity at age vector
         
         if (Species=='Shark')
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
           
           OptNTZSize<-   optim(mguess,f=FindOptimalMPASize,lower=0,upper=0.999,FTemp=FScenarios[f],StartPop=StartPop,FleetSpill=FleetSpill) #You need a better optimization here, gets really stuck with any kind of stochasticity
           
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
               
               if (MPANames[m]=='CatchShareEqNTZ'& Year>1)
               {
                 #             FVec<- MoveFleet(FTemp,CurrentMPA,0,0)            
                 FVec<- MoveFleet(Fmsy$par,CurrentMPA,0,0)            
                 
               }
               
               PassPop<- GrowPopulation(TempPop,FVec,1,0,'TEST') #grow the new population
               
               TempPop<- PassPop$FinalNumAtAge
               
               #           NPBTimeSeries<- Discount(ResultStorage$Yield[1:y]-BaseConditions$Yield[f],Fleet$YieldDiscount,y)$NPV#Yield balance through 10
               
               TotalStorage[c,]<-data.frame(SpeciesList[s],MPANames[m],f,FScenarios[f],CurrentMPA,y,PassPop$Performance$MeanYield,PassPop$Performance$MeanBiomass,PassPop$Performance$MeanNumbers, BaseConditions$Yield[f],BaseConditions$Numbers[f], BaseConditions$Biomass[f], 
                                            OptimalConditions$Yield[f], OptimalConditions$Numbers[f], OptimalConditions$Biomass[f], OptNTZSize$par,stringsAsFactors=F)
               
               TotalStorage$YieldBalance<- TotalStorage$Yield-TotalStorage$SQYield
               
               TotalStorage$ScenId<- with(TotalStorage,paste(Species,m,f,sep='-'))
               
               TotalStorage<- ddply(TotalStorage,c('ScenId'),mutate,
                                    PresentYield=Yield*(1+Fleet$YieldDisc)^-(Year-1),PresentBalance=(Yield-SQYield)*(1+Fleet$YieldDisc)^-(Year-1),
                                    NPY=cumsum(PresentYield),NPB=cumsum(PresentBalance),RequestedLoan = sum(PresentBalance[YieldBalance<0])) %>% subset(m!='StatusQuo')
               #   quartz()
               
             } #Close TimeToRun loop   
             
           } #Close loop over MPA proposals
           
         } #Close loop over fishing scenarios
         
       } #Close species list loop

  
} #Close function