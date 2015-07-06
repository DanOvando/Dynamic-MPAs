# Prepare Populations -----------------------------------------------------
# Function to prepare life history data for each species in operating model
PreparePopulations<- function(s,LifeHistories,LifeVars,BaseLife,BasePatches,BatchFolder,TargetPop)
{
  
  Species<- LifeHistories$CommName[gsub(' ','',LifeHistories$CommName,fixed = T) == s ]
  
  lh <- BaseLife
  
  StoreRun <- paste(BatchFolder,Species,sep = '') #1 if you want to create a seeded folder to store results, 0 if you want it in the generic working folder
  
  if ( StoreRun == 1 )
  {
    SeedFolder<- paste(Species,'Created-',Sys.time(),'/',sep="") #Folder to store outputs with a timestamp
  }
  if (is.character(StoreRun))
  {
    SeedFolder<- paste(StoreRun,'/',sep='')
  }
  
 if (dir.exists(SeedFolder)==F){ dir.create(SeedFolder)}
  
  FigureFolder<- paste(SeedFolder,'Figures/',sep='')
  
  ResultFolder<- paste(SeedFolder,'Results/',sep='')
  
  if (dir.exists(FigureFolder)==F){ dir.create(FigureFolder)}
  
  if (dir.exists(ResultFolder)==F){ dir.create(ResultFolder)}
  
  save.image(file=paste(SeedFolder,'ModelSettings.Rdata',sep=''))
  
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
  
  lh$MovementArray<- movArray(NumPatches,((lh$Range)*NumPatches)/5,'Wrap',FigureFolder)
  
  lh$Bmsy<- -999
  
  lh$LengthAtAge<- Length(1:lh$MaxAge,lh) #Calculate length at age vector
  
  lh$WeightAtAge<- Weight(lh$LengthAtAge,lh$WeightForm,lh) #Calculate weight at age vector
  
  lh$FecundityAtAge<- Fecundity(lh$WeightAtAge,'Weight',lh) #Calculate Fecundity at age vector
  
  if (lh$NoFecundRelate==1)
  {
    lh$FecundityAtAge<- lh$WeightAtAge
  }
  
  lh$MaturityMode<- as.character(LifeHistory$MaturityMode)
  
  lh$MaturityAtAge<- Maturity(1:lh$MaxAge, lh$MaturityMode,lh) # Calculate maturity at age vector
  
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
  

# Create Unfished Population ----------------------------------------------

  EQPopulation<- GrowPopulation(1000,rep(0,NumPatches),'EQ',0,'EQ Run',Species=Species,lh=lh,Patches=Patches,FigureFolder=FigureFolder) #Run the population out to unfished equilibrium
  
  lh<- EQPopulation$lh
  
  lh$CarryingCapacityWeight<- (colSums(EQPopulation$FinalNumAtAge*lh$WeightAtAge)) #Calculate carrying capacity in weight
  
  UnfishedPopulation<- EQPopulation$FinalNumAtAge #Unfished numbers at age

  lh$UnfishedPopulation<- EQPopulation
  

# Estimate Reference Points -----------------------------------------------


  Fmsy<- optimize(log(lh$m[1]),f=FindReferencePoint,Target='FMSY',TargetValue=NA,lower=-10,upper=4,Species=Species,lh=lh,UnfishedPopulation=UnfishedPopulation,Patches=Patches) #Find FMSY  
  
  Fmsy$par<- exp(Fmsy$minimum) 
  
  BmsyPopulation<- GrowPopulation(UnfishedPopulation,rep(Fmsy$par,NumPatches),'EQ',0,'Bmsy Run',Species=Species,lh=lh,Patches=Patches,FigureFolder=FigureFolder) #Bmsy Population
  
  lh$Fmsy<- Fmsy
  
  lh$Nmsy<- (colSums(BmsyPopulation$FinalNumAtAge)) #Nmsy
  
  lh$Bmsy<- colSums(BmsyPopulation$FinalNumAtAge*lh$WeightAtAge) #Bmsy
  
  FatTarget<- optimize(log(2*Fmsy$par),f=FindReferencePoint,Target='BvBmsy',TargetValue=TargetPop,lower=-10,upper=4,Species=Species,lh=lh,UnfishedPopulation=UnfishedPopulation,Patches=Patches) #Find F that results in target B/Bmsy
  
  FatTarget$par<- exp(FatTarget$minimum)
  
  lh$FatTarget<- FatTarget
  
  short_species<- gsub(' ','',Species,fixed = TRUE)
  
  eval(parse(text=paste(short_species,'<- lh',sep='')))

  return( list(eval(parse(text = paste(short_species, sep = '') ) ) ))

#   return( eval(parse(text = paste('list(', short_species, ' = ', short_species,')', sep = '') ) ) )
  
#   return( eval(parse(text = paste('unlist(', short_species, ' = ', short_species,')', sep = '') ) ) )
  
#   return( list( eval( parse( text = paste( short_species, ' = ', short_species, sep = '') ) ),wtf=2))
}

