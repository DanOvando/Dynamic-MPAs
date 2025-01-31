
####### SET OPERATING PARAMETERS ########

InputFolder<- 'Inputs/' #Folder where the population parameter files are
SeedFolder<- 'Working/'#Folder where results will be stored (working version)


FontSize<- 14 #Font size for figures
Font<- "Helvetica" #Font type for figures

####### SET OPERATING PARAMETERS ########
PopTolerance<- .1 #the cutoff for identifying point where population stops changing
InitialPopulation<- 1000 #Seed for initial population 
CollapseThreshold<- 0.1
LookAtLengths<- 0
ReservePosition<- 'Center'
OptTime<- 65 #Time Horizon to optimize over
Alpha<- 0.5

####### Load in Population Parameters ########

# source(paste(InputFolder,Species,'LifeHistory.R',sep='')) #load in population parameters
source(paste(InputFolder,'GenericLifeHistory.R',sep='')) #load in population parameters

# source('GASP.R') #source GASP functions
lh$Bmsy<- -999
lh$SSB0<- NA
# LengthAtAge<- Length(1:lh$MaxAge) #Calculate length at age vector
# WeightAtAge<- Weight(LengthAtAge,lh$WeightForm) #Calculate weight at age vector
# FecundityAtAge<- Fecundity(LengthAtAge,'Length') #Calculate Fecundity at age vector
# if (lh$NoFecundRelate==1)
# {
#   FecundityAtAge<- WeightAtAge
# }
# MaturityAtAge<- Maturity(LengthAtAge,'Length',lh$LengthMa50, lh$LengthMa95) # Calculate maturity at age vector
# MaturityAtAge<- Maturity(1:lh$MaxAge, MaturityMode) # Calculate maturity at age vector


# ####### Set Fishing Fleet/Management Parameters ########
Fleet<- NULL
Fleet$YieldDiscount<- 0.1
Fleet$BiomassDiscount<- 0
Fleet$SizeLimit<- .001
Fleet$s50<- .001
Fleet$s95<- .002
FleetSpill<- 1


####### LOAD IN Habitat Structure and MPAs  ########
NumPatches<- 100 #number of patches
Patches<- NULL
Patches$PatchSizes<- rep(1,NumPatches) #Automatically sets up the right number of patches
Patches$SizeLocations<- 0 #Place the sizes where you want them. Entering 0 here will turn this off (leaving each patch size at default)
# Patches$c[Patches$SizeLocations]<- Patches$PatchSizes
Patches$MPALocations<- rep(0,NumPatches) #Locations of MPAs in the patches. 
Patches$HabQuality<- rep(1,NumPatches) #Habitat quality, change to vector of number of your choice if needed

lh$R0<- 5000

lh$R0<- lh$R0*Patches$HabQuality

if (max(Patches$SizeLocations)>NumPatches | max(Patches$MPALocations)>NumPatches)
{
  warning('PatchSizes, PatchLocations, or MPALocations are not the correct length; fix in controlfile')
}

# save.image(file=paste(SeedFolder,'ModelSettings.Rdata',sep=''))