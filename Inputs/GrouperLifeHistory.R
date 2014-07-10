####### Grouper Life History ########

#### MPASHRINK  Life history Parameters

###Lutjanus malabaricus from Blaber et al. 2005
lh<- NULL
lh$Nmsy<- -999
lh$MaxAge <- 25 
lh$m <- rep(0.2,lh$MaxAge)  # blaber et al estimate  Z= .11    so i cut this in half for m???
lh$CarryingCapacity<- -999
 
###VonBert Growth Params
lh$Linf <- 55.3  
lh$k <- 0.1 
lh$t0 <- -.6 
lh$VBSD<- 0.1 #Base standard deviation of the length curve
lh$VBErrorMean<- 0
lh$VBErrorSlope <- 0.1 #Sets the responsiveness of the SD parameter to age 0 nullifies effect, 1 makes increase in SD decrease with age

###Weight at Age
lh$wa <- 0.017  
lh$wb <- 3.012  
   
###Fecundity at Age
# lh$fa <- 0.229  
lh$fa <- 0.829  
# lh$fb <- 3.620  
lh$fb <- 1.62
lh$NoFecundRelate<- 1

###Maturity at Age
 # lh$ma50 <- -0.487  
# lh$ma50 <- -0.227  
# lh$ma95 <- 13.523  
 lh$ma50 <- 35 #The length at which 50% of individuals are sexually mature
lh$maBeta<- 3 # steepness parameters, the lower the value, the more knife edged the relationship is

###Movement
lh$Range<- 0.2
lh$MoveType<- 'Simple'


###Recruitment
lh$DDForm<- 'BH'
lh$B0<- 10000 #Virgin biomass
lh$R0<- lh$B0*(1-exp(-lh$m[1])) #Recruitment must equal deaths at EQ
lh$SexRatio<- 0.5
lh$LarvalChoice<- 0 #Whether larvae can choose to target good habitat. 1 means that eggs move towards better habitat, 0 means that recruitment is affected by habitat
lh$BH.Steepness<- 0.45 #Beverton-holt steepness parameter, estimated from several sebastes in RAM
lh$RecDevMean<- 0 #Mean log recruitment deviate (basically leave at 0)
lh$RecDevSTD<- 0 #Standard deviation of recruitment deviations

 ####### Set Fishing Fleet/Management Parameters ########
 Fleet<- NULL
 Fleet$YieldDiscount<- 0.05
 Fleet$BiomassDiscount<- 0.05
 Fleet$SizeLimit<- 1
 Fleet$s50<- .1
 Fleet$s95<- 1
