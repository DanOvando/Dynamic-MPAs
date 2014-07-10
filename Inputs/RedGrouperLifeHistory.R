####### Red Grouper Life History ########


lh<- NULL
lh$Nmsy<- -999

lh$m <- 0.14 

lh$MaxAge <- 29

lh$m <- rep(lh$m,lh$MaxAge)  # blaber et al estimate  Z= .11    so i cut this in half for m???

lh$CarryingCapacity<- -999
 
###VonBert Growth Params
lh$Linf <- 854
lh$k <- 0.16
lh$t0 <- -.19 
lh$VBSD<- 0.1 #Base standard deviation of the length curve
lh$VBErrorMean<- 0
lh$VBErrorSlope <- 0.1 #Sets the responsiveness of the SD parameter to age 0 nullifies effect, 1 makes increase in SD decrease with age

###Weight at Age
lh$wa <- 6e-9
lh$wb <- 3.14
   
   
###Fecundity at Age
 lh$fa <- 0.229  
#lh$fa <- 0.829  
 lh$fb <- 3.620  
#lh$fb <- 1.62
lh$NoFecundRelate<- 0

###Maturity at Age

lh$LengthMa50<- 341
lh$LengthMa95<- 342
MaturityMode<- 'Age'
lh$AgeMa50<- 3
lh$AgeMa95<- 3.1

###Movement
lh$Range<- 0.3
lh$MoveType<- 'Simple'



###Recruitment
lh$DDForm<- 'BH'
lh$B0<- 1000 #Virgin biomass
lh$R0<- lh$B0*(1-exp(-lh$m[1])) #Recruitment must equal deaths at EQ
lh$SexRatio<- 0.5
lh$LarvalChoice<- 0 #Whether larvae can choose to target good habitat. 1 means that eggs move towards better habitat, 0 means that recruitment is affected by habitat
 lh$BH.Steepness<- 0.84 #Beverton-holt steepness parameter average of rockfish in RAM
#lh$BH.Steepness<- 0.7 #Beverton-holt steepness parameter average of rockfish in RAM

lh$RecDevMean<- 0 #Mean log recruitment deviate (basically leave at 0)
lh$RecDevSTD<- 0 #Standard deviation of recruitment deviations