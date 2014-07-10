#Two-patch certification model
#patch A is open access
#patch B is NTZ

#######################
#model input parameters
#######################


B<-.1776			#percent of fish habitat in NTZ
A<-(1-B)			#percent of fish habitat in open access (out of fishery/range??)


Fmort<-.05                      #instantanious fishing mortality
#u<-1-exp(-Fmort)               #harvest rate in open access (iff equal pressure applied by fleet)
u<-0.0558                       #actual harvest rate
x<-30		     	   	         #maximum age considered
am<-4		    		           #age at maturity
af<-7                      #age enter fishery
v<-c(rep(0,af),rep(1,(x+1-af))) #vector of vulnerability at age 
im<-.25                     #instantanious natural mortality
m<-1-exp(-im)
surv<-c(1,rep((1-m),x))		     	 #vector of survival at age
f<-c(rep(0,am),rep(1,x+1-am))	 #vector of fecundity at age

Linf<-83.86                 #B-H parameter, asymptotic length
kpar<-0.068                 #B-H parameter growth rate
t0par<-0                    #B-H parameter, size at first settlement
age<-seq(0,x,1)             #vector of ages considered
len<-Linf*(1-exp(-kpar*(age-t0par)))   #vector of lengths

wtpar1<-.000026935         #weight at length parameter
wtpar2<-2.857              #weight at length parameter
wt<-wtpar1*len^wtpar2       #vector of mass at age

disp_ab<-.2275                #percentage larvae dispersing from A to B
disp_ba<-.7654               #percentage larvae dispersing from B to A
alpha1<-NA                  #BH parameter pre-dispersal
alpha2<-4.565334                  #BH parameter pre-dispersal 
beta1<-NA                   #BH parameter post-dispersal
beta2<-.00000298                  #BH parameter post-dispersal
m_ab<-0                   #pre-harvest movement from A to B
m_ba<-0                   #pre-harvest movement from B to A
n_ab<-0                   #post-harvest movement from A to B
n_ba<-0                   #post-harvest movement from B to A



#######################
#simulation parameters
#######################

years<-100			           #number of simulated years
Ninit<-100
initAbun<-rep(Ninit,x+1)		 #initial number of 


NstartA<-matrix(NA,nrow=years+1,ncol=x+1)
NstartA[1,]<-initAbun
NfishedA<-matrix(NA,nrow=years,ncol=x+1)
NsurvivedA<-matrix(NA,nrow=years,ncol=x+1)
NmovedA<-matrix(NA,nrow=years,ncol=x+1)


Eggs_A<-NA
Rec_pre_A<-NA
Rec_post_A<-NA
Settle_A<-NA

NstartB<-matrix(NA,nrow=years+1,ncol=x+1)
NstartB[1,]<-initAbun
NfishedB<-matrix(NA,nrow=years,ncol=x+1)
NsurvivedB<-matrix(NA,nrow=years,ncol=x+1)
NmovedB<-matrix(NA,nrow=years,ncol=x+1)

Eggs_B<-NA
Rec_pre_B<-NA
Rec_post_B<-NA
Settle_B<-NA

tempAgeA<-NA
tempAgeB<-NA
CatchAN<-NA               #catch in numbers
CatchAB<-NA               #catch in biomass

#######################
#SIMULATION
#######################

#sim<-function(param1,param2,u){

#loop over each year

for(y in 1:years){


	 
   
   #Produce Eggs
	 Eggs_A[y]=sum(f*NstartA[y,]*wt)    #f=0 until maturity, recruits and young don't produce eggs
	 Eggs_B[y]=sum(f*NstartB[y,]*wt)
	
	
	 #pre-dispersal density dependence IGNORE FOR FIRST ROUND
	 Rec_pre_A[y]=Eggs_A[y]
	 Rec_pre_B[y]=Eggs_B[y]
	
	 #Settlement
	 Settle_A[y]=Rec_pre_A[y]*(1-disp_ab)+Rec_pre_B[y]*disp_ba
   Settle_B[y]=Rec_pre_B[y]*(1-disp_ba)+Rec_pre_A[y]*disp_ab
	
	 #post-dispersal density dependence
	 Rec_post_A[y]=(alpha2*Settle_A[y])/(1+beta2*(1/A)*Settle_A[y])
	 Rec_post_B[y]=(alpha2*Settle_B[y])/(1+beta2*(1/B)*Settle_B[y])
	
   #Fish the fish
   NfishedA[y,]<-NstartA[y,]-NstartA[y,]*u*v
   CatchAN[y]<-sum(NstartA[y,]*u*v)
   CatchAB[y]<-sum(NstartA[y,]*u*v*wt)
   NfishedB[y,]<-NstartB[y,]
   
   
   #Survival
   NsurvivedA[y,]<-NfishedA[y,]*surv
   NsurvivedB[y,]<-NfishedB[y,]*surv
   
   #Move fish
   NmovedA[y,1]<-NsurvivedA[y,1]
	 NmovedA[y,2:(x+1)]=NsurvivedA[y,2:(x+1)]*(1-m_ab)+NsurvivedB[y,2:(x+1)]*m_ba
   NmovedB[y,1]<-NsurvivedB[y,1]
   NmovedB[y,2:(x+1)]=NsurvivedB[y,2:(x+1)]*(1-m_ba)+NsurvivedA[y,2:(x+1)]*m_ab
  

   #Age the fish
   NstartA[y+1,1]<-Rec_post_A[y]
   NstartA[y+1,x+1]<-NmovedA[y,x]+NmovedA[y,x+1]
   NstartA[y+1,2:x]<-NmovedA[y,1:(x-1)]

   NstartB[y+1,1]<-Rec_post_B[y]
   NstartB[y+1,x+1]<-NmovedB[y,x]+NmovedB[y,x+1]
   NstartB[y+1,2:x]<-NmovedB[y,1:(x-1)]
   
   Biomass<-sum(NstartA[y+1,(am+1):(x+1)]*wt[(am+1):(x+1)])+sum(NstartB[y+1,(am+1):(x+1)]*wt[(am+1):(x+1)])
   CR<-alpha2/((Rec_post_A[y]+Rec_post_B[y])/Biomass)
   
}
    plot((seq(1:(years+1))),(rowSums(NstartA)+rowSums(NstartB)))
     lines((seq(1:(years+1))),(rowSums(NstartA)),col=2)
  lines((seq(1:(years+1))),(rowSums(NstartB)),col=3)
#    print(Biomass)
#    priagxxnt(CR)
#    minopt<-((Biomass-1000000)^2+(CR-4)^6)
#    return(minopt)
#    
#    }
#    
     #optim(c(4,.000002),sim)

