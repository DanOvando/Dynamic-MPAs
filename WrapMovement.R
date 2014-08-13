

P<-30                         #Number of homogeneous patches, differ only in closures

##################################
#Movement
##################################

#dispersal:common larval pool
	DL1=matrix(1/(P),nrow=P,ncol=P)

#dispersal:gaussian movement
	#create movement probability matrix, start with distance matrix
	loc<-seq(-(P-1),2*P,1)
	area<-seq(1,P,1)

	#create distance matrix
	dist<-matrix(NA,nrow=P,ncol=P*3)
	for(i in 1:P)
		{
			for(j in 1:(P*3))
			{
			dist[i,j]<-area[i]-loc[j]
			}
		}
	#now create the movement matrix of probabilities of movement
	p.init<-round(exp(-((dist)^2)/(2*(sigmaL^2))),2)

	#now add matrices on ends to wrap movement and normalize so movement from any one area sums to one
	p.all<-matrix(NA,nrow=P,ncol=(3*P))
	for(i in 1:P)
	{
		for(j in 1:(3*P))
		{
			p.all[i,j]<-(p.init[i,j])/sum(p.init[i,])
		}
	}

	p1<-p.all[,1:P]
	p2<-p.all[,(2*P+1):(3*P)]
	parea<-p.all[,(P+1):(2*P)]
	DL2<-p1+p2+parea




#calculates larval dispersal among patches, so far can choose between common larval pool (DL1) and gaussian dispersal (DL2)
settle<-function(EggV,Disp=DL2){
        settlers=EggV%*%Disp
        return(settlers)
        }
