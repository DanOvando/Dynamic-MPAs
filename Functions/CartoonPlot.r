CartoonPlot<- function(Theme)
{

r=0.24
K=100
Nvec=NA
Nvec[1]=K
u1=0.16
u2=0.12
yvec=NA

for(t in 2:150){
  if(t<100)u=u1 else u=u2
  Nvec[t]=Nvec[t-1]+r*Nvec[t-1]*(1-Nvec[t-1]/K)-u*Nvec[t-1]
  yvec[t-1]=u*Nvec[t-1]

  }



Years=1:149
Nvec=Nvec[1:149]
Yield=(yvec-yvec[98])/max(yvec)                 #creat standardized deviations from status quo
Biomass=(Nvec-Nvec[98])/max(Nvec)
dat=data.frame(Years,Biomass,Yield)

datA=subset(dat,Years>90 & Years<=149)
 datA$Years=-7:51

 datA$Balance=cumsum(datA$Yield)

#find year balance is positive for shaded area cutoff
test=subset(datA,Years>0)
test$flag=0
test$flag[test$Balance>0]=1
test$dup=duplicated(test$flag)
breakeven=which(test$dup==FALSE)[2]
breakYear=test$Years[breakeven]



#
#long format
longdat <- reshape(datA[,1:3],
  varying = c("Biomass","Yield"),
  v.names = "value",
  timevar = "Objectives",
  times = c("Biomass","Yield"),
  direction = "long")




p <- ggplot(longdat, aes(x=Years,y=value))

pp <- p +
  geom_area(data=subset(longdat,Years<=12 & Objectives=="Yield"),
              fill="grey47") +
  geom_area(data=subset(longdat,Years>=12 &Years<=breakYear& Objectives=="Yield"),
                                            fill="grey74")+ geom_line(aes(linetype=Objectives),size=1.1) +
  geom_vline(xintercept = 30,linetype=2) +
  ylab("Deviation from Status Quo")  +
  scale_x_discrete(breaks=c(0,30),labels=c("0","Yield Balance")) +
  scale_y_continuous(limits=c(-0.1, 0.2),breaks=0) +
  Theme


# p + geom_area(data=subset(longdat,Years<=12 & Objectives=="Yield"),
#               fill="darkorange")+ geom_area(data=subset(longdat,Years>=12 &Years<=breakYear& Objectives=="Yield"),
#                                             fill="cornflowerblue")+ geom_line(aes(colour=Objectives),size=1.1)+
#   geom_vline(xintercept = 30,linetype=2) + scale_color_manual(values=c("aquamarine3", "black"))  + ylab("Deviation from Status Quo")  + scale_x_discrete(breaks=c(0,30),labels=c("0","Yield Balance"))     + scale_y_continuous(limits=c(-0.1, 0.2),breaks=0)

}
