# load('Completed Workspace.rdata')
# 
# 
# ###Not as slick as xyplot, but gives you full control over what is plotted
# 
# plotyears=30  #no. years to plot
# 
# #subset data for each scenario (m)
# m1=subset(PlotStorage,m=="EQ Optimal")
# m2=subset(PlotStorage,m=="SNTZ")
# m3=subset(PlotStorage,m=="GNTZ")
# m4=subset(PlotStorage,m=="OptNTZ")
# 
# #get points where yield equals status quo
# y1=PlotStorage[PlotStorage$f=="F/Fmsy=1.5" & PlotStorage$PositiveYields==1,c(6,18)]
# y2=PlotStorage[PlotStorage$f=="F/Fmsy=3" & PlotStorage$PositiveYields==1,c(6,18)]
# 
# 
#  windows()
# par(mfrow=c(2,1),mar=c(0,0,0,0),oma=c(4,6,2,6))
#  
# plot(m1$Year[as.vector(m1$f)=="F/Fmsy=1.5"],m1$NPB[as.vector(m1$f)=="F/Fmsy=1.5"],
#   type="l",ylab="",xlab="",ylim=c(min(PlotStorage$NPB),max(PlotStorage$NPB)),
#   xlim=c(1,plotyears),col="dodgerblue",lwd=2,las=1,xaxt="n" )
# mtext("F/Fmsy = 1.5",line=-1)
# text(1,max(PlotStorage$NPB)*0.99,"A.")
# lines( m2$Year[as.vector(m2$f)=="F/Fmsy=1.5"],m2$NPB[as.vector(m2$f)=="F/Fmsy=1.5"],col=2,lwd=2,lty=3)
# lines( m3$Year[as.vector(m3$f)=="F/Fmsy=1.5"],m3$NPB[as.vector(m3$f)=="F/Fmsy=1.5"],col="medium orchid",lwd=2,lty=4)
# lines( m4$Year[as.vector(m4$f)=="F/Fmsy=1.5"],m4$NPB[as.vector(m4$f)=="F/Fmsy=1.5"],col="darkgreen",lwd=2,lty=5)
# abline(h=0,lty=2)
# points(y1[,1],y1[,2],col=c("dodgerblue",2,"medium orchid","darkgreen"),pch=16,cex=1.5)
# 
# 
# plot(m1$Year[as.vector(m1$f)=="F/Fmsy=3"],m1$NPB[as.vector(m1$f)=="F/Fmsy=3"],
#   type="l",ylab="",xlab="Year",ylim=c(min(PlotStorage$NPB),max(PlotStorage$NPB)),
#   xlim=c(1,plotyears),col="dodgerblue",lwd=2 ,las=1)
# mtext("F/Fmsy = 3.0",line=-1)
# text(1,max(PlotStorage$NPB)*0.99,"B.")
# lines( m2$Year[as.vector(m2$f)=="F/Fmsy=3"],m2$NPB[as.vector(m2$f)=="F/Fmsy=3"],col=2,lwd=2,lty=3)
# lines( m3$Year[as.vector(m3$f)=="F/Fmsy=3"],m3$NPB[as.vector(m3$f)=="F/Fmsy=3"],col="medium orchid",lwd=2,lty=4)
# lines( m4$Year[as.vector(m4$f)=="F/Fmsy=3"],m4$NPB[as.vector(m4$f)=="F/Fmsy=3"],col="darkgreen",lwd=2,lty=5)
# abline(h=0,lty=2)         
# points(y2[,1],y2[,2],col=c("dodgerblue",2,"medium orchid","darkgreen"),pch=16,cex=1.5)
# 
# par(new=TRUE, mfrow=c(1,1), oma=c(0,0,0,0))
# 
# plot(1,
#   xlim=c(0,1),
#   ylim=c(0,1),
#   xaxs="i",
#   yaxs="i",
#   type="n",
#   axes=FALSE)
# 
# 
# legend(0.8,0.9, legend=c("EQ Opt","SNTZ","GNTZ","Opt NTZ"), lty=c(1,3,4,5),
# col=c("dodgerblue",2,"medium orchid","darkgreen"),lwd=2,bty="n",cex=0.8)
# 
#  text(.5,.05,"Year",font=2)
# 
#  text(.03,.5,"Net Present Balance",srt=90,font=2)
# 
