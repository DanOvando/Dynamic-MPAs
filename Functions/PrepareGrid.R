
PrepareGrid<- function(Species,Fs,ReserveInc,InterceptInc,SlopeInc,Discs)
{
  ReserveSeq<- seq(from=ReserveInc,to=0.6,by=ReserveInc)
  
  InterceptSeq<- seq(from=0,to=1,by=InterceptInc)
  
  SlopeSeq<- seq(from=0,to=1,by=SlopeInc)
  
  RunMat<- as.data.frame(matrix(NA,nrow= length(Species)*length(Fs)*length(ReserveSeq)*length(InterceptSeq)*length(SlopeSeq)*length(Discs),ncol=6))
  
  colnames(RunMat)<- c('Species','FLevel','ReserveSize','Intercept','Slope','DiscountRate')
  
  b<- 0
  
  for (s in 1:length(Species))
  { 
    for (f in 1:length(Fs))
    {  
      for (r in 1:length(ReserveSeq))
      {
        for (i in 1:length(InterceptSeq))
        {
          for (m in 1:length(SlopeSeq))
          {
            for(d in 1:length(Discs))
            {
              b<- b+1
              RunMat[b,]<- data.frame(Species[s],Fs[f],ReserveSeq[r],InterceptSeq[i],SlopeSeq[m],Discs[d],stringsAsFactors=F)
            }
          }
        }  
      }
    }
  }
  return(RunMat)
}