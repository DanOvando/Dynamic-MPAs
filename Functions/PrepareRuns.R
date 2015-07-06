
PrepareRuns<- function(Species,Fs,MPAs,Discs)
{
  
  RunMat<- as.data.frame(matrix(NA,nrow= length(Species)*length(Fs)*length(MPAs)*length(Discs),ncol=4))
  
  colnames(RunMat)<- c('Species','FLevel','ReserveStrategy','DiscountRate')
  
  b<- 0
  
  for (s in 1:length(Species))
  { 
    for (f in 1:length(Fs))
    {  
      for (m in 1:length(MPAs))
      {
        for(d in 1:length(Discs))
        {
          b<- b+1
          RunMat[b,]<- data.frame(Species[s],Fs[f],MPAs[m],Discs[d],stringsAsFactors=F)
        }
      }
    }  
  }
  
  return(RunMat)
}