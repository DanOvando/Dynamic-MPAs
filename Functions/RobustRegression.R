RobustRegression<- function(Model,Data,ClusterVar)
{
  ##Created by Dan Ovando
  ## April 14 2015

  ## This function calculates heteroskedastic robust standard errors in the manner of the 'robust'
  #option in STATA. Robust SEs and associated Variance-Covariance matrix are stored in the model
  # Tidy and summary versions of results are produced and stored
  #Outputs:
  # Model: The regression object
  # TidyModel: A cleaned up version of the mode
  # GlanceModel: summary statistics of the model
  # AugModel: Data with augmented statistics


  library(sandwich,quietly=T)
  library(lmtest,quietly=T)
  library(broom,quietly=T)
  library(tidyr,quietly=T)
  library(dplyr,quietly=T)
  library(zoo,quietly=T)

  Model$VCOV<- vcovHC(Model,type='HC1')

  if (ClusterVar !='None')
  {
    Model$VCOV<- ClusteredVCOV(Model,data=Data,cluster=ClusterVar)
  }


  SEs<- data.frame(t(sqrt(diag(Model$VCOV))),stringsAsFactors=F) %>% gather('variable','RobustSE')

  SEs$variable<- as.character(SEs$variable)

  SEs$variable[SEs$variable=='X.Intercept.']<- '(Intercept)'

  Model$RobustSEs<- SEs

  RobustTest<- (coeftest(Model,vcov.=Model$VCOV))

  StatNames<- colnames(RobustTest)

  VarNames<- rownames(RobustTest)

  RobustMat<- as.data.frame(matrix(NA,nrow=length(VarNames),ncol=2))

  colnames(RobustMat)<- c('variable','RobustPvalue')

  for (i in 1:length(VarNames))
  {
    RobustMat[i,]<- data.frame(as.character(VarNames[i]),RobustTest[i,'Pr(>|t|)'],stringsAsFactors=F)
  }
  TidyModel<- tidy(Model) %>% dplyr::rename(variable=term) %>% join(SEs,by='variable') %>% join(RobustMat,by='variable')

  AugModel<- augment(Model)

  GlanceModel<- glance(Model)

  TidyModel$variable<- as.factor(TidyModel$variable)

  TidyModel$variable <- reorder(TidyModel$variable, TidyModel$RobustPval)

  TidyModel$ShortPval<- pmin(TidyModel$RobustPval,0.2)


  RegPlot<- (ggplot(data=TidyModel,aes(x=variable,y=estimate,fill=ShortPval))+
               geom_bar(position='dodge',stat='identity',color='black')+
               scale_fill_gradient2(high='black',mid='gray99',low='red',midpoint=0.1,
                                    breaks=c(0.05,0.1,0.15,0.2),labels=c('0.05','0.10','0.15','>0.20')
                                    ,name='P-Value',guide=guide_colorbar(reverse=T))
             +theme(axis.text.x=element_text(angle=45,hjust=0.9,vjust=0.9))+
               geom_errorbar(mapping=aes(ymin=estimate-1.96*RobustSE,ymax=estimate+1.96*RobustSE))+
               xlab('Variable')+
               ylab(paste('Marginal Effect on ',names(Model$model)[1],sep='')))

  TidyModel$ShortPval<- NULL

  TCrit<-(qt(c(0.025,0.975),df=Model$df.residual)[2])

  TidyModel$LCI95<- TidyModel$estimate-TCrit*TidyModel$RobustSE

  TidyModel$UCI95<- TidyModel$estimate+TCrit*TidyModel$RobustSE

  return(list(Model=Model,TidyModel=TidyModel,AugModel=AugModel,GlanceModel=GlanceModel,RegPlot=RegPlot))

}