
# Shrinking Nuts ----------------------------------------------------------

rm(list = ls())
require(stats)
library(proftools)
library(grid)
library(gridExtra)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(broom)
library(RColorBrewer) 
# library(parallel)
library(snowfall)

sapply(list.files(pattern = "[.]R$", path = "Functions", full.names = TRUE), source)

# Set Controls ------------------------------------------------------------

source('BetterNutzControlfile.R')

PopTolerance <- .1 #the cutoff for identifying point where population stops changing
NumCores<- 8
InitialPopulation <- 1000 #Seed for initial population 
CollapseThreshold <- 0.1
LookAtLengths <- 0
ReservePosition <- 'Center'
OptTime <- 10 #Time Horizon to optimize over
Alpha <- 0.5


TimeToRun <- OptTime + 1

DiscRates <- 0.1

BasePatches <- Patches

BatchFolder <- 'Results/Full Grid/'

RunAnalysis <- TRUE


# Run Analysis ---------------------------------------------------------------


if (RunAnalysis == TRUE) {
  
  dir.create(paste(BatchFolder,sep = ''))
  
  LifeHistories <- read.csv('Inputs/Life History Inputs.csv',stringsAsFactors = F)
  
  LifeColumns <- colnames(LifeHistories)
  
  LifeVars <- c('Range','MaxAge','m','k','Linf','t0','AgeMa50','LengthMa50','MaturityMode','BH.Steepness','wa','wb','WeightForm','fa','fb','DDForm')
  
  SpeciesList <- LifeHistories$CommName
  
  #   SpeciesList<- SpeciesList[1:2]
  
  SystemBmsyStorage <- as.data.frame(matrix(NA,nrow=length(SpeciesList),ncol=2))
  
  colnames(SystemBmsyStorage) <- c('Species','Bmsy')
  
  SystemBmsyStorage$Species <- as.character(SystemBmsyStorage$Species)
  
  NumFs <- 1
  
  DefaultLifeHistory <- lh
  
  BasePatches <- Patches
  
  Populations<- sapply(gsub(' ','',SpeciesList,fixed = T),TargetPop=0.25,PreparePopulations,LifeHistories=LifeHistories,BaseLife=lh,LifeVars = LifeVars,BasePatches = BasePatches,BatchFolder = BatchFolder,USE.NAMES = T)
  
  
  RunMatrix <- PrepareGrid(SpeciesList,Fs='Set to 0.25',ReserveInc = 0.25,InterceptInc = 0.25,SlopeInc = 0.25,DiscRates)
  
  save(file = paste(BatchFolder,'run_matrix.Rdata'),RunMatrix)
  
  #       Rprof(tmp <- tempfile(),line.profiling=T)
  
  #   dim(RunMatrix)[1]
  if (file.exists('NutsProgress.txt')) {file.remove('NutsProgress.txt') }

  if (Sys.info()[1]=='Windows')
  {
    
    sfInit( parallel=TRUE, cpus=NumCores)
    
#     sfExport("RunMatrix","BasePatches","Populations","BatchFolder","TimeToRun",local=FALSE)
    
    sfExportAll()
    sfLibrary(dplyr)
    
    ReserveResults <- sfClusterApplyLB(1:dim(RunMatrix)[1], RunGridReserve,RunMatrix=RunMatrix,BasePatches = BasePatches,
                                       Populations = Populations, BatchFolder = BatchFolder,TimeToRun=TimeToRun) %>% ldply()
    sfStop()
  }
  
  if (Sys.info()[1]!='Windows')
  {
    ReserveResults <- (mclapply(1:dim(RunMatrix)[1], RunGridReserve, RunMatrix = RunMatrix,BasePatches = BasePatches,
                                Populations = Populations, BatchFolder = BatchFolder,TimeToRun=TimeToRun, mc.cores = NumCores )) %>% ldply()
  }
  
  save(file = paste(BatchFolder,'Reserve Results.Rdata'),ReserveResults)

    ReserveResults$YieldBalance <- ReserveResults$Yield - ReserveResults$SQYield
  
  #   ReserveResults$ScenId <- with(ReserveResults,paste(Species,m,f,sep = '- '))
  
  ReserveResults <- ReserveResults %>%
    group_by(Run) %>% 
    arrange(Year) %>%
    mutate(PresentYield=Yield*(1+Fleet$YieldDisc)^-(Year-1),PresentBalance=(Yield-SQYield)*(1+Fleet$YieldDisc)^-(Year-1)) %>%
    mutate(NPB=cumsum(PresentBalance),NPY=cumsum(PresentYield),RequestedLoan = sum(PresentBalance[YieldBalance<0]),
           ScaledNPB=NPB/max(NPB,na.rm=T))
  
#   quartz()
#   (ggplot(subset(ReserveResults,Year==max(Year)),aes(Intercept,FinalReserve))+
#     geom_tile(aes(fill=NPB))+facet_wrap(~Species,scales='free')+
#     scale_fill_gradientn(colours=RColorBrewer::brewer.pal(name='RdYlGn',n=9)))
  
  
  save.image(file=paste(BatchFolder,'NutsResults.Rdata',sep=''))
}
if (RunAnalysis==F)
{
  load(file=paste(BatchFolder,'NutsResults.Rdata',sep=''))
}


# Process Results ---------------------------------------------------------



# pdf(file='colorpal.pdf')
# display.brewer.all()
# dev.off()

FontColor<- 'black'


KeynoteTheme<- theme(legend.position='top',plot.background=element_rect(color=NA),rect=element_rect(fill='transparent',color=NA),text=element_text(size=22,family=Font,color=FontColor),
                     axis.text=element_text(color=FontColor),axis.title.y=element_text(size=12,hjust=0.5,angle=0),axis.text.y=element_text(size=12),axis.text.x=element_text(angle=35, vjust=0.9,hjust=0.9,color=FontColor,size=12),
                     legend.text=element_text(size=10,color='black'),legend.background=element_rect(fill="gray90"),legend.title=element_blank())

PaperTheme<- theme(legend.position='top',text=element_text(size=22,family=Font,color=FontColor),
                   axis.text=element_text(color=FontColor),axis.title.y=element_text(size=25,hjust=0.5,angle=0),axis.text.y=element_text(size=30),axis.text.x=element_text(angle=35, vjust=0.9,hjust=0.9,color=FontColor,size=22),
                   legend.text=element_text(size=14,color='black'),legend.title=element_text(size=16,color='black'))


Breaks<- c(min(ReserveResults$YieldBalance),0,max(ReserveResults$YieldBalance))

pdf(file=paste(BatchFolder,'NPB Trajectory.pdf',sep=''))
NPBPlot<- (ggplot(data=subset(ReserveResults,m=='OptNTZ' | m=='EqNTZ'),aes(x=Year,y=NPB,linetype=m))+geom_line(size=1.2)+
             geom_point(aes(fill=YieldBalance),size=3,shape=21) + facet_wrap(~Species,scales='free_y')+
             scale_fill_gradient2(low='red',mid='yellow',high='green',breaks=Breaks)+geom_hline(aes(yintercept=0)))
print(NPBPlot)
dev.off()

# pdf(file=paste(BatchFolder,'Scaled NPB Trajectory.pdf',sep=''))
# ScaledNPBPlot<- (ggplot(data=subset(ReserveResults,m=='OptNTZ' | m=='EqNTZ'),aes(x=Year,y=ScaledNPB,linetype=m))+geom_line(size=1.2)+
#              geom_point(aes(fill=YieldBalance),size=3,shape=21)+facet_wrap(~Species,scale='free_y')+
#              scale_fill_gradient2(low='red',mid='yellow',high='green',breaks=Breaks)+geom_hline(aes(yintercept=0)))
# print(ScaledNPBPlot)
# dev.off()


ResSummary<- ReserveResults %>%
  group_by(Species,m,f) %>%
  summarize(TimeToNPB=which(NPB>=0)[1],
            NegativeYields=Discount(pmin(0,YieldBalance),Fleet$YieldDiscount,length(YieldBalance))$NPV,
            StatusQuoYields=Discount(SQYield*(YieldBalance<=0),Fleet$YieldDiscount,length(YieldBalance))$NPV,
            PriceInc=100*(StatusQuoYields/(StatusQuoYields+NegativeYields)-1),
            AvailableSurplus=Discount(pmax(0,YieldBalance),Fleet$YieldDiscount,length(YieldBalance))$NPV) %>%
  ungroup() %>%
  mutate(RunName=paste(Species,m,f,sep='-'))


RunNames<- unique(ResSummary$RunName)
for (i in seq_len(length(RunNames)))
{
  Where<- ResSummary$RunName==RunNames[i]
  ResSummary$MaxInterestRate[i]<- 100*exp(optim(-4,FindMaxInterestRate,LoanTime=10,LoanAmount=-ResSummary$NegativeYields[Where],Surplus=ResSummary$AvailableSurplus[Where],lower=-10,upper=10,method='Brent')$par)
}


pdf(file=paste(BatchFolder,'Value Gain Needed.pdf',sep=''))
PriceGainPlot<- (ggplot(ResSummary,aes(x=Species,y=PriceInc,fill=m))+
                   geom_bar(stat='identity',position='dodge')+KeynoteTheme+theme(legend.title=element_blank()))
print(PriceGainPlot)
dev.off()

pdf(file=paste(BatchFolder,'Time to Positive NPB.pdf',sep=''))
TimeToNPBPlot<- (ggplot(ResSummary,aes(x=Species,y=TimeToNPB,fill=m))+
                   geom_bar(stat='identity',position='dodge')+KeynoteTheme)
print(TimeToNPBPlot)
dev.off()

pdf(file= paste(BatchFolder,'Max Interest Rate.pdf',sep= ''))
MaxInterestRatePlot<- (ggplot(ResSummary,aes(x= Species,y=MaxInterestRate,fill=m))+
                         geom_bar(stat='identity',position='dodge')+KeynoteTheme)
print(MaxInterestRatePlot)
dev.off()

