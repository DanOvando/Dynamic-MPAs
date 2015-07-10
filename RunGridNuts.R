
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
library(parallel)
library(snowfall)
library(scales)

sapply(list.files(pattern = "[.]R$", path = "Functions", full.names = TRUE), source)

# Set Controls ------------------------------------------------------------

source('BetterNutzControlfile.R')

PopTolerance <- 1 #the cutoff for identifying point where population stops changing
NumCores<- 2
InitialPopulation <- 1000 #Seed for initial population 
CollapseThreshold <- 0.1
LookAtLengths <- 0
ReservePosition <- 'Center'
OptTime <- 30 #Time Horizon to optimize over
Alpha <- 0.5
fig_width <- 7
fig_height <- 5
text_size <- 10
TimeToRun <- OptTime + 1

DiscRates <- 0.1

BasePatches <- Patches

BatchFolder <- 'Results/Full Grid 2.0/'

RunAnalysis <- FALSE

OptMode <- 'Utility'

Scale_Yields <- T

long_term_objective <- 0 # 0 to 1, 0 means optimal reserve defined by NPB, 1 by total yields (i.e. discount = 0)

discount_rate <- .1

fixed_slope <- 0.2

# Run Analysis ---------------------------------------------------------------


if (RunAnalysis == TRUE) {
  
  dir.create(paste(BatchFolder,sep = ''))
  
  LifeHistories <- read.csv('Inputs/Life History Inputs.csv',stringsAsFactors = F)
  
  LifeColumns <- colnames(LifeHistories)
  
  LifeVars <- c('Range','MaxAge','m','k','Linf','t0','AgeMa50','LengthMa50','MaturityMode','BH.Steepness','wa','wb','WeightForm','fa','fb','DDForm')
  
  SpeciesList <- LifeHistories$CommName
  
  #     SpeciesList<- SpeciesList[1:2]
  
  SystemBmsyStorage <- as.data.frame(matrix(NA,nrow=length(SpeciesList),ncol=2))
  
  colnames(SystemBmsyStorage) <- c('Species','Bmsy')
  
  SystemBmsyStorage$Species <- as.character(SystemBmsyStorage$Species)
  
  NumFs <- 1
  
  DefaultLifeHistory <- lh
  
  BasePatches <- Patches
  
  
  Populations<- sapply(gsub(' ','',SpeciesList,fixed = T),TargetPop=0.25,PreparePopulations,LifeHistories=LifeHistories,BaseLife=lh,LifeVars = LifeVars,BasePatches = BasePatches,BatchFolder = BatchFolder,USE.NAMES = T)
  
  #   Rprof(tmp <- tempfile(),line.profiling=T)
  #   
  #   BatTargetPopulation<- GrowPopulation(100, 0,'EQ',0,'BatTarget Run',Species='Yellowtail Snapper',lh=Populations$YellowtailSnapper,Patches=Patches,FigureFolder=FigureFolder)
  #   
  #       Rprof() 
  #       summaryRprof(tmp)
  #       unlink(tmp)
  
  RunMatrix <- PrepareGrid(SpeciesList,Fs='Set to 0.25',ReserveInc = 0.05,InterceptInc = 0.1,SlopeInc = 0.1,DiscRates)
  
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
    Rprof(tmp <- tempfile(),line.profiling=T)
    
    ReserveResults <- (mclapply(1:dim(RunMatrix)[1], RunGridReserve, RunMatrix = RunMatrix,BasePatches = BasePatches,
                                Populations = Populations, BatchFolder = BatchFolder,TimeToRun=TimeToRun, mc.cores = NumCores )) %>% ldply()
    Rprof() 
    summaryRprof(tmp)
    unlink(tmp)
    
  }
  
  save(file = paste(BatchFolder,'Reserve Results.Rdata'),ReserveResults)
  
  save.image(file=paste(BatchFolder,'NutsResults.Rdata',sep=''))
}
if (RunAnalysis==F)
{
  load(file = paste(BatchFolder,'Reserve Results.Rdata'))
}

# Prepare Figure Themes ---------------------------------------------------
# text_size <- 9


FontColor <- 'black'


KeynoteTheme <- theme(legend.position = 'top', plot.background = element_rect(color = NA), rect=element_rect(fill = 'transparent', color = NA),
                      text = element_text(size = 22,family = Font, color= FontColor),
                      axis.text = element_text(color = FontColor),axis.title.y = element_text(size = 12,hjust= 0.5,angle = 0),
                      axis.text.y = element_text(size = 12),axis.text.x = element_text(angle = 35, vjust = 0.9,hjust = 0.9,color= FontColor,size = 12),
                      legend.text = element_text(size = 10,color= 'black'),legend.background=element_rect(fill = "gray90"),legend.title=element_blank())

PaperTheme <- theme(legend.position = 'top',text = element_text(size= 22,family= Font,color= FontColor),
                    axis.text = element_text(color= FontColor),axis.title.y=element_text(size = 25,hjust = 0.5,angle = 0)
                    ,axis.text.y = element_text( size = 30), axis.text.x = element_text(angle = 35, vjust = 0.9,hjust = 0.9,color = FontColor, size = 22),
                    legend.text = element_text(size = 14,color = 'black'),legend.title = element_text(size = 16,color = 'black'))

SimpleTheme <- theme(legend.position = 'top',text = element_text(size = text_size,family = Font,color = FontColor),
                     panel.background =  element_rect(fill = "white", colour = NA),
                     panel.border =      element_rect(fill = NA, colour = "grey50"), 
                     panel.grid.major =  element_line(colour = "grey60", size = 0.2),
                     panel.grid.minor =  element_line(colour = "grey78", size = 0.5))

# Process Results ---------------------------------------------------------

Fleet$YieldDiscount <- discount_rate

if (Scale_Yields == T)
{
  YieldFactor<- ReserveResults$Yield/ReserveResults$SQYield
  
  ReserveResults$SQYield <- 100
  
  ReserveResults$Yield <- ReserveResults$SQYield * YieldFactor
  
  ReserveResults$YieldBalance <- ReserveResults$Yield - ReserveResults$SQYield
  
}

ReserveResults$YieldBalance <- ReserveResults$Yield - ReserveResults$SQYield

ReserveResults <- ReserveResults %>%
  group_by(Run) %>% 
  arrange(Year) %>%
  mutate(PresentYield=Yield*(1+Fleet$YieldDisc)^-(Year-1),PresentBalance=(Yield-SQYield)*(1+Fleet$YieldDisc)^-(Year-1)) %>%
  mutate(NPB=cumsum(PresentBalance),NPY=cumsum(PresentYield),RequestedLoan = sum(PresentBalance[YieldBalance<0])) %>%
  ungroup() %>%
  group_by(Species) %>% 
  mutate(s_Yield = Yield/max(Yield,na.rm=T),s_PresentYield=s_Yield*(1+Fleet$YieldDisc)^-(Year-1)) %>%
  ungroup() %>%
  group_by(Run) %>%
  arrange(Year) %>%
  mutate(net_yields = cumsum(Yield), s_NPY=cumsum(s_PresentYield),s_Balance=(Yield/SQYield),
         s_NPB=cumsum(s_Balance*(1+Fleet$YieldDisc)^-(Year-1))) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(Utility = long_term_objective * (net_yields / max(net_yields)) + (1 - long_term_objective) * (NPY/max(NPY)) ) %>%
  ungroup()



OptimalRun<- filter(ReserveResults,Year==max(Year)) %>%
  group_by(Species) %>%
  summarize(OptNPB = Run[NPB == max(NPB)], OptEQ = Run[net_yields == max(net_yields)], OptUtility = Run[Utility == max(Utility)])

ReserveResults$Scenario <- with(ReserveResults,paste(FinalReserve,Intercept,Slope,sep = '-'))

UnifiedRuns<- filter(ReserveResults,Year==max(Year)) %>%
  group_by(Scenario) %>%
  summarize(TotalUtility = sum(Utility))

BestUnifiedScenario <- UnifiedRuns$Scenario[UnifiedRuns$TotalUtility == max(UnifiedRuns$TotalUtility)]

ReserveResults$BestUnifiedRun <- ReserveResults$Scenario %in% BestUnifiedScenario


if (OptMode == 'NPB')
{
  ReserveResults$BestRun <- ReserveResults$Run %in% OptimalRun$OptNPB
}

if (OptMode == 'Yields')
{
  ReserveResults$BestRun <- ReserveResults$Run %in% OptimalRun$OptEQ
}

if (OptMode == 'Utility')
{
  ReserveResults$BestRun <- ReserveResults$Run %in% OptimalRun$OptUtility
}

ResSummary <- ReserveResults %>%
  group_by(Run) %>%
  summarize(TimeToNPB = which(NPB >= 0)[1],
            Species = unique(Species), 
            BestRun = unique(BestRun), 
            FinalNPB = last(NPB), 
            ReserveSize = mean(FinalReserve), 
            Intercept = mean(Intercept), 
            Slope = mean(Slope),
            NegativeYields = Discount(pmin(0,YieldBalance),Fleet$YieldDiscount,length(YieldBalance))$NPV,
            StatusQuoYields = Discount(SQYield * (YieldBalance <= 0),Fleet$YieldDiscount,length(YieldBalance))$NPV,
            PriceInc = 100 * (StatusQuoYields / (StatusQuoYields + NegativeYields) - 1),
            AvailableSurplus = Discount(pmax(0, YieldBalance),Fleet$YieldDiscount,length(YieldBalance))$NPV)

RunNames <- unique(ResSummary$Run)
for (i in seq_len(length(RunNames)))
{
  Where <- ResSummary$Run ==  RunNames[i]
  ResSummary$MaxInterestRate[i] <- 100 * exp( optim(-4, FindMaxInterestRate, LoanTime = 10, LoanAmount = -ResSummary$NegativeYields[Where],
                                                    Surplus = ResSummary$AvailableSurplus[Where],lower = -10, upper = 10,method = 'Brent')$par)
}

ResSummary$TimeToNPB[is.na(ResSummary$TimeToNPB)] <- max(ReserveResults$Year)

# Analyze Optimal Reserve -------------------------------------------------

OptRun<- filter(ReserveResults,BestRun == T)

opt_npb_plot<- (ggplot(data = OptRun ,aes(x = Year, y = NPB))+geom_line(size = 1.2)+
                  geom_point(aes(fill = s_Balance-1,size = CurrentReserve), shape = 21) + 
                  facet_wrap(~Species,scales = 'free_y') +
                  scale_fill_gradient2(labels = percent, name = '% Change from SQ Yields',low = 'red', mid = 'yellow', high = 'green',midpoint = 0)
                + geom_hline(aes(yintercept = 0), linetype = 'longdash', color = 'grey2') + 
                  SimpleTheme + 
                  scale_size_continuous(labels = percent, name = '% Reserve') + 
                  theme(legend.position = 'right'))

ggsave(file=paste(BatchFolder,'Optimal NPB Trajectory.pdf',sep=''),plot=opt_npb_plot,width=fig_width,height=fig_height)



opt_reserve_plot<- (ggplot(data = OptRun ,aes(x = Year, y = CurrentReserve))+geom_line(size = 1.2)+
                      geom_point(aes(fill = NPB), shape = 21) + 
                      facet_wrap(~Species,scales  = 'fixed') +
                      scale_fill_gradient2(name = 'NPB',low = 'red', mid = 'yellow', high = 'green',midpoint = 1) + 
                      SimpleTheme + 
                      geom_hline(aes(yintercept = 0, linetype = 'longdash')) + 
                      xlab('Year') + 
                      scale_y_continuous(labels = percent) + 
                      theme(legend.position = 'right'))
ggsave(file=paste(BatchFolder,'Optimal Reserve Trajectory.pdf',sep=''),plot=opt_reserve_plot,width=fig_width,height=fig_height)


unified_npb_plot<- (ggplot(data = filter(ReserveResults,BestUnifiedRun == T) ,aes(x = Year, y = NPB))+geom_line(size = 1.2)+
                  geom_point(aes(fill = s_Balance-1,size = CurrentReserve), shape = 21) + 
                  facet_wrap(~Species,scales = 'free_y') +
                  scale_fill_gradient2(labels = percent, name = '% Change from SQ Yields',low = 'red', mid = 'yellow', high = 'green',midpoint = 0)
                + geom_hline(aes(yintercept = 0), linetype = 'longdash', color = 'grey2') + 
                  SimpleTheme + 
                  scale_size_continuous(labels = percent, name = '% Reserve') + 
                  theme(legend.position = 'right'))

ggsave(file=paste(BatchFolder,'Unified NPB Trajectory.pdf',sep=''),plot=unified_npb_plot,width=fig_width,height=fig_height)


unified_reserve_plot<- (ggplot(data = filter(ReserveResults,BestUnifiedRun == T) ,aes(x = Year, y = CurrentReserve))+geom_line(size = 1.2)+
                      geom_point(aes(fill = NPB), shape = 21) + 
                      facet_wrap(~Species,scales  = 'fixed') +
                      scale_fill_gradient2(name = 'NPB',low = 'red', mid = 'yellow', high = 'green',midpoint = 1) + 
                      SimpleTheme + 
                      geom_hline(aes(yintercept = 0, linetype = 'longdash')) + 
                      xlab('Year') + 
                      scale_y_continuous(labels = percent) + 
                      theme(legend.position = 'right'))
ggsave(file=paste(BatchFolder,'Unified Reserve Trajectory.pdf',sep=''),plot=unified_reserve_plot,width=fig_width,height=fig_height)



SpeciesList<- unique(ReserveResults$Species)

species_comparison<- as.data.frame(matrix(NA,nrow=length(SpeciesList)*length(SpeciesList),ncol=6))

colnames(species_comparison)<- c('Species1','Species2','NPB','TimeToNPB','PrinceInc','MaxInterestRate')

cc <- 0

ResSummary$RunDescription<- paste(ResSummary$ReserveSize,ResSummary$Intercept,ResSummary$Slope,sep='-')

for (i in 1:length(SpeciesList))
{
  Where<- ResSummary$Species == SpeciesList[i] & ResSummary$BestRun == T
  
  species_opt_res<- ResSummary$RunDescription[Where]
  
  for (j in 1:length(SpeciesList))
  {
    cc <- cc+1
    
    Where2<- ResSummary$Species == SpeciesList[j] & ResSummary$RunDescription == species_opt_res
    
    species_comparison[cc,]<- data.frame(SpeciesList[i], SpeciesList[j], 
                                         ResSummary$FinalNPB[Where2],
                                         ResSummary$TimeToNPB[Where2], 
                                         ResSummary$PriceInc[Where2], ResSummary$MaxInterestRate[Where2],stringsAsFactors = F)
    
  }
}

species_comparison<- species_comparison %>%
  ungroup() %>%
  mutate(pos_NPB = NPB + max(NPB)) %>%
  group_by(Species2) %>%
  mutate(Percent_of_best_NPB = sign(max(NPB)) * (NPB / max(NPB) - 1 ) ) %>%
  arrange(Species2)

species_comp_plot<- (ggplot(data = species_comparison,aes(Species1,Species2,fill = Percent_of_best_NPB),color='white') + 
                       geom_tile(color = 'black') + 
                       scale_fill_gradient(name = '% of Opt. NPB',labels = percent, low = 'red', high = 'green') + 
                       theme_classic() + 
                       xlab('Reserve Used') + 
                       ylab('Species Managed') + 
                       theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.x = element_text(angle = 35, vjust = 0.9,hjust = 0.9)))

ggsave(file = paste(BatchFolder,'Opt Species Comparison.pdf',sep = ''),plot = species_comp_plot,width = fig_width,height = fig_height)

species_compnoshark_plot<- (ggplot(data = filter(species_comparison,Species2 != 'Blacknose Shark' & Species2 != 'Silver Trevally'),aes(Species1,Species2,fill = Percent_of_best_NPB)) + 
                              geom_tile(color = 'black') + 
                              scale_fill_gradient(name = '% of Opt. NPB',labels = percent, low = 'red', high = 'green') + 
                              theme_classic() + 
                              xlab('Reserve Used') + 
                              ylab('Species Managed') + 
                              theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.x = element_text(angle = 35, vjust = 0.9,hjust = 0.9)))

ggsave(file = paste(BatchFolder,'Opt Species Comparison no shark.pdf',sep = ''),plot = species_compnoshark_plot,width = fig_width,height = fig_height)



species_comp_bar_plot<- (ggplot(data = species_comparison,aes(x = Species1,y = Percent_of_best_NPB),color='white') + 
                           geom_bar(stat = 'identity', position = 'dodge') + 
                           facet_wrap(~Species2) +
                           scale_y_continuous(labels = percent) + 
                           SimpleTheme + 
                           geom_hline(aes(yintercept = 0)) + 
                           xlab('Reserve Used') + 
                           ylab('Species Managed') + 
                           theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.x = element_text(angle = 35, vjust = 0.9,hjust = 0.9)))

ggsave(file = paste(BatchFolder,'Species Comparison Barchart.pdf',sep = ''),plot = species_comp_bar_plot,width = fig_width,height = fig_height)



# Analyze Static Reserve -----------------------------------------

StaticReserve <- filter(ReserveResults, Year == max(Year) & Intercept == 0 & Slope == 0)

StaticSummary <- filter(ResSummary, Intercept == 0 & Slope == 0)

static_NPB_plot <- (ggplot(data = StaticReserve  ,
                           aes(x = FinalReserve,y = NPB,fill = NPB)) +
                      geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
                      facet_wrap(~Species,scales = 'fixed') + 
                      geom_hline(aes(yintercept = 0), linetype = 'longdash') + 
                      scale_fill_gradient2(low = 'red', mid = 'yellow', high = 'darkgreen', midpoint = 0, name = 'Relative Yields') +
                      #                  scale_fill_brewer(palette = 'Dark2') +
                      xlab('Size of Reserve') + 
                      scale_x_continuous(labels = percent) + 
                      SimpleTheme + 
                      theme(legend.position = 'none'))
ggsave(file = paste(BatchFolder,'Static NPB.pdf',sep = ''),plot = static_NPB_plot,width = fig_width,height = fig_height)

EQ_yieldbalance_plot <- (ggplot(data = StaticReserve,
                                aes(x = FinalReserve,y = s_Balance,fill = s_Balance)) +
                           geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
                           facet_wrap(~Species,scales = 'fixed') + 
                           geom_hline(aes(yintercept = 1), linetype = 'longdash') + 
                           scale_fill_gradientn(colours = RColorBrewer::brewer.pal(name = 'RdYlGn', n = 8)) +
                           #                  scale_fill_brewer(palette = 'Dark2') +
                           xlab('Size of Reserve') + 
                           ylab('% of SQ Yields') + 
                           scale_x_continuous(labels = percent) + 
                           scale_y_continuous(labels = percent) + 
                           SimpleTheme + 
                           theme(legend.position = 'none'))
ggsave(file = paste(BatchFolder,'EQ Yield Balance.pdf',sep = ''),plot = EQ_yieldbalance_plot,width = fig_width,height = fig_height)


static_netbenefit_plot <- (ggplot(data = StaticSummary  ,
                                  aes(x = ReserveSize,y = TimeToNPB,fill = TimeToNPB)) +
                             geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
                             facet_wrap(~Species,scales = 'free_y') + 
                             geom_hline(aes(yintercept = 0), linetype = 'longdash') + 
                             scale_fill_gradient(low = 'green', high = 'red') +
                             #                  scale_fill_brewer(palette = 'Dark2') +
                             xlab('Size of Reserve') + 
                             ylab('Years to Net Benefit') + 
                             scale_x_continuous(labels = percent) + 
                             SimpleTheme + 
                             theme(legend.position = 'none'))
ggsave(file = paste(BatchFolder,'Static Years to Net Benefit.pdf',sep = ''),plot = static_netbenefit_plot,width = fig_width,height = fig_height)


static__priceinc_plot <- (ggplot(data = StaticSummary  ,
                                 aes(x = ReserveSize,y = PriceInc/100,fill = PriceInc)) +
                            geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
                            facet_wrap(~Species, scales = 'fixed') + 
                            geom_hline(aes(yintercept = 0), linetype = 'longdash') + 
                            scale_fill_gradient(low = 'green', high = 'red') +
                            #                  scale_fill_brewer(palette = 'Dark2') +
                            xlab('Size of Reserve') + 
                            ylab('Price Increase Needed') + 
                            scale_x_continuous(labels = percent) + 
                            scale_y_continuous(labels = percent) + 
                            SimpleTheme + 
                            theme(legend.position = 'none'))
ggsave(file = paste(BatchFolder,'Static Price Increase Needed.pdf',sep = ''),plot = static__priceinc_plot,width = fig_width,height = fig_height)


static__maxinterest_plot <- (ggplot(data = StaticSummary  ,
                                    aes(x = ReserveSize,y = MaxInterestRate/100,fill = MaxInterestRate)) +
                               geom_bar(stat = 'identity', position = 'dodge', color = 'black') +
                               facet_wrap(~Species, scales = 'fixed') + 
                               scale_fill_gradient(high = 'green', low = 'red') +
                               #                  scale_fill_brewer(palette = 'Dark2') +
                               xlab('Size of Reserve') + 
                               ylab('Maximum Interest Rate') + 
                               scale_x_continuous(labels = percent) + 
                               scale_y_continuous(labels = percent) + 
                               SimpleTheme + 
                               theme(legend.position = 'none'))
ggsave(file = paste(BatchFolder,'Static Interest Rate Possible.pdf',sep = ''),plot = static__maxinterest_plot,width = fig_width,height = fig_height)

# Analyze Dynamic Reserves ------------------------------------------------

surface_results <- subset(ReserveResults,Year == max(Year) & Slope == fixed_slope)

surface_summary <- subset(ResSummary, Slope == fixed_slope)   

yield_surface_plot<- (ggplot(surface_results,aes(Intercept,FinalReserve))+
                        geom_raster(aes(fill=net_yields),interpolate = T)+facet_wrap(~Species,scales='fixed')+
                        scale_fill_gradient(name = 'Cumulative Yields',low = 'red', high = 'green') + 
                        scale_y_continuous(labels = percent, name = '% in Reserve') + 
                        SimpleTheme + 
                        theme(axis.line = element_blank())
)
ggsave(file=paste(BatchFolder,'Cumulative Yield Surface.pdf',sep=''),plot = yield_surface_plot, height = fig_height, width = fig_width)

npb_surface_plot<- (ggplot(surface_results,aes(Intercept,FinalReserve))+
                      geom_raster(aes(fill=NPB),interpolate = T)+facet_wrap(~Species,scales='fixed')+
                      scale_fill_gradient2(name = 'NPB',low = 'red', mid = 'yellow', high = 'green', midpoint = 0) + 
                      scale_y_continuous(labels = percent, name = '% in Reserve') + 
                      SimpleTheme + 
                      theme(axis.line = element_blank())
)
ggsave(file=paste(BatchFolder,'NPB Surface.pdf',sep=''),plot = npb_surface_plot, height = fig_height, width = fig_width)

priceinc_surface_plot<- (ggplot(surface_summary,aes(Intercept,ReserveSize))+
                           geom_raster(aes(fill=pmin(100,PriceInc)/100),interpolate = T)+facet_wrap(~Species,scales='fixed')+
                           scale_fill_gradient2(name = '% Price Gain Needed', labels = percent,low = 'green', mid = 'yellow', high = 'red', midpoint = .25) + 
                           scale_y_continuous(labels = percent, name = '% in Reserve') + 
                           SimpleTheme + 
                           theme(axis.line = element_blank())
)
ggsave(file=paste(BatchFolder,'Prince Inc Surface.pdf',sep=''),plot = priceinc_surface_plot, height = fig_height, width = fig_width)


surface_summary$LoanType[surface_summary$MaxInterestRate <= 5] <- 'Philanthropy' 

surface_summary$LoanType[surface_summary$MaxInterestRate > 5 & surface_summary$MaxInterestRate <=15] <- 'Social Investment' 

surface_summary$LoanType[surface_summary$MaxInterestRate > 15] <- 'Commercial' 

loan_surface_plot<- (ggplot(surface_summary,aes(Intercept,ReserveSize))+
                       geom_raster(aes(fill=pmin(100,MaxInterestRate)/100),interpolate = T)+facet_wrap(~Species,scales='fixed')+
                       scale_fill_gradient2(name = 'Max Loan Rate', labels = percent,low = 'red', mid = 'yellow', high = 'green', midpoint = .25) + 
                       scale_y_continuous(labels = percent, name = '% in Reserve') + 
                       SimpleTheme + 
                       theme(axis.line = element_blank())
)
ggsave(file=paste(BatchFolder,'Loan Surface.pdf',sep=''),plot = loan_surface_plot, height = fig_height, width = fig_width)


loantype_surface_plot<- (ggplot(surface_summary,aes(Intercept,ReserveSize))+
                           geom_raster(aes(fill=LoanType))+facet_wrap(~Species,scales='fixed')+
                           scale_fill_manual(values = c('green','red','steelblue2'),name = element_blank()) + 
                           scale_y_continuous(labels = percent, name = '% in Reserve') + 
                           SimpleTheme + 
                           theme(axis.line = element_blank())
)
ggsave(file=paste(BatchFolder,'Loan Type Surface.pdf',sep=''),plot = loantype_surface_plot, height = fig_height, width = fig_width)


NPB_model<- RobustRegression(lm(NPB ~ FinalReserve+Intercept+Slope+as.factor(Species),data=filter(ReserveResults,Year==max(Year))),
                             Data=filter(ReserveResults,Year==max(Year)),ClusterVar = 'None')

save(file = paste(BatchFolder,'NutsPlots.Rdata',sep = ''),list=ls(pattern = '_plot'))

