
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
NumCores<- 1
InitialPopulation <- 1000 #Seed for initial population 
CollapseThreshold <- 0.1
LookAtLengths <- 0
ReservePosition <- 'Center'
OptTime <- 50 #Time Horizon to optimize over
Alpha <- 0.5

Font <- 'Helvetica'
fig_width <- 6
fig_height <- 4
text_size <- 12
TimeToRun <- OptTime + 1
LoanTime <- 11
DiscRates <- 0.1

BasePatches <- Patches

BatchFolder <- 'Results/Full Grid 2.2 BvBmsy _15/'

RunAnalysis <- FALSE

OptMode <- 'Utility'

Scale_Yields <- T

npb_focus <- .9 # 0 to 1, 0 means optimal reserve defined by NPB, 1 by total yields (i.e. discount = 0)

discount_rate <- .22

fixed_slope <- 0.2

# Run Analysis ---------------------------------------------------------------

LifeHistories <- read.csv('Inputs/Life History Inputs.csv',stringsAsFactors = F)

LifeColumns <- colnames(LifeHistories)

LifeVars <- c('Range','MaxAge','m','k','Linf','t0','AgeMa50','LengthMa50','MaturityMode','BH.Steepness','wa','wb','WeightForm','fa','fb','DDForm')

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
                     axis.text.x = element_text(size = 9, angle = 35, vjust = 0.9, hjust = 0.9),
                     panel.background =  element_rect(fill = "white", colour = NA),
                     panel.border =      element_rect(fill = NA, colour = "grey50"), 
                     panel.grid.major =  element_line(colour = "grey60", size = 0.2),
                     panel.grid.minor =  element_line(colour = "grey78", size = 0.5))

# Process Results ---------------------------------------------------------

Fleet$YieldDiscount <- discount_rate

ProcessedNuts <- ProcessNuts(ReserveResults = ReserveResults, Fleet = Fleet, Scale_Yields = Scale_Yields, npb_focus = npb_focus)

ReserveResults <- ProcessedNuts$ReserveResults

species_comparison <- ProcessedNuts$species_comparison

ResSummary <- ProcessedNuts$ResSummary

# Analyze Optimal Reserve -------------------------------------------------

OptRun <- filter(ReserveResults,BestRun == T)

opt_npb_plot <- opt_npb_plot_fun(OptRun = OptRun, Theme = SimpleTheme)

ggsave(file=paste(BatchFolder,'Optimal NPB Trajectory.pdf',sep=''),plot=opt_npb_plot,width=8,height=6)

opt_biomass_plot <- opt_biomass_plot_fun(OptRun = OptRun, Theme = SimpleTheme)

ggsave(file=paste(BatchFolder,'Optimal Biomass Trajectory.pdf',sep = ''),plot = opt_biomass_plot, width = fig_width, height = fig_height)

opt_reserve_plot <- opt_reserve_plot_fun(OptRun = OptRun, Theme = SimpleTheme)

ggsave(file=paste(BatchFolder,'Optimal Reserve Trajectory.pdf',sep=''),plot=opt_reserve_plot,width=fig_width,height=fig_height)


unified_npb_plot<- opt_npb_plot_fun(OptRun = filter(ReserveResults,BestUnifiedRun == T) , Theme = SimpleTheme )

ggsave(file=paste(BatchFolder,'Unified NPB Trajectory.pdf',sep=''),plot=unified_npb_plot,width=fig_width,height=fig_height)


unified_reserve_plot<- opt_reserve_plot_fun( OptRun = filter(ReserveResults,BestUnifiedRun == T) , Theme = SimpleTheme )

ggsave(file=paste(BatchFolder,'Unified Reserve Trajectory.pdf',sep=''),plot=unified_reserve_plot,width=fig_width,height=fig_height)

species_comp_plot <- species_comp_plot_fun(PlotData = species_comparison, Theme = theme_classic())

ggsave(file = paste(BatchFolder,'Opt Species Comparison.pdf',sep = ''),plot = species_comp_plot,width = fig_width,height = fig_height)


# Different Discount Rates ------------------------------------------------

Discounts <- c(0,0.05,0.2)

Opt_by_Discount <- lapply(Discounts, Best_Run_By_Discount, Runs = subset(ReserveResults, Intercept == 0 & Slope ==0)
                          , npb_focus = 1) %>% ldply()

discount_rate_plot <- discount_npb_plot_fun(filter(Opt_by_Discount, Species == 'Yellowtail Snapper'), SimpleTheme)

ggsave(file = paste(BatchFolder,'Discount Rate NPB.pdf',sep = ''),plot = discount_rate_plot,width = 5,height = 5)

# Analyze Static Reserve -----------------------------------------

StaticReserve <- filter(ReserveResults, Year == max(Year) & Intercept == 0 & Slope == 0)

StaticSummary <- filter(ResSummary, Intercept == 0 & Slope == 0)

static_NPB_plot <- static_NPB_plot_fun(PlotData = StaticSummary, Theme = SimpleTheme)

ggsave(file = paste(BatchFolder,'Static NPB.pdf',sep = ''),plot = static_NPB_plot,width = fig_width,height = 6)

static_NB_plot <- static_NB_plot_fun(PlotData = StaticSummary, Theme = SimpleTheme)

ggsave(file = paste(BatchFolder,'Static NB.pdf',sep = ''),plot = static_NB_plot,width = fig_width,height = fig_height)

static_NB2_plot <- static_NB2_plot_fun(PlotData = StaticSummary, Theme = SimpleTheme)

ggsave(file = paste(BatchFolder,'Static NB2.pdf',sep = ''),plot = static_NB2_plot,width = fig_width,height = fig_height)


static_comp_plot <- static_comp_plot_fun(PlotData = StaticSummary, Theme = SimpleTheme, DiscRate = Fleet$YieldDiscount)

ggsave(file = paste(BatchFolder,'Static Comp Plot.pdf',sep = ''),plot = static_comp_plot,width = fig_width,height = fig_height)

static_nb_v_npb_plot <- arrangeGrob(static_NPB_plot, static_NB_plot + theme(legend.position = 'none'), nrow = 1, ncol =2)

static_netbenefit_plot <- static_netbenefit_plot_fun(PlotData = StaticSummary, Theme = SimpleTheme)

ggsave(file = paste(BatchFolder,'Static Years to Net Benefit.pdf',sep = ''),plot = static_netbenefit_plot,width = fig_width,height = fig_height)

static__priceinc_plot <- static__priceinc_plot_fun(PlotData = StaticSummary, Theme = SimpleTheme)

ggsave(file = paste(BatchFolder,'Static Price Increase Needed.pdf',sep = ''),plot = static__priceinc_plot,width = fig_width,height = fig_height)

StaticSummary$LoanType[StaticSummary$MaxInterestRate <= 1] <- 'Philanthropy' 

StaticSummary$LoanType[StaticSummary$MaxInterestRate > 1 & StaticSummary$MaxInterestRate <=5] <- 'Social Investment' 

StaticSummary$LoanType[StaticSummary$MaxInterestRate > 5] <- 'Commercial' 

static__maxinterest_plot <- static__maxinterest_plot_fun(PlotData = StaticSummary, Theme = SimpleTheme)

ggsave(file = paste(BatchFolder,'Static Interest Rate Possible.pdf',sep = ''),plot = static__maxinterest_plot,width = fig_width,height = fig_height)

# Analyze Dynamic Reserves ------------------------------------------------

surface_results <- subset(ReserveResults,Year == max(Year) & Slope == fixed_slope)

surface_summary <- subset(ResSummary, Slope == fixed_slope)   

yield_surface_plot <- yield_surface_plot_fun(PlotData = surface_results, Theme = SimpleTheme)

ggsave(file=paste(BatchFolder,'Cumulative Yield Surface.pdf',sep=''),plot = yield_surface_plot, height = fig_height, width = fig_width)

npb_surface_plot <- npb_surface_plot_fun(PlotData = surface_results, Theme = SimpleTheme)

ggsave(file=paste(BatchFolder,'NPB Surface.pdf',sep=''),plot = npb_surface_plot, height = fig_height, width = fig_width)

priceinc_surface_plot <- priceinc_surface_plot_fun(PlotData = surface_summary, Theme = SimpleTheme)

ggsave(file=paste(BatchFolder,'Prince Inc Surface.pdf',sep=''),plot = priceinc_surface_plot, height = fig_height, width = fig_width)

surface_summary$LoanType[surface_summary$MaxInterestRate <= 5] <- 'Philanthropy' 

surface_summary$LoanType[surface_summary$MaxInterestRate > 5 & surface_summary$MaxInterestRate <=15] <- 'Social Investment' 

surface_summary$LoanType[surface_summary$MaxInterestRate > 15] <- 'Commercial' 

loan_surface_plot <- loan_surface_plot_fun(PlotData = surface_summary, Theme = SimpleTheme)

ggsave(file=paste(BatchFolder,'Loan Surface.pdf',sep=''),plot = loan_surface_plot, height = fig_height, width = fig_width)

loantype_surface_plot <- loantype_surface_plot_fun(PlotData = surface_summary, Theme = SimpleTheme)

ggsave(file=paste(BatchFolder,'Loan Type Surface.pdf',sep=''),plot = loantype_surface_plot, height = fig_height, width = fig_width)


NPB_model<- RobustRegression(lm(NPB ~ FinalReserve+Intercept+Slope+as.factor(Species),data=filter(ReserveResults,Year==max(Year))),
                             Data=filter(ReserveResults,Year==max(Year)),ClusterVar = 'None')

# Analyze Multispecies Scenario -------------------------------------------


multispecies <- ReserveResults %>%
  ungroup() %>%
  select(-Run) %>%
  rename(Run = Scenario) %>%
  group_by(Run, Year) %>%
  summarize( 
    Species = 'Multispecies',        
    SQYield = sum(SQYield), 
    Yield = sum(Yield), 
    CurrentReserve = mean(CurrentReserve), 
    FinalReserve = mean(FinalReserve), 
    Intercept = mean(Intercept), 
    Slope = mean(Slope))

save(file = paste(BatchFolder,'NutsPlots.Rdata',sep = ''),list=c('SimpleTheme',ls(pattern = '_plot')))



