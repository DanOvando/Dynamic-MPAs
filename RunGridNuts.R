
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
NumCores<- 4
InitialPopulation <- 1000 #Seed for initial population 
CollapseThreshold <- 0.1
LookAtLengths <- 0
ReservePosition <- 'Center'
OptTime <- 50 #Time Horizon to optimize over
Alpha <- 0.5

Font <- 'Helvetica'
FontColor <- 'black'
fig_width <- 8
fig_height <- 6
text_size <- 14
figure_use <- 'paper'

TimeToRun <- OptTime + 1
LoanTime <- 11
DiscRates <- seq(.05,.2,by = .05) #.1

BasePatches <- Patches

BatchFolder <- 'Results/Revised Version/'

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
  
  LifeHistories <- read.csv('Inputs/Life History Inputs.csv',stringsAsFactors = F) # Read in life history data
  
  LifeColumns <- colnames(LifeHistories)
  
  LifeVars <- c('Range','MaxAge','m','k','Linf','t0','AgeMa50','LengthMa50','MaturityMode','BH.Steepness','wa','wb','WeightForm','fa','fb','DDForm')
  
  SpeciesList <- LifeHistories$CommName
  
  SystemBmsyStorage <- as.data.frame(matrix(NA,nrow=length(SpeciesList),ncol=2))
  
  colnames(SystemBmsyStorage) <- c('Species','Bmsy')
  
  SystemBmsyStorage$Species <- as.character(SystemBmsyStorage$Species)
  
  NumFs <- 1
  
  DefaultLifeHistory <- lh #Base life history object
  
  BasePatches <- Patches #Base patch structure
  
  # Set up populations
  
  Populations<- sapply(gsub(' ','',SpeciesList,fixed = T),TargetPop=0.25,PreparePopulations,LifeHistories=LifeHistories,BaseLife=lh,LifeVars = LifeVars,BasePatches = BasePatches,BatchFolder = BatchFolder,USE.NAMES = T)

  RunMatrix <- PrepareGrid(SpeciesList,Fs='Set to 0.25',ReserveInc = 0.05,InterceptInc = 0.1,SlopeInc = 0.1,DiscRates)
  
  unique(RunMatrix$DiscountRate)

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
#     Rprof(tmp <- tempfile(),line.profiling=T)
    
    ReserveResults <- (mclapply(1:dim(RunMatrix)[1], RunGridReserve, RunMatrix = RunMatrix,BasePatches = BasePatches,
                                Populations = Populations, BatchFolder = BatchFolder,TimeToRun=TimeToRun, mc.cores = NumCores )) %>% ldply()
#     Rprof() 
#     summaryRprof(tmp)
#     unlink(tmp)
    
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




KeynoteTheme <- theme(legend.position = 'top', plot.background = element_rect(color = NA), rect=element_rect(fill = 'transparent', color = NA),
                      text = element_text(size = 22,family = Font, color= FontColor),
                      axis.text = element_text(color = FontColor),axis.title.y = element_text(size = 12,hjust= 0.5,angle = 0),
                      axis.text.y = element_text(size = 12),axis.text.x = element_text(angle = 35, vjust = 0.9,hjust = 0.9,color= FontColor,size = 12),
                      legend.text = element_text(size = 10,color= 'black'),legend.background=element_rect(fill = "gray90"),legend.title=element_blank())

PaperTheme <- theme(legend.position = 'top',text = element_text(size= 22,family= Font,color= FontColor),
                    axis.text = element_text(color= FontColor),axis.title.y=element_text(size = 25,hjust = 0.5,angle = 0)
                    ,axis.text.y = element_text( size = 30), axis.text.x = element_text(angle = 35, vjust = 0.9,hjust = 0.9,color = FontColor, size = 22),
                    legend.text = element_text(size = 14,color = 'black'),legend.title = element_text(size = 16,color = 'black'))

if (figure_use == 'paper')
{
  
  SimpleTheme <- theme(legend.position = 'top',text = element_text(size = text_size,family = Font,color = FontColor),
                       axis.text.x = element_text(size = 9, angle = 35, vjust = 0.9, hjust = 0.9),
                       panel.background =  element_rect(fill = "white", colour = NA),
                       panel.border =      element_rect(fill = NA, colour = "grey50"), 
                       panel.grid.major =  element_line(colour = "grey60", size = 0.2),
                       panel.grid.minor =  element_line(colour = "grey78", size = 0.5))
}                 
if (figure_use == 'presentation')
{
  SimpleTheme <- theme(legend.position = 'top',text = element_text(size = text_size,family = Font,color = FontColor),
                       axis.text.x = element_text(size = text_size, angle = 35, vjust = 0.9, hjust = 0.9,color = FontColor),
                       axis.text.y = element_text(color = FontColor),
                       plot.background=element_rect(color=NA),
                       rect = element_rect(fill='transparent',color=NA),
                       panel.background =  element_rect(fill = "white", colour = NA),
                       panel.border =      element_rect(fill = NA, colour = "grey50"), 
                       panel.grid.major =  element_line(colour = "grey60", size = 0.2),
                       panel.grid.minor =  element_line(colour = "grey78", size = 0.5),
                       strip.text = element_text(color = 'black'))
}
# Process Results ---------------------------------------------------------

Fleet$YieldDiscount <- discount_rate

ProcessedNuts <- ProcessNuts(ReserveResults = ReserveResults, Fleet = Fleet, Scale_Yields = Scale_Yields, npb_focus = npb_focus)

ReserveResults <- ProcessedNuts$ReserveResults

species_comparison <- ProcessedNuts$species_comparison

ResSummary <- ProcessedNuts$ResSummary


# Cartoon Plots -----------------------------------------------------------

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
  
  #long format
  longdat <- reshape(datA[,1:3], 
                     varying = c("Biomass","Yield"), 
                     v.names = "value",
                     timevar = "Objectives", 
                     times = c("Biomass","Yield"), 
                     direction = "long")
  
  
  
  
  p <- ggplot(longdat, aes(x=Years,y=value))
  
  p + geom_area(data=subset(longdat,Years<=12 & Objectives=="Yield"),  fill="darkorange")+ geom_area(data=subset(longdat,Years>=12 &Years<=breakYear& Objectives=="Yield"),  fill="cornflowerblue")+ geom_line(aes(colour=Objectives),size=1.1)+ geom_vline(xintercept = 30,linetype=2) + scale_color_manual(values=c("aquamarine3", "black"))  + ylab("Deviation from Status Quo")  + scale_x_discrete(breaks=c(0,30),labels=c("0","Yield Balance"))     + scale_y_continuous(limits=c(-0.1, 0.2),breaks=0) + Theme
  
  
#   p + geom_area(data=subset(longdat,Years<=12 & Objectives=="Yield"),  fill="grey26")+ geom_area(data=subset(longdat,Years>=12 &Years<=breakYear& Objectives=="Yield"),  fill="grey68")+ geom_line(aes(linetype=Objectives),size=1.1)+ geom_vline(xintercept = 30,linetype=2) + scale_color_manual(values=c("black", "black"))  + ylab("Deviation from Status Quo")  + scale_x_discrete(breaks=c(0,30),labels=c("0","Yield Balance"))     + scale_y_continuous(limits=c(-0.1, 0.2),breaks=0) + Theme
  
  
  
}

cartoon_plot <- CartoonPlot(SimpleTheme)

ggsave(file=paste(BatchFolder,'Cartoon Plots.pdf',sep=''),plot=cartoon_plot,width=fig_width,height=fig_height)


best_static_run <- filter(ReserveResults, Slope == 0 & Intercept == 0 & Year == max(Year)) %>%
  group_by(Species) %>%
  summarize(best_npb_run = unique(Run[NPB == max(NPB,na.rm = T)]),
            best_yield_run = unique(Run[Yield == max(Yield,na.rm = T)]))
  

ReserveResults <- ReserveResults %>%
  group_by(Run) %>%
  mutate(BiomassBalance = Biomass-Biomass[1])


npb_illustration_plot <- (ggplot(subset(ReserveResults, Species == 'Yellowtail Snapper' 
              & Run == best_static_run$best_npb_run[best_static_run$Species == 'Yellowtail Snapper'] & Year < 20),
       aes(factor(Year-1),YieldBalance/max(YieldBalance), fill = NPB/max(NPB))) + geom_bar(color = 'black',stat = 'identity',position = 'dodge') + 
  scale_fill_gradient2(name = 'NPB',low = 'red',mid = 'white',high = 'grey0', midpoint = 0) + 
  SimpleTheme + 
  geom_point(aes(x = factor(Year -1), y = BiomassBalance/max(BiomassBalance), size = Biomass/max(Biomass)),fill = 'green',shape = 21)
  + scale_size_continuous(name = 'Biomass') + 
  xlab('Year') + 
  ylab('Change From Status Quo') + 
  geom_hline(aes(yintercept = 0)))

ggsave(file=paste(BatchFolder,'NPB Illustration.pdf',sep=''),plot=npb_illustration_plot,width=fig_width,height=fig_height)

# Analyze Optimal Reserve -------------------------------------------------

OptRun <- filter(ReserveResults,BestRun == T)

# One option. Add in the metrics for EQ reserve of size finalreserve
# Another option Add in metrics for best EQ reserve by species


opt_npb_plot <- opt_npb_plot_fun(OptRun = OptRun, Theme = SimpleTheme)

ggsave(file=paste(BatchFolder,'Optimal NPB Trajectory.pdf',sep=''),plot=opt_npb_plot,width=fig_width,height=fig_height)

opt_biomass_plot <- opt_biomass_plot_fun(OptRun = OptRun, Theme = SimpleTheme)

ggsave(file=paste(BatchFolder,'Optimal Biomass Trajectory.pdf',sep = ''),plot = opt_biomass_plot, width = fig_width, height = fig_height)

browser()
opt_reserve_plot <- opt_reserve_plot_fun(OptRun = OptRun, Theme = SimpleTheme)

ggsave(file=paste(BatchFolder,'Optimal Reserve Trajectory.pdf',sep=''),plot=opt_reserve_plot,width=fig_width,height=fig_height)


unified_npb_plot<- opt_npb_plot_fun(OptRun = filter(ReserveResults,BestUnifiedRun == T) , Theme = SimpleTheme )

ggsave(file=paste(BatchFolder,'Unified NPB Trajectory.pdf',sep=''),plot=unified_npb_plot,width=fig_width,height=fig_height)


unified_reserve_plot<- opt_reserve_plot_fun( OptRun = filter(ReserveResults,BestUnifiedRun == T) , Theme = SimpleTheme )

ggsave(file=paste(BatchFolder,'Unified Reserve Trajectory.pdf',sep=''),plot=unified_reserve_plot,width=fig_width,height=fig_height)

species_comp_plot <- species_comp_plot_fun(PlotData = species_comparison, Theme = SimpleTheme)

# species_comp_plot <- species_comp_plot_fun(PlotData = species_comparison, Theme = theme_classic())

ggsave(file = paste(BatchFolder,'Opt Species Comparison.pdf',sep = ''),plot = species_comp_plot,width = fig_width,height = fig_height)


# Different Discount Rates ------------------------------------------------

Discounts <- c(0,0.2)

# Discounts <- c(0,0.2)


Opt_by_Discount <- lapply(Discounts, Best_Run_By_Discount, Runs = subset(ReserveResults, Intercept == 0 & Slope ==0)
                          , npb_focus = 1) %>% ldply()

discount_rate_plot <- discount_npb_plot_fun(filter(Opt_by_Discount, Species == 'Yellowtail Snapper'), SimpleTheme)

ggsave(file = paste(BatchFolder,'Discount Rate NPB.pdf',sep = ''),plot = discount_rate_plot,width = fig_width,height = fig_height)

ConsRun <- unique(Opt_by_Discount$FinalReserve[Opt_by_Discount$DiscountRate == 0 & Opt_by_Discount$Species == 'Yellowtail Snapper'])

Discounts <- c(0,0.2)

target_by_Discount <- lapply(Discounts, Target_Run_By_Discount, Runs = subset(ReserveResults, Intercept == 0 & Slope ==0)
                          , npb_focus = 1,target_reserve = ConsRun) %>% ldply()

conserve_discount_rate_plot <- discount_npb_plot_fun(filter(target_by_Discount, Species == 'Red Grouper'), SimpleTheme)

ggsave(file = paste(BatchFolder,'Conserve Discount Rate NPB.pdf',sep = ''),plot = conserve_discount_rate_plot,width = fig_width,height = fig_height)




# Analyze Static Reserve -----------------------------------------

StaticReserve <- filter(ReserveResults, Year == max(Year) & Intercept == 0 & Slope == 0)

StaticSummary <- filter(ResSummary, Intercept == 0 & Slope == 0)

static_NPB_plot <- static_NPB_plot_fun(PlotData = StaticSummary, Theme = SimpleTheme)

ggsave(file = paste(BatchFolder,'Static NPB.pdf',sep = ''),plot = static_NPB_plot,width = fig_width,height = fig_height)

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



