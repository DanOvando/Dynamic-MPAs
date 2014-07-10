##-----------------------------------------
## Regression trees: deterministic results
## CM, AC: Thu Jun 13 13:43:20 2013
## Notes:
## Awaiting final results
## Accounts for grouping and time series autocorrelation 
##-----------------------------------------
rm(list=ls())
setwd('/Users/danovando/Dropbox/Shrinking NTZ')
library(REEMtree)
library(rpart.plot) #need different version

## note: not final (old)
modeldata <- read.csv("Deterministic_performance_results.csv")
head(modeldata)

table(modeldata$method_id)

## attempted to code to ease transition to
## a function that will use various time ranges
## and summary statistics

## take proportional error in last 5 years as statistic first

stat.name<-"prop.error"
year.sub<-56:60
method.name<-"Costello"
## subset the data
sub.modeldata<-subset(modeldata, year%in%year.sub & method_id==method.name)

my.formula<-formula(bquote(.(as.name(stat.name))~LH+ID+ED+TS))

## ignoring grouping and autocorrelation
rtree.nocorr <- rpart(my.formula, data=sub.modeldata,control=list(minsplit=100))
rtree.nocorr <- rpart(my.formula, data=sub.modeldata,control=list(maxdepth=2))
rtree.nocorr <- rpart(my.formula, data=sub.modeldata,control=list(cp=.1))

plot(rtree.nocorr) ## note some fancy plots below
print(rtree.nocorr)
plot(predict(rtree.nocorr),residuals(rtree.nocorr))

## incorporating random effects
rtree.re <- REEMtree(my.formula, 
                     random = ~1|stock_id, 
                     data=sub.modeldata,
                     verbose=TRUE,
                     #tree.control=list(minsplit=100),
                     #tree.control=list(minbucket=100),
                     #tree.control=list(mincut=100),
                     #tree.control=list(cp=.1),
                     cpmin=.04,
                     method="REML")
plot(rtree.re)
print(rtree.re)
plot(fitted(rtree.re),residuals(rtree.re))
## note variance increase

## random effect and ar1 on resids
rtree.re.ar1 <- REEMtree(my.formula, 
                      random = ~1|stock_id, 
                      correlation=corAR1(,~year),
                      data=sub.modeldata,
                      verbose=TRUE,
                      method="ML")

plot(rtree.re.ar1)
print(rtree.re.ar1)
plot(fitted(rtree.re.ar1),residuals(rtree.re.ar1))


AIC(rtree.re.ar1,rtree.re) ## ballpark - note difference in length of y

## some fancy plots

## color palette
blue2red<-colorRampPalette(c("blue","white","red"))
cols<-blue2red(100)
max.stat<-max(abs(modeldata[,stat.name]), na.rm=TRUE)
## breaks will need to be worked on for others
breaks<-seq(-1.1*max.stat, 1.1*max.stat, length=101)


pdf("regression_tree_example.pdf", height=8, width=10)
col.index<-cut(rtree.nocorr$frame$yval,breaks)
##
## no corr
par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(2,2,1,1))
prp(rtree.nocorr, box.col=cols[col.index],faclen=0, do.par=FALSE)
legend("topleft", legend="Independent", bty="n")
box()
## random effect by stock id
col.index<-cut(rtree.re$Tree$frame$yval,breaks)
prp(rtree.re$Tree, box.col=cols[col.index],faclen=0, do.par=FALSE)
legend("topleft", legend="RE by stock", bty="n")
box()
##
## random effect by stock id and AR(1) on resids
col.index<-cut(rtree.re.ar1$Tree$frame$yval,breaks)
prp(rtree.re.ar1$Tree, box.col=cols[col.index],faclen=0, do.par=FALSE)
legend("topleft", legend="RE and AR(1) by stock", bty="n")
box()
dev.off()

