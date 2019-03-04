setwd("/Users/royr/UCSF/WallaceMarshall")

####################################################################
####################################################################

load("data.RData")

medW=apply(cellW2,2,median,na.rm=T)
medM=apply(cellM2,2,median,na.rm=T)
sdW=apply(cellW2,2,sd,na.rm=T)
sdM=apply(cellM2,2,sd,na.rm=T)
medWM=apply(rbind(cellW2,cellM2),2,median,na.rm=T)
sdWM=apply(rbind(cellW2,cellM2),2,sd,na.rm=T)

summary(medW)
summary(medM)
summary(medWM)

summary(sdW)
summary(sdM)
summary(sdWM)

"
> summary(medW)
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
-936.00      0.15      0.93   5882.31     66.00 213000.00
> summary(medM)
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
-1031.0      0.1      1.0  17260.1     57.0 868880.0
> summary(medWM)
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
-967.00      0.15      0.93   6876.40     66.23 268000.00
>
> summary(sdW)
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
0.0      0.1      0.4  12783.6     71.4 374417.1
> summary(sdM)
Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
0.0       0.0       0.4   29702.8      55.7 1036531.5
> summary(sdWM)
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
0.0      0.1      0.4  21220.4     71.3 659977.3
"

png("densityPlot_median.png",width=480,height=3*480)
par(mfcol=c(3,1))
xlim=c(-1031,268000); ylim=c(0,1)
plot(density(medW),na.rm=T)
plot(density(medM),na.rm=T)
plot(density(medWM),na.rm=T)
dev.off()

png("densityPlot_sd.png",width=480,height=3*480)
par(mfcol=c(3,1))
xlim=c(0,1036531.5); ylim=c(0,1)
plot(density(sdW),na.rm=T)
plot(density(sdM),na.rm=T)
plot(density(sdWM),na.rm=T)
dev.off()


png("plot_med_sd.png",width=480,height=480)
par(mfcol=c(2,2))
lim=range(c(medW,medM),na.rm=T)
plot(medW,medM,xlim=lim,ylim=lim); abline(c(0,1),col="red")
lim=range(c(sdW,sdM),na.rm=T)
plot(sdW,sdM,xlim=lim,ylim=lim); abline(c(0,1),col="red")
dev.off()

####################################################################
####################################################################
## Split wt into 2 groups, train 1. See how well group 2 fall & compare the myc/ras
##

set.seed(53453)
n=nrow(cellM2)
j=sample(1:nrow(cellW2),2*n,replace=F)
j1=j[1:n]; j2=j[(n+1):2*n]



####################################################################
####################################################################
