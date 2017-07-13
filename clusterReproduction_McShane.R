## Use clustRep.R to get cluster designations of perturbed data

computerFlag="cluster"
computerFlag=""

## ----------------------------------------------
if (computerFlag=="cluster") {
    setwd("/home/royr/project/WallaceMarshall")
} else {
    dirSrc="/Users/royr/UCSF/"
    dirSrc2=dirSrc
    setwd(paste(dirSrc2,"WallaceMarshall",sep=""))
}

## ----------------------------------------------

####################################################################
####################################################################

repFlag="_fold"
repFlag="_perm"

cohort="_mycRasWt"; typeFlag=""; type2Flag="_reducedBioFeatPC"; subsetFlag=""
distMethod="spearman"
distMethod="pearson"
distMethod="kendall"
linkMethod="ward.D2"

centrFlag=""
scaleFlag=""

orderFlag="_wtOrd"
featureFlag="sample"

stdDev="_sd0.2"
stdDev="_sd0.05"
stdDev="_sd0.1"

switch(cohort,
    "_mycRas"={
        cohortName="Myc/Ras"
    },
    "_wt"={
        cohortName="Wildtype"
    },
    "_mycRasWt"={
        cohortName="Myc/Ras & Wildtype"
    }
)

load(file=paste("arrayData2",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,".RData",sep=""))

datadir2=paste("results/heatmap/",sub("_","",sub("Ord","",orderFlag)),"/",sub("_","",orderFlag),"/",sub("_","",type2Flag),"/",sep="")
datadir2=paste("results/w3830/heatmap/",sub("_","",cohort),"/",sub("_","",orderFlag),"/",sub("_","",type2Flag),"/",sep="")
datadir2=""
fNameOut4=paste(cohort,orderFlag,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,centrFlag,scaleFlag,"_",distMethod,sep="")
clustInfo=read.table(paste(datadir2,fNameOut4,"/clusterInfoSample",fNameOut4,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
clustInfo$clustId=as.integer(sub("cluster","",clustInfo$clustId))

heading=paste(cohortName,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),subsetFlag,centrFlag,scaleFlag,", ",distMethod,ifelse(repFlag=="_perm",paste(", sd ",sub("_sd","",stdDev),sep=""),", 2fold"),sep="")

datadir3="clusterReproduction/"
fNameOut=paste(cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,centrFlag,scaleFlag,"_",distMethod,sep="")
if (repFlag=="_perm") {
    load(paste(datadir3,"clustObs",fNameOut,".RData",sep=""))
    fileList=dir(datadir3)
    fileList=fileList[grep("clustPerm_",fileList)]
    fileList=fileList[grep(paste(fNameOut,stdDev,sep=""),fileList)]
    clustPerm1=NULL
    for (fId in 1:length(fileList)) {
        load(paste(datadir3,fileList[fId],sep=""))
        clustPerm1=cbind(clustPerm1,clustPerm)
    }
    clustPerm=clustPerm1
    rm(clustPerm1)
} else {
    load(paste(datadir3,"clustFold",fNameOut,"_2fold.RData",sep=""))
    clustObs=clustFold[,1]
    clustPerm=clustFold[,c(2,2)]
}

clustObs=clustObs[match(clustInfo$id,names(clustObs))]
clustPerm=clustPerm[match(clustInfo$id,rownames(clustPerm)),]

## Check the number of samples in a cluster
## Eliminate singletons
xxx=rep(0,14)
for (nClust in 2:15) {
    cat("\n\n========= nClust ",nClust," ==========\n",sep="")
    x0=sapply(clustObs,function(x) as.integer(strsplit(x,",")[[1]][nClust-1]),USE.NAMES=F)
    xx=sapply(1:nClust,function(iClust) {
        j0=which(x0==iClust)
        length(j0)
    },USE.NAMES=F)
    xxx[nClust-1]=min(xx)
}
## Clusters with least number of samples
summary(xxx)

## ----------------------------------------------
fNameOut4=paste(fNameOut,ifelse(repFlag=="_perm",stdDev,"_2fold"),sep="")
if (F) {
    ## Calculate r-index
    rIndexMat=matrix(nrow=14,ncol=ncol(clustPerm))
    nClust=2
    timeStamp=Sys.time()
    print(format(timeStamp, "%x %X"))
    for (nClust in 2:15) {
        cat("\n\n========= nClust ",nClust," ==========\n",sep="")
        x0=sapply(clustObs,function(x) as.integer(strsplit(x,",")[[1]][nClust-1]),USE.NAMES=F)
        rIndexVec=rep(NA,ncol(clustPerm))
        for (rId in 1:length(rIndexVec)) {
            x=sapply(clustPerm[,rId],function(x) as.integer(strsplit(x,",")[[1]][nClust-1]),USE.NAMES=F)
            y=sapply(1:nClust,function(iClust) {
                j0=which(x0==iClust)
                zzz=0
                for (j1 in 1:(length(j0)-1)) {
                    j=(j1+1):length(j0)
                    zzz=zzz+sum(x[j0[j]]==x[j0[j1]])
                }
                c(zzz,length(j0)*(length(j0)-1)/2)
            },USE.NAMES=F)
            zzz=apply(y,1,sum)
            rIndexVec[rId]=zzz[1]/zzz[2]
        }
        rIndexMat[nClust-1,]=rIndexVec
    }
    timeStamp=c(timeStamp,Sys.time())
    print(format(timeStamp[2], "%x %X"))
    print(diff(timeStamp))
    save(rIndexMat,file=paste("rIndexMat",fNameOut4,".RData",sep=""))
}

datadir3="clusterReproduction/"
load(file=paste(datadir3,"rIndexMat",fNameOut4,".RData",sep=""))

rIndex=apply(rIndexMat,1,mean)

rIndex
order(rIndex,decreasing=T)+1

"

## w5969
## Kendall, perm sd=.2
> rIndex
[1] 0.7166948 0.5210458 0.4658965 0.4875449 0.4696809 0.4565775 0.4534397 0.4371416 0.4264433 0.4230367
[11] 0.4096585 0.3980088 0.3906103 0.3955113
> order(rIndex,decreasing=T)+1
[1]  2  3  5  6  4  7  8  9 10 11 12 13 15 14

## Kendall, perm sd=.1
> rIndex
[1] 0.7155110 0.5307929 0.4868207 0.5046960 0.4771619 0.4701034 0.4673336 0.4505908 0.4480373
[10] 0.4447264 0.4321384 0.4189885 0.4078637 0.4104430
> order(rIndex,decreasing=T)+1
[1]  2  3  5  4  6  7  8  9 10 11 12 13 15 14

## Kendall, perm sd=.05
> rIndex
[1] 0.7088255 0.5307728 0.4960546 0.5069154 0.4836639 0.4718761 0.4749769 0.4584177 0.4531488
[10] 0.4468273 0.4329869 0.4231333 0.4151060 0.4198464
> order(rIndex,decreasing=T)+1
[1]  2  3  5  4  6  8  7  9 10 11 12 13 15 14


## Pearson, perm sd=.2
> rIndex
[1] 0.7826496 0.6250258 0.5745302 0.5203695 0.4999606 0.4655163 0.4636076 0.4603992 0.4507791
[10] 0.4722170 0.4764931 0.4770515 0.4692578 0.4616856
> order(rIndex,decreasing=T)+1
[1]  2  3  4  5  6 13 12 11 14  7  8 15  9 10

## Pearson, perm sd=.1
> rIndex
[1] 0.7857816 0.6241487 0.5901042 0.5330351 0.5112914 0.4778175 0.4767687 0.4752938 0.4653046
[10] 0.4893916 0.4920757 0.4966562 0.4884966 0.4801348
> order(rIndex,decreasing=T)+1
[1]  2  3  4  5  6 13 12 11 14 15  7  8  9 10

## Pearson, perm sd=.05
> rIndex
[1] 0.8033804 0.6278978 0.5965995 0.5445826 0.5204504 0.4869816 0.4866617 0.4859698 0.4758320
[10] 0.4977957 0.5037811 0.5087227 0.5015185 0.4944858
> order(rIndex,decreasing=T)+1
[1]  2  3  4  5  6 13 12 14 11 15  7  8  9 10

## Kendall, 2fold
> rIndex
[1] 0.8125293 0.8510211 0.5610775 0.5765489 0.5827885 0.4980182 0.3169000 0.3206685 0.3203645 0.3242540
[11] 0.3405379 0.2686577 0.3873180 0.3903995
order(rIndex,decreasing=T)+1
[1]  3  2  6  5  4  7 15 14 12 11  9 10  8 13


## w3830
## Kendall, perm sd=.2
> rIndex
[1] 0.7150903 0.5459151 0.5464257 0.5035394 0.4744973 0.4719460 0.4555063 0.4480241 0.4350241
[10] 0.4224851 0.4195369 0.4121903 0.4120911 0.4110929
order(rIndex,decreasing=T)+1
[1]  2  4  3  5  6  7  8  9 10 11 12 13 14 15

## Spearman, perm sd=.1
> rIndex
[1] 0.6671929 0.6001567 0.5183934 0.4958563 0.4810826 0.4663948 0.4573841 0.4497517 0.4573107
[10] 0.4614046 0.4523854 0.4474784 0.4448312 0.4481154
# Best no. of clusters in decreasing order
order(rIndex,decreasing=T)+1
# 2  3  4  5  6  7 11  8 10 12  9 15 13 14

## Spearman, perm sd=.2
> rIndex
[1] 0.6735515 0.5934907 0.5294806 0.5099775 0.4896137 0.4804254 0.4753499 0.4715280 0.4816780
[10] 0.4860125 0.4776211 0.4734514 0.4713723 0.4751431
# Best no. of clusters in decreasing order
order(rIndex,decreasing=T)+1
[1]  2  3  4  5  6 11 10  7 12  8 15 13  9 14

"

## ----------------------------------------------
fNameOut4=paste(fNameOut,ifelse(repFlag=="_perm",stdDev,"_2fold"),sep="")
if (F) {
    ## Calculate d-index
    dIndexMatO=dIndexMatB=matrix(nrow=14,ncol=ncol(clustPerm))
    nClust=2
    timeStamp=Sys.time()
    print(format(timeStamp, "%x %X"))
    for (nClust in 2:15) {
        cat("\n\n========= nClust ",nClust," ==========\n",sep="")
        x0=sapply(clustObs,function(x) as.integer(strsplit(x,",")[[1]][nClust-1]),USE.NAMES=F)
        nm0=names(x)
        for (dId in 1:ncol(dIndexMatO)) {
            x=sapply(clustPerm[,dId],function(x) as.integer(strsplit(x,",")[[1]][nClust-1]),USE.NAMES=F)
            nm=names(x)
            clust0=table(x)
            y=sapply(1:nClust,function(iClust) {
                j0=which(x0==iClust)
                clust=table(x[j0])
                clustId=as.integer(names(clust)[which(clust==max(clust))])
                if (length(clustId)!=1) {
                    clustId=clustId[which.min(clust0[clustId])]
                }
                j=which(x==clustId)
                #zO=sum(!names(x0)[j0]%in%names(x)[j])
                #zB=sum(!names(x)[j]%in%names(x0)[j0])
                zO=sum(!nm0[j0]%in%nm[j])
                zB=sum(!nm[j]%in%nm0[j0])
                c(zO,zB)
            },USE.NAMES=F)
            zzz=apply(y,1,sum)
            dIndexMatO[nClust-1,dId]=zzz[1]
            dIndexMatB[nClust-1,dId]=zzz[2]
        }
    }
    timeStamp=c(timeStamp,Sys.time())
    print(format(timeStamp[2], "%x %X"))
    print(diff(timeStamp))
    save(dIndexMatO,dIndexMatB,file=paste("dIndexMat",fNameOut4,".RData",sep=""))
}

datadir3="clusterReproduction/"
load(file=paste(datadir3,"dIndexMat",fNameOut4,".RData",sep=""))

dIndexO=apply(dIndexMatO,1,mean)
dIndexB=apply(dIndexMatB,1,mean)

dIndexO
order(dIndexO,decreasing=F)+1
dIndexB
order(dIndexB,decreasing=F)+1

"

## w5969

## Kendall, perm sd=.2
> dIndexO
[1]  969.668 1802.020 1965.188 1814.288 1936.804 2022.936 2118.436 2210.268 2248.020 2233.472 2282.344
[12] 2328.680 2310.948 2271.092
> order(dIndexO,decreasing=F)+1
[1]  2  3  5  6  4  7  8  9 11 10 15 12 14 13
> dIndexB
[1] 1216.684 2364.252 2432.164 2075.864 2335.876 2409.432 2657.708 2983.532 2905.988 2805.216 2777.944
[12] 2810.692 2773.812 2799.180
> order(dIndexB,decreasing=F)+1
[1]  2  5  6  3  7  4  8 14 12 15 11 13 10  9

## Kendall, perm sd=.1
> dIndexO
[1]  934.6636 1698.5818 1867.3909 1770.8273 1920.4727 1989.2273 2075.8545 2159.1364 2169.9818
[10] 2153.8455 2203.6091 2242.6091 2232.5636 2198.1727
> order(dIndexO,decreasing=F)+1
[1]  2  3  5  4  6  7  8 11  9 10 15 12 14 13
> dIndexB
[1] 1010.845 2243.955 2553.764 2006.345 2358.309 2516.064 2703.855 2897.709 2818.745 2696.000
[11] 2689.018 2664.755 2665.791 2714.082
> order(dIndexB,decreasing=F)+1
[1]  2  5  3  6  7  4 13 14 12 11  8 15 10  9


## Kendall, perm sd=.05
> dIndexO
[1]  968.08 1699.41 1817.50 1762.46 1910.31 1972.98 2025.89 2120.00 2133.81 2131.91 2181.07 2204.13
[13] 2184.23 2130.75
> order(dIndexO,decreasing=F)+1
[1]  2  3  5  4  6  7  8  9 15 11 10 12 14 13
> dIndexB
[1] 1055.66 2240.06 2632.51 2036.89 2341.00 2381.21 2570.33 2917.34 2833.38 2742.24 2660.20 2630.19
[13] 2639.43 2686.92
> order(dIndexB,decreasing=F)+1
[1]  2  5  3  6  7  8 13  4 14 12 15 11 10  9


## Pearson, perm sd=.2
> dIndexO
[1] 1066.031 1623.733 1886.755 2373.496 2359.796 2560.992 2691.168 2785.890 2847.789 2858.411
[11] 2856.029 2923.888 2907.549 2923.622
> order(dIndexO,decreasing=F)+1
[1]  2  3  4  6  5  7  8  9 10 12 11 14 15 13
> dIndexB
[1] 1774.493 1776.427 1975.724 3495.304 2811.244 3145.576 3225.603 3407.159 3481.153 3855.499
[11] 3938.848 4180.807 3825.052 3685.078
> order(dIndexB,decreasing=F)+1
[1]  2  3  4  6  7  8  9 10  5 15 14 11 12 13

## Pearson, perm sd=.1
> dIndexO
[1] 1030.115 1613.234 1803.293 2303.186 2294.451 2485.712 2614.018 2699.888 2757.780 2758.431
[11] 2758.706 2801.563 2790.311 2805.890
> order(dIndexO,decreasing=F)+1
[1]  2  3  4  6  5  7  8  9 10 11 12 14 13 15
> dIndexB
[1] 1560.493 1844.164 1888.759 3379.916 2704.184 3047.583 3120.679 3287.533 3350.277 3718.804
[11] 3777.326 4024.617 3676.988 3534.335
> order(dIndexB,decreasing=F)+1
[1]  2  3  4  6  7  8  9 10  5 15 14 11 12 13

## Pearson, perm sd=.05
> dIndexO
[1]  948.022 1590.848 1762.640 2239.158 2230.704 2419.966 2545.269 2627.542 2687.908 2703.640
[11] 2686.167 2723.764 2712.726 2724.184
> order(dIndexO,decreasing=F)+1
[1]  2  3  4  6  5  7  8  9 12 10 11 14 13 15
> dIndexB
[1] 1369.616 1826.549 1837.570 3226.082 2628.997 2956.067 3016.905 3261.491 3289.299 3611.507
[11] 3683.220 3884.403 3561.610 3432.480
> order(dIndexB,decreasing=F)+1
[1]  2  3  4  6  7  8  5  9 10 15 14 11 12 13


## Kendall, 2fold
> dIndexO
[1]  955  596 1600 1508 1562 1761 2732 2935 2802 2732 2814 3125 2684 2653
> order(dIndexO,decreasing=F)+1
[1]  3  2  5  6  4  7 15 14  8 11 10 12  9 13
> dIndexB
[1] 3922  596 2527 3514 3216 2900 3237 3902 3126 3686 4712 4916 5275 4225
> order(dIndexB,decreasing=F)+1
[1]  3  4  7 10  6  8  5 11  9  2 15 12 13 14



## w3830
## Kendall, perm sd=.2
> dIndexO
[1] 2341.752 1674.892 1624.060 1771.784 1900.800 1955.300 2051.636 2117.508 2124.164 2177.752 2238.372
[12] 2245.704 2285.776 2322.980
> order(dIndexO,decreasing=F)+1
[1]  4  3  5  6  7  8  9 10 11 12 13 14 15  2
> dIndexB
[1] 3135.736 1872.972 1755.328 2363.488 2470.164 2478.408 2858.332 2909.572 2791.720 2748.416 2770.656
[12] 2739.352 2854.568 2953.404
> order(dIndexB,decreasing=F)+1
[1]  4  3  5  6  7 13 11 12 10 14  8  9 15  2

"

png(paste("rIndex_dIndex_mcshane",fNameOut4,".png",sep=""))
par(mfcol=c(2,2))
lim=range(c(dIndexO,dIndexB))
plot((1:length(dIndexO))+1, dIndexO, ylim=lim, type="b", main=paste(heading,sep=""), xlab="Number of clusters",ylab="D-index (mismatched observed)")
plot((1:length(dIndexB))+1, dIndexB, ylim=lim, type="b", main=paste(heading,sep=""), xlab="Number of clusters",ylab="D-index (mismatched best cluster)")
lim=c(0,1)
plot((1:length(rIndex))+1, rIndex, ylim=lim, type="b", main=paste(heading,sep=""), xlab="Number of clusters",ylab="R-index")
dev.off()

## --------------------------
clustFlag=""
clustFlag="_with2Cl"

limR=limDo=limDb=c()
for (distMethod in c("pearson","kendall")) {
    for (stdDev in c("_sd0.2","_sd0.1","_sd0.05")) {
        fNameOut4=paste(fNameOut,ifelse(repFlag=="_perm",stdDev,"_2fold"),sep="")
        datadir3="clusterReproduction/"
        load(file=paste(datadir3,"rIndexMat",fNameOut4,".RData",sep=""))
        load(file=paste(datadir3,"dIndexMat",fNameOut4,".RData",sep=""))
        rIndex=apply(rIndexMat,1,mean)
        dIndexO=apply(dIndexMatO,1,mean)
        dIndexB=apply(dIndexMatB,1,mean)
        if (clustFlag=="") {
            i=2:length(rIndex)
        } else {
            i=1:length(rIndex)
        }
        limR=c(limR,range(rIndex[i]))
        limDo=c(limDo,range(dIndexO[i]))
        limDb=c(limDb,range(dIndexB[i]))
    }
}
limR=range(limR)
limDo=range(limDo)
limDb=range(limDb)

for (distMethod in c("pearson","kendall")) {
    pdf(paste("rIndex_dIndex_mcshane_",distMethod,clustFlag,".pdf",sep=""),width=7*2, height=7)
    par(mfcol=c(3,3))
    for (stdDev in c("_sd0.2","_sd0.1","_sd0.05")) {
        heading=paste(cohortName,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),subsetFlag,centrFlag,scaleFlag,", ",distMethod,ifelse(repFlag=="_perm",paste(", sd ",sub("_sd","",stdDev),sep=""),", 2fold"),sep="")
        fNameOut4=paste(fNameOut,ifelse(repFlag=="_perm",stdDev,"_2fold"),sep="")
        datadir3="clusterReproduction/"
        load(file=paste(datadir3,"rIndexMat",fNameOut4,".RData",sep=""))
        load(file=paste(datadir3,"dIndexMat",fNameOut4,".RData",sep=""))
        rIndex=apply(rIndexMat,1,mean)
        rIndex
        order(rIndex,decreasing=T)+1
        dIndexO=apply(dIndexMatO,1,mean)
        dIndexB=apply(dIndexMatB,1,mean)
        dIndexO
        order(dIndexO,decreasing=F)+1
        dIndexB
        order(dIndexB,decreasing=F)+1
        if (clustFlag=="") {
            i=2:length(rIndex)
        } else {
            i=1:length(rIndex)
        }
        offset=i[1]
        lim=limDo
        plot((1:length(dIndexO[i]))+offset, dIndexO[i], ylim=lim, type="b", main=paste(heading,sep=""), xlab="Number of clusters",ylab="D-index (mismatched observed)")
        lim=limDb
        plot((1:length(dIndexB[i]))+offset, dIndexB[i], ylim=lim, type="b", main=paste(heading,sep=""), xlab="Number of clusters",ylab="D-index (mismatched best cluster)")
        lim=limR
        plot((1:length(rIndex[i]))+offset, rIndex[i], ylim=lim, type="b", main=paste(heading,sep=""), xlab="Number of clusters",ylab="R-index")
    }
    dev.off()
}


####################################################################
####################################################################

set.seed(5354)
reproduceCluster=function(dat,nClust,stdDev,nPerm) {

    permSetFlag="1"; set.seed(5354)
    
    cohort="_mycRasWt"; typeFlag=""; type2Flag="_reducedBioFeatPC"; subsetFlag=""
    load(file=paste("arrayData2",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,".RData",sep=""))

    dat=arrayData2[1:4,1:5]; nClust=2; stdDev=0.1; nPerm=5
    dat=arrayData2[1:40,1:100]; nClust=2; stdDev=0.2; nPerm=5
    dat=arrayData2[1:40,1:100]; nClust=2; stdDev=0.2; nPerm=5
    dat=arrayData2[1:40,1:100]; nClust=2; stdDev=0.5; nPerm=5
    
    id=rownames(dat)
    
    datThis=dat
    #distMat=as.dist(1 - cor(datThis,method=distMethod,use="complete.obs"))

    distMat=dist(datThis)
    clustThis=hclust(distMat, method=linkMethod)
    clustObs=cutree(clustThis,k=nClust)
    
    if (F) {
        n=length(clustObs)*(length(clustObs)-1)/2
        specPair=matrix(nrow=n,ncol=nPerm)
        tmp=rep(NA,n)
        annSpecPair=data.frame(spec1=tmp,spec2=tmp,stringsAsFactors=F)
        for (j1 in 1:length(clustObs)) {
            for (j1 in 1:length(clustObs)) {
                
            }
        }

        lapply(clustObs,function(x) {})
    }

    clustPerm=matrix(nrow=length(clustObs),ncol=nPerm)
    rownames(clustPerm)=names(clustObs)
    set.seed(5354)
    for (iPerm in 1:nPerm) {
        datThis2=datThis+rnorm(n=nrow(datThis)*ncol(datThis),mean=0,sd=stdDev)
        #print(summary(c(datThis2)))
        distMat=dist(datThis2)
        clustThis=hclust(distMat, method=linkMethod)
        clustPerm[,iPerm]=cutree(clustThis,k=nClust)
    }
    save(clustPerm,file=paste("clustPerm",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,"_",permSetFlag,".RData",sep=""))
    
    out=matrix(nrow=length(clustObs),ncol=length(clustObs))
    for (iPerm in 1:nPerm) {
    }

    if (F) {
        clustId=cutree(clustThis,k=nClust)[clustThis$order]
        k1=which(!duplicated(clustId))
        for (k in 1:length(k1)) {
            clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
        }
    }

    ## R-index: Measures proportion of pairs of specimens within a cluster
    ## for which the members of the pair remain together in the re-clustered
    ## purturbed data

    ## For a singleton cluster in the original data,
    ## it can be informative to record the proportion
    ## of times it remains a singleton as opposed to
    ## being merged into another cluster in the
    ## perturbed data
}


####################################################################
####################################################################
## Get cluster designations

cohort="_mycRasWt"; typeFlag=""; type2Flag="_reducedBioFeatPC"; subsetFlag=""
distMethod="spearman"
distMethod="kendall"
distMethod="pearson"
linkMethod="ward.D2"
load(file=paste("arrayData2",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,".RData",sep=""))

datThis=arrayData2[1:40,1:100]; nClust=15
datThis=arrayData2; nClust=15

cat("\n\n=========",paste("clustObs",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,"_",distMethod,".RData",sep=""),"========\n")

clustObs=rep(NA,ncol(datThis))
names(clustObs)=colnames(datThis)
#set.seed(5354)
timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))
datThis2=datThis
#print(summary(c(datThis2)))
distMat=as.dist(1 - cor(datThis2,method=distMethod,use="complete.obs"))
clustThis=hclust(distMat, method=linkMethod)
x=rep("",length(clustObs))
for (k in 2:nClust) {
    x=paste(x,cutree(clustThis,k=k),sep=",")
}
x=sub(",","",x)
clustObs=x
names(clustObs)=colnames(datThis)
timeStamp=c(timeStamp,Sys.time())
print(format(timeStamp[2], "%x %X"))
print(diff(timeStamp))
save(clustObs,file=paste("clustObs",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,"_",distMethod,".RData",sep=""))


####################################################################
####################################################################
## Get cluster designations of perturbed data

permSetFlag="2"; set.seed(2546)
permSetFlag="3"; set.seed(3677)
permSetFlag="4"; set.seed(4232)
permSetFlag="5"; set.seed(5243)
permSetFlag="6"; set.seed(6923)
permSetFlag="7"; set.seed(7019)
permSetFlag="8"; set.seed(8204)
permSetFlag="9"; set.seed(9283)
permSetFlag="10"; set.seed(10657)
permSetFlag="1"; set.seed(1354)

cohort="_mycRasWt"; typeFlag=""; type2Flag="_reducedBioFeatPC"; subsetFlag=""
distMethod="kendall"
distMethod="spearman"
linkMethod="ward.D2"
load(file=paste("arrayData2",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,".RData",sep=""))

datThis=arrayData2[1:4,1:5]; nClust=2; stdDev=0.1; nPerm=5
datThis=arrayData2[1:40,1:100]; nClust=2; stdDev=0.2; nPerm=5
datThis=arrayData2; nClust=15; stdDev=0.2; nPerm=250
datThis=arrayData2; nClust=15; stdDev=0.1; nPerm=250
datThis=arrayData2[1:40,1:100]; nClust=5; stdDev=0.2; nPerm=5
datThis=arrayData2; nClust=15; stdDev=0.2; nPerm=25

cat("\n\n=========",paste("clustPerm",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,"_",distMethod,"_sd",stdDev,"_",permSetFlag,".RData",sep=""),"========\n")

clustPerm=matrix(nrow=ncol(datThis),ncol=nPerm)
rownames(clustPerm)=colnames(datThis)
#set.seed(5354)
timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))
for (iPerm in 1:nPerm) {
    datThis2=datThis+rnorm(n=nrow(datThis)*ncol(datThis),mean=0,sd=stdDev)
    #print(summary(c(datThis2)))
    distMat=as.dist(1 - cor(datThis2,method=distMethod,use="complete.obs"))
    clustThis=hclust(distMat, method=linkMethod)
    x=rep("",nrow(clustPerm))
    for (k in 2:nClust) {
        x=paste(x,cutree(clustThis,k=k),sep=",")
    }
    x=sub(",","",x)
    clustPerm[,iPerm]=x
}
timeStamp=c(timeStamp,Sys.time())
print(format(timeStamp[2], "%x %X"))
print(diff(timeStamp))
save(clustPerm,file=paste("clustPerm",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,"_",distMethod,"_sd",stdDev,"_",permSetFlag,".RData",sep=""))

x1=apply(clustPerm,1,function(x) {mean(x!=x[1])})
summary(x1) ## Some should be != 0

####################################################################
####################################################################
## Get cluster designations of subset of observed data

nFold=2

cohort="_mycRasWt"; typeFlag=""; type2Flag="_reducedBioFeatPC"; subsetFlag=""
distMethod="kendall"
distMethod="spearman"
linkMethod="ward.D2"

load(file=paste("arrayData2",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,".RData",sep=""))
load(file=paste("clustObs",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,"_",distMethod,".RData",sep=""))

datThis=arrayData2; nClust=15
datThis=arrayData2; nClust=5

cat("\n\n=========",paste("clustFold",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,"_",distMethod,"_",nFold,"fold",sep=""),"========\n")

clustFold=matrix(nrow=ncol(datThis),ncol=nFold)
rownames(clustFold)=colnames(datThis)
#set.seed(5354)
timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))
for (iFold in 1:nFold) {
    x=rep("",nrow(clustFold))
    for (k in 2:nClust) {
        x0=sapply(clustObs,function(x) as.integer(strsplit(x,",")[[1]][nClust-1]),USE.NAMES=F)
        j0=c()
        for (iClust in 1:k) {
            j=which(x0==iClust)
            samId=round(seq(1,length(j)+1,length=nFold+1)); samId=cbind(samId[1:(length(samId)-1)],samId[2:length(samId)]-1)
            set.seed(1354)
            j=sample(j,replace=F)
            j0=c(j0,j[samId[iFold,1]:samId[iFold,2]])
        }
        datThis2=datThis[,j0]
        #print(summary(c(datThis2)))
        distMat=as.dist(1 - cor(datThis2,method=distMethod,use="complete.obs"))
        clustThis=hclust(distMat, method=linkMethod)
        x=paste(x,cutree(clustThis,k=k),sep=",")
    }
    x=sub(",","",x)
    clustFold[,iFold]=x
}
timeStamp=c(timeStamp,Sys.time())
print(format(timeStamp[2], "%x %X"))
print(diff(timeStamp))
save(clustFold,file=paste("clustFold",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,"_",distMethod,"_",nFold,"fold.RData",sep=""))

x1=apply(clustFold,1,function(x) {mean(x!=x[1])})
summary(x1) ## Some should be != 0

####################################################################
####################################################################
