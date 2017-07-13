argv=commandArgs(TRUE)
cohort=argv[1]
distMethod=argv[2]
stdDev=argv[3]
nPerm=argv[4]
permSetFlag=argv[5]

stdDev=as.numeric(stdDev)
nPerm=as.integer(nPerm)

if (F) {
    cohort="_mycRasWt"
    stdDev=0.2
    distMethod="kendall"
    distMethod="pearson"
    permSetFlag="1"
}


####################################################################
####################################################################
## Get cluster designations of perturbed data


switch(permSetFlag,
"1"={set.seed(1354)},
"2"={set.seed(2546)},
"3"={set.seed(3677)},
"4"={set.seed(4232)},
"5"={set.seed(5243)},
"6"={set.seed(6923)},
"7"={set.seed(7019)},
"8"={set.seed(8204)},
"9"={set.seed(9283)},
"10"={set.seed(10657)}
)

typeFlag=""; type2Flag="_reducedBioFeatPC"; subsetFlag=""
linkMethod="ward.D2"
load(file=paste("arrayData2",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,".RData",sep=""))

if (F) {
datThis=arrayData2[1:4,1:5]; nClust=2; stdDev=0.1; nPerm=5
datThis=arrayData2; nClust=15; stdDev=0.2; nPerm=250
datThis=arrayData2[1:40,1:100]; nClust=5; stdDev=0.2; nPerm=1
datThis=arrayData2; nClust=15; stdDev=0.2; nPerm=25
}
datThis=arrayData2; nClust=15
#datThis=arrayData2[1:40,1:100]; nClust=5

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

####################################################################
####################################################################
