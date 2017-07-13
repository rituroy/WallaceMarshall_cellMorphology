####################################################################
####################################################################
## Get cluster designations

cohort="_mycRasWt"; typeFlag=""; type2Flag="_reducedBioFeatPC"; subsetFlag=""
distMethod="spearman"
distMethod="kendall"
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
