minCPS=10; maxCPS=15; nPermPS=100
minCPS=10; maxCPS=15; nPermPS=10
minCPS=10; maxCPS=15; nPermPS=20
minCPS=2; maxCPS=3; nPermPS=3
minCPS=2; maxCPS=15; nPermPS=100


####################################################################
####################################################################

#The prediction strength for each of the clusterings is the mean (in prediction.strength(), it is minimum) (taken over all clusters) relative frequency of correctly predicted pairs of points of that cluster.
prediction.strengthMean=function (xdata, Gmin = 2, Gmax = 10, M = 50, clustermethod = kmeansCBI,
classification = "centroid", cutoff = 0.8, nnk = 1, distances = inherits(xdata,
"dist"), count = FALSE, ...)
{
    xdata <- as.matrix(xdata)
    n <- nrow(xdata)
    nf <- c(floor(n/2), n - floor(n/2))
    indvec <- clcenters <- clusterings <- jclusterings <- classifications <- list()
    corrpred <- list()
    for (k in Gmin:Gmax) {
        if (count)
        cat(k, " clusters\n")
        corrpred[[k]] <- numeric(0)
        for (l in 1:M) {
            nperm <- sample(n, n)
            if (count)
            cat(" Run ", l, "\n")
            indvec[[l]] <- list()
            indvec[[l]][[1]] <- nperm[1:nf[1]]
            indvec[[l]][[2]] <- nperm[(nf[1] + 1):n]
            for (i in 1:2) {
                if (distances)
                clusterings[[i]] <- clustermethod(as.dist(xdata[indvec[[l]][[i]],
                indvec[[l]][[i]]]), k, ...)
                else clusterings[[i]] <- clustermethod(xdata[indvec[[l]][[i]],
                ], k, ...)
                jclusterings[[i]] <- rep(-1, n)
                jclusterings[[i]][indvec[[l]][[i]]] <- clusterings[[i]]$partition
                centroids <- NULL
                if (classification == "centroid") {
                    if (identical(clustermethod, kmeansCBI))
                    centroids <- clusterings[[i]]$result$centers
                    if (identical(clustermethod, claraCBI))
                    centroids <- clusterings[[i]]$result$medoids
                }
                j <- 3 - i
                if (distances)
                classifications[[j]] <- classifdist(as.dist(xdata),
                jclusterings[[i]], method = classification,
                centroids = centroids, nnk = nnk)[indvec[[l]][[j]]]
                else classifications[[j]] <- classifnp(xdata,
                jclusterings[[i]], method = classification,
                centroids = centroids, nnk = nnk)[indvec[[l]][[j]]]
            }
            ps <- matrix(0, nrow = 2, ncol = k)
            for (i in 1:2) {
                #ctable <- table(clusterings[[i]]$partition, classifications[[i]])
                ctable <- xtable(clusterings[[i]]$partition, classifications[[i]])
                for (kk in 1:k) {
                    ps[i, kk] <- sum(ctable[kk, ]^2 - ctable[kk,
                    ])
                    cpik <- clusterings[[i]]$partition == kk
                    nik <- sum(cpik)
                    if (nik > 1)
                    ps[i, kk] <- ps[i, kk]/(nik * (nik - 1))
                    else ps[i, kk] <- 1
                }
            }
            corrpred[[k]][l] <- mean(c(mean(ps[1, ]), mean(ps[2,
            ])))
        }
    }
    mean.pred <- numeric(0)
    if (Gmin > 1)
    mean.pred <- c(1)
    if (Gmin > 2)
    mean.pred <- c(mean.pred, rep(NA, Gmin - 2))
    for (k in Gmin:Gmax) mean.pred <- c(mean.pred, mean(corrpred[[k]]))
    optimalk <- max(which(mean.pred > cutoff))
    out <- list(predcorr = corrpred, mean.pred = mean.pred, optimalk = optimalk,
    cutoff = cutoff, method = clusterings[[1]]$clustermethod,
    Gmax = Gmax, M = M)
    class(out) <- "predstr"
    out
}

####################################################################
####################################################################
load("tmp.RData")

library(marray)
library(fpc)
#source(paste(dirSrc,"functions/heatmap.5.R",sep=""))
#source(paste(dirSrc,"functions/heatmapAcgh.7.R",sep=""))

cohort="_mycRas"
cohort="_mycRasWt"
cohort="_wt"
cohort="_wt906"

classifFlag="averagedist"
classifFlag="centroid"

predStrFlag=""
predStrFlag="_mean"

if (predStrFlag=="_mean") {
    predictionStrengthThis=prediction.strengthMean
} else {
    predictionStrengthThis=prediction.strength
}

datadir="results/"
datadir=""

typeFlag="cell"

sepFClustFlag=F
sepFClustFlag=T ## NOT USED for typeFlag=""

centrFlag="_noCentering"
centrFlag=""

scaleList=c("","_noScaling")
scaleList=c("")
scaleFlag="_noScaling"
scaleFlag=""

type2Flag="_noisyFeat"; typeList=""; centrFlag=""; scaleList=""
type2Flag="_reducedFeatPC_PC"; typeList=sort(unique(ann$type)); centrFlag=""; scaleList=""
type2Flag="_reducedFeatPC_PC"; typeList=sort(unique(ann$type)); centrFlag="_noCentering"; scaleList="_noScaling"
type2Flag="_reducedFeatPC"; typeList=sort(unique(ann$type)); centrFlag="_noCentering"; scaleList="_noScaling"
type2Flag="_reducedFeatMean"; typeList=sort(unique(ann$type)); centrFlag="_noCentering"; scaleList="_noScaling"
type2Flag="_reducedFeatMean"; typeList=sort(unique(ann$type)); centrFlag=""; scaleList=""
type2Flag="_reducedFeatPC"; typeList=sort(unique(ann$type)); centrFlag=""; scaleList=""
type2Flag="_zernike1Only_reducedFeatPC"; typeList=sort(unique(ann$type)); centrFlag=""; scaleList=""
type2Flag="_zernike1Only_reducedFeatMean"; typeList=sort(unique(ann$type)); centrFlag=""; scaleList=""
## NOT working. Check!!! type2Flag=""; typeList=sort(unique(ann$type)); centrFlag=""; scaleList=""
type2Flag="_reducedBioFeatMean"; typeList=c("",sort(unique(ann$type))); centrFlag=""; scaleList=""
type2Flag="_reducedBioFeatPC"; typeList=c("",sort(unique(ann$type))); centrFlag=""; scaleList=""
typeList=""
typeList=c("",sort(unique(ann$type)))
typeList=sort(unique(ann$type))
type2Flag=""; typeList=c("",sort(unique(ann$type)))

outFormat="pdf"
outFormat="png"

sampleBar=""
sampleBar="cluster"

distMethod="spearman"
distMethod="pearson"
distMethod="kendall"

linkMethod="ward.D2"

colList=c("skyblue","blue","yellow","purple","red")
colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
colList2=c("skyblue","blue")
colHM=c("red","blue","grey")

limSdSam=c(20,270)
limSdSam=c(20,300)

subsetFlag=""
subsetList=c("_zernike1Only")
subsetList=c("")

getCorFlag=T
getCorFlag=F

orderFlag="_mycRasOrd"
orderFlag="_wtOrd"
orderFlag=""
datadir2=paste("results/heatmap/",sub("_","",sub("Ord","",orderFlag)),"/",sub("_","",orderFlag),"/",sub("_","",type2Flag),"/",sep="")

#for (cohort in c("_mycRas","_mycRasWt","_wt")) {
#for (cohort in c("_mycRasWt")) {
for (cohort in c("_wt906")) {
    for (sepFClustFlag in c(F)) {
    #for (sepFClustFlag in c(F,T)) {
        if (sepFClustFlag) {
            typeList=sort(unique(ann$type))
        } else {
            typeList=""
        }

        switch(cohort,
        "_mycRas"={
            cohortName="Myc/Ras"
            cell=cellM
            cell_rbf=cell_rbfM
            cell_rbfm=cell_rbfmM
            sdSamPc=sdSamPcM
            sdSamMn=sdSamMnM
            sdFeatPc=sdFeatPcM
            sdFeatMn=sdFeatMnM
        },
        "_wt"={
            cohortName="Wildtype"
            cell=cellW
            cell_rbf=cell_rbfW
            cell_rbfm=cell_rbfmW
            sdSamPc=sdSamPcW
            sdSamMn=sdSamMnW
            sdFeatPc=sdFeatPcW
            sdFeatMn=sdFeatMnW
        },
        "_mycRasWt"={
            cohortName="Myc/Ras & Wildtype"
            cell=cellMW
            cell_rbf=cell_rbfMW
            cell_rbfm=cell_rbfmMW
            sdSamPc=sdSamPcMW
            sdSamMn=sdSamMnMW
            sdFeatPc=sdFeatPcMW
            sdFeatMn=sdFeatMnMW
        },
        "_wt906"={
            cohortName="Wildtype 906"
            cell=cellW2
            annCell=annCellW2
            cell_rbf=cell_rbfW2
            cell_rbfm=cell_rbfmW2
            sdSamPc=sdSamPcW2
            sdSamMn=sdSamMnW2
            sdFeatPc=sdFeatPcW2
            sdFeatMn=sdFeatMnW2
        },
        "_mycRas289"={
            cohortName="Myc/Ras 289"
            cell=cellM2
            annCell=annCellM2
            cell_rbf=cell_rbfM2
            cell_rbfm=cell_rbfmM2
            sdSamPc=sdSamPcM2
            sdSamMn=sdSamMnM2
            sdFeatPc=sdFeatPcM2
            sdFeatMn=sdFeatMnM2
        }
        )
        if (type2Flag=="_reducedBioFeatPC") {
            sdFeat=sdFeatPc
            sdSam=sdSamPc
        } else if (type2Flag=="_reducedBioFeatMean") {
            sdFeat=sdFeatMn
            sdSam=sdSamMn
        } else {
            sdFeat=sdSam=rep(NA,nrow(cell))
        }
        if (cohort%in%c("_mycRas","_wt","_mycRasWt")) {
            #ann=annW
            #ann_rbf=ann_rbf1
            #ann_rbfm=ann_rbfm1
        } else {
            #ann=annW2
            #ann_rbf=ann_rbf2
            #ann_rbfm=ann_rbfm2
            featId1=1:nrow(ann)
        }
        
        for (typeFlag in typeList) {
            if (typeFlag=="" & sepFClustFlag) {
                cat("Cannot run typeFlag=='' & sepFClustFlag!\n")
                next
            }
            for (scaleFlag in scaleList) {
                for (subsetFlag in subsetList) {
                    fNameOut=paste(cohort,orderFlag,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,centrFlag,scaleFlag,"_",distMethod,sep="")
                    header=paste(cohortName,typeFlag,sep="")
                    if (typeFlag=="") {
                        header=paste(cohortName,": All samples clustered together",sep="")
                    } else {
                        header=paste(cohortName,": Samples clustered within ",typeFlag,sep="")
                    }
                    #phen2=data.frame(id=paste("sam",1:nrow(cell),sep=""),sd=sdSam,stringsAsFactors=F)
                    phen2=data.frame(id=rownames(cell),sd=sdSam,stringsAsFactors=F)
                    
                    varFList="featRatio"
                    varFName=paste(c("featRatio")," ",sep="")
                    varFList=varFName=NULL
                    varList=c("sd")
                    varName=paste(c("stdDev")," ",sep="")
                    
                    if (cohort=="_mycRasWt") {
                        #phen2=cbind(phen2,cohort=c(rep("Myc/Ras",nrow(cellM)),rep("Wildtype",nrow(cellW))))
                        phen2$cohort="Myc/Ras"
                        phen2$cohort[substr(phen2$id,1,1)=="w"]="Wildtype"
                        varList=c(varList,"cohort")
                        varName=c(varName,paste(c("cohort")," ",sep=""))
                    }
                    
                    if (type2Flag%in%c("_reducedFeatPC_PC","_reducedFeatPC","_reducedFeatMean","_noisyFeat")) {
                        limit1=c(-2,2)
                        nClust=c(NA,NA)
                    } else {
                        if (scaleFlag=="") {
                            limit1=c(-2,2)
                        } else {
                            limit1=c(-4*10^5,4*10^5)
                            limit1=c(-200,200)
                        }
                        nClust=c(3,NA)
                        nClust=c(10,NA)
                        nClust=c(15,NA)
                    }


                    geneBar=""
                    geneBar="clusterPr"
                    if (type2Flag=="") {
        				#featId=featId1[ann$type[featId1]==typeFlag]
                        featId=featId1

                        arrayData=t(cell[,featId])
                        annCol=ann[featId,]

                        if (subsetFlag=="_zernike1Only") {
                            header=paste(header,", no zernike2...",sep="")
                            i=grep("Zernike",annCol$feature)
                            i=i[which(!annCol$feature[i]%in%c("CellZernike1","NucZernike1"))]
                            if (length(i)==0) next
                            arrayData=arrayData[-i,]
                            annCol=annCol[-i,]
                        }
                    } else {
                        switch(type2Flag,
                            "_reducedFeatPC_PC"={
                               header=paste(header,", based on prin comps after PCA based feature reduction",sep="")
                               tbl=read.table(paste("results/pca/prinComp",cohort,"_",typeFlag,"_reducedFeatPCA.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                               arrayData=t(as.matrix(tbl[,-1]))
                               annCol=data.frame(feature=rownames(arrayData),cohort=cohort,type=typeFlag,stringsAsFactors=F)
                               nClust=c(3,NA)
                            },
                            "_reducedFeatPC"={
                               header=paste(header,", PC based feature reduction",sep="")
                               arrayData=t(as.matrix(cell_rf[,which(ann_rf$type==typeFlag)]))
                               if (centrFlag=="") {
                                i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
                                arrayData[i2,]=t(cell[,i1])
                                limit1=c(-2,2)
                               } else {
                               limit1=c(-0.5,0.5)
                               if (typeFlag%in%c("mitochondria","nucleus")) {
                                limit1=c(-0.0005,0.0005)
                               }
                               }
                               annCol=data.frame(feature=rownames(arrayData),cohort=cohort,type=typeFlag,stringsAsFactors=F)
                               nClust=c(3,NA)
                            },
                            "_zernike1Only_reducedFeatPC"={
                               if (typeFlag%in%c("mitochondria")) next
                               header=paste(header,", no zernike2..., prin comp based feature reduction",sep="")
                               arrayData=t(as.matrix(cell_nzrf[,which(ann_nzrf$type==typeFlag)]))
                               if (centrFlag=="") {
                               i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
                               arrayData[i2,]=t(cell[,i1])
                               limit1=c(-2,2)
                               } else {
                               limit1=c(-0.5,0.5)
                               if (typeFlag%in%c("mitochondria","nucleus")) {
                               limit1=c(-0.0005,0.0005)
                               }
                               }
                               annCol=data.frame(feature=rownames(arrayData),cohort=cohort,type=typeFlag,stringsAsFactors=F)
                               nClust=c(3,NA)
                            },
                            "_reducedFeatMean"={
                               header=paste(header,", mean based feature reduction",sep="")
                               arrayData=t(as.matrix(cell_rfm[,which(ann_rfm$type==typeFlag)]))
                               if (centrFlag=="") {
                               i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
                               arrayData[i2,]=t(cell[,i1])
                               limit1=c(-2,2)
                               } else {
                               limit1=c(-0.5,0.5)
                               if (typeFlag%in%c("mitochondria","nucleus")) {
                               limit1=c(-0.0005,0.0005)
                               }
                               }
                               annCol=data.frame(feature=rownames(arrayData),cohort=cohort,type=typeFlag,stringsAsFactors=F)
                               nClust=c(3,NA)
                            },
                            "_zernike1Only_reducedFeatMean"={
                               if (typeFlag%in%c("mitochondria")) next
                               header=paste(header,", no zernike2..., mean based feature reduction",sep="")
                               arrayData=t(as.matrix(cell_nzrfm[,which(ann_nzrfm$type==typeFlag)]))
                               if (centrFlag=="") {
                               i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
                               arrayData[i2,]=t(cell[,i1])
                               limit1=c(-2,2)
                               } else {
                               limit1=c(-0.5,0.5)
                               if (typeFlag%in%c("mitochondria","nucleus")) {
                               limit1=c(-0.0005,0.0005)
                               }
                               }
                               annCol=data.frame(feature=rownames(arrayData),cohort=cohort,type=typeFlag,stringsAsFactors=F)
                               nClust=c(3,NA)
                            },
                            "_reducedBioFeatPC"={
                                #header=paste(header,", prin comp based bio feature reduction",sep="")
                                #header=paste(header,", PC based bio feature reduction",sep="")
                               header=paste(header,", features reduced by PCA",sep="")
                               if (sepFClustFlag | typeFlag=="") {
                                    geneBar="clusterPr"
                               } else {
                                    geneBar=""
                               }
                               arrayData=t(cell_rbf)
                               if (centrFlag=="") {
                                   i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
                                   arrayData[i2,]=t(cell[,i1])
                                    limit1=c(-2,2)
                               } else {
                                    limit1=c(-0.5,0.5)
                               }
                               annCol=ann_rbf
                               nClust=c(ifelse(geneBar=="",NA,3),NA)
                               if (is.null(varFList)) {
                                   varFList=c("type")
                                   varFName=paste(c("type")," ",sep="")
                               } else {
                                   varFList=c(varFList,c("type"))
                                   varFName=c(varFName,paste(c("type")," ",sep=""))
                               }
                               },
                            "_reducedBioFeatMean"={
                                #header=paste(header,", mean based bio feature reduction",sep="")
                               header=paste(header,", features reduced by averaging",sep="")
                               if (sepFClustFlag | typeFlag=="") {
                                geneBar="clusterPr"
                               } else {
                                geneBar=""
                               }
                               arrayData=t(cell_rbfm)
                               if (centrFlag=="") {
                                   i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
                                   arrayData[i2,]=t(cell[,i1])
                                   limit1=c(-2,2)
                               } else {
                                    limit1=c(-0.5,0.5)
                               }
                               annCol=ann_rbfm
                               nClust=c(ifelse(geneBar=="",NA,3),NA)
                               if (is.null(varFList)) {
                                    varFList=c("type")
                                    varFName=paste(c("type")," ",sep="")
                               } else {
                                   varFList=c(varFList,c("type"))
                                   varFName=c(varFName,paste(c("type")," ",sep=""))
                               }
                            },
                            "_noisyFeat"={
                               header=paste("noisy features",sep="")
                               i=which(substr(ann_rf$feature,1,nchar("cluster"))!="cluster")
                               arrayData=t(as.matrix(cell_rf[,i]))
                               annCol=ann_rf[i,]
                               varFList=c("type")
                               varFName=paste(c("type")," ",sep="")
                            }
                        )
                    }
                    if (orderFlag!="") {
                        geneBar=""
                        nClust[1]=NA
                    }

                    annCol$featRatio=annCol$feature
                    i=match(annCol$feature,ann$feature); i2=which(!is.na(i)); i1=i[i2]
                    annCol$featRatio[i2][which(ann$complex[i1]==0)]=NA
                    annCol$featRatio[is.na(i)]=NA

                    annColAll=annCol
                    phenAll=phen2
                    nClustAll=nClust

                    if (centrFlag=="") {
                        centr=apply(arrayData,1,median,na.rm=T)
                        for (i in 1:nrow(arrayData)) {
                            arrayData[i,]=arrayData[i,]-centr[i]
                        }
                    }
                    
                    if (scaleFlag=="") {
                        #header=paste(header,", scaled data",sep="")
                        scal=apply(arrayData,1,sd,na.rm=T)
                        for (i in 1:nrow(arrayData)) {
                            arrayData[i,]=arrayData[i,]/scal[i]
                        }
                    }
                    
                    limit2=limit1
                    if (quantile(c(arrayData),probs=.9,na.rm=T)<10^-3) {
                        arrayData=arrayData*10^4
                        limit2=limit2*10^4
                    }

                    varFListAll=varFList
                    varFNameAll=varFName
                    varListAll=varList
                    varNameAll=varName

                    cat("\n\n======================",fNameOut,"======================\n")
                    cat("\n\n============",paste("predictionStrength_",classifFlag,predStrFlag,"_",minCPS,"clustMin_",maxCPS,"clustMax_",nPermPS,"perms",sep=""),"============\n")


                    ## -------------------
                    if (typeFlag=="") arrayData2=arrayData else arrayData2=arrayData[which(annCol$type==typeFlag),]
                    
                    if (sampleBar=="cluster") {
                        fNameOut4=sub(orderFlag,"",fNameOut)
                        switch(distMethod,
                            "pearson"={
                               distMat=as.dist(1 - cor(arrayData2,method=distMethod,use="complete.obs"))
                           },
                           "spearman"={
                               distMat=as.dist(1 - cor(arrayData2,method=distMethod,use="complete.obs"))
                           },
                           "kendall"={
                               if (getCorFlag) {
                                    corMatSam2=cor(arrayData2,method=distMethod,use="complete.obs")
                                    save(corMatSam2,file=paste("corMatSam",fNameOut4,".RData",sep=""))
                               } else {
                                    load(file=paste(datadir,"corMatSam",fNameOut4,".RData",sep=""))
                               }
                               distMat=as.dist(1 - corMatSam2)
                           },
                           "euclidean"={
                               distMat=dist(arrayData2, method=distMethod)
                           }
                        )
                        clustC=hclust(distMat, method=linkMethod)
                    } else {
                        clustC=NA
                        nClust[2]=NA
                    }
                    arrayDataAll=arrayData
                    clustCAll=clustC
                    nClustAll=nClust
                    
                    
                    if (F) {
                        x=arrayData2
                        x=arrayData2[1:10,1:4]
                        y=dist(x)
                        timeStamp=Sys.time()
                        print(format(timeStamp, "%x %X"))
                        set.seed(98765)
                        res=predictionStrengthThis(,2,3,M=3,clustermethod=hclustCBI,classification=classifFlag,method="ward.D2",scaling=FALSE)
                        timeStamp=c(timeStamp,Sys.time())
                        print(format(timeStamp[2], "%x %X"))
                        print(diff(timeStamp))
                        res=hclustCBI(data=x,k=3,cut="number",method="ward.D2",scaling=F,noisecut=0)
                        res=predictionStrengthThis(x,2,3,M=3,clustermethod=hclustCBI,method="ward.D2",scaling=TRUE)
                        res=predictionStrengthThis(x,2,3,M=3,clustermethod=claraCBI)
                    }
                    
                    timeStamp=Sys.time()
                    print(format(timeStamp, "%x %X"))
                    set.seed(98765)
                    #distMat=dist(t(arrayData2[1:20,1:50]))
                    res=predictionStrengthThis(distMat,minCPS,maxCPS,M=nPermPS,clustermethod=hclustCBI,classification=classifFlag,method="ward.D2",scaling=FALSE)
                    timeStamp=c(timeStamp,Sys.time())
                    print(format(timeStamp[2], "%x %X"))
                    print(diff(timeStamp))
                    save(res,file=paste("predictionStrengthSam_",classifFlag,predStrFlag,"_",minCPS,"clustMin_",maxCPS,"clustMax_",nPermPS,"perms",fNameOut4,".RData",sep=""))


                    ## -------------------
                    if (cohort!="_mycRasWt") {
                        if (sepFClustFlag) {
                            type3List=sort(unique(annColAll$type))
                        } else {
                            type3List=typeFlag
                        }
                        #par(mfrow=c(length(type3List),1))
                        for (type3Flag in type3List) {
                            if (sepFClustFlag & type2Flag%in%c("_reducedBioFeatPC","_reducedBioFeatMean")) {
                                fNameOut2=paste(fNameOut,"_",type3Flag,sep="")
                            } else {
                                fNameOut2=fNameOut
                            }
                            fNameOut4=fNameOut2
                            if (typeFlag=="") i=1:nrow(annColAll) else i=which(annColAll$type==type3Flag)
                            arrayData=arrayData2=arrayDataAll[i,]
                            annCol=annColAll[i,]
                            #if (!sepFClustFlag & type2Flag%in%c("_reducedBioFeatPC","_reducedBioFeatMean")) {
                            if (!sepFClustFlag & typeFlag!="" & type2Flag%in%c("","_reducedBioFeatPC","_reducedBioFeatMean")) {
                                arrayData=arrayDataAll
                                annCol=annColAll
                                x=c()
                                for (type4Flag in sort(unique(annColAll$type))) {
                                    fNameOut4=paste(fNameOut,"_",type4Flag,sep="")
                                    arrayData2=arrayData[which(annCol$type==type4Flag),]
                                    if (orderFlag%in%c("_mycRasOrd","_wtOrd")) {
                                        #fNameOut4=sub(cohort,"",sub(orderFlag,"",fNameOut2))
                                        fNameOut4=sub(cohort,"",sub("Ord","",fNameOut2))
                                        #fNameOut4=sub("spearman","kendall",sub(cohort,"",sub(orderFlag,"",fNameOut2)))
                                        clustInfo=read.table(paste(datadir2,fNameOut4,"/clusterInfoFeature",paste(fNameOut4,"_",type4Flag,sep=""),".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                        x=c(x,annCol$feature[match(clustInfo$feature,annCol$feature)])
                                    } else {
                                    switch(distMethod,
                                        "pearson"={
                                           distMat=as.dist(1 - cor(t(arrayData2),method=distMethod,use="complete.obs"))
                                        },
                                        "spearman"={
                                           distMat=as.dist(1 - cor(t(arrayData2),method=distMethod,use="complete.obs"))
                                        },
                                        "kendall"={
                                           if (getCorFlag) {
                                                corMat2=cor(t(arrayData2),method=distMethod,use="complete.obs")
                                                save(corMat2,file=paste("corMat",fNameOut4,".RData",sep=""))
                                           } else {
                                                load(file=paste(datadir,"corMat",fNameOut4,".RData",sep=""))
                                           }
                                           distMat=as.dist(1 - corMat2)
                                        },
                                        "euclidean"={
                                           distMat=dist(t(arrayData2), method=distMethod)
                                        }
                                    )
                                    clustR=hclust(distMat, method=linkMethod)
                                    x=c(x,clustR$labels[clustR$order])
                                    }
                                }
                                i=match(x,annCol$feature)
                                arrayData=arrayData[i,]
                                annCol=annCol[i,]
                                clustR=NA
                                nClust[1]=NA
                            } else {
                                if (orderFlag%in%c("_mycRasOrd","_wtOrd")) {
                                    #fNameOut4=sub(cohort,"",sub(orderFlag,"",fNameOut))
                                    fNameOut4=sub(cohort,"",sub("Ord","",fNameOut))
                                    clustInfo=read.table(paste(datadir2,fNameOut4,"/clusterInfoFeature",fNameOut4,ifelse(type3Flag=="","","_"),type3Flag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                    i=match(clustInfo$feature,annColAll$feature)
                                    arrayData=arrayDataAll[i,]
                                    annCol=annColAll[i,]
                                    clustR=NA
                                    nClust[1]=NA
                                    if (type3Flag==typeFlag) {
                                        clustC=clustCAll
                                        nClust[2]=nClustAll[2]
                                    } else {
                                        j=clustCAll$order
                                        arrayData=arrayData[,j]
                                        clustC=NA
                                        nClust[2]=NA
                                    }
                                } else {
                                    if (sepFClustFlag & type2Flag%in%c("_reducedBioFeatPC","_reducedBioFeatMean")) {
                                        if (type3Flag==typeFlag) {
                                            clustC=clustCAll
                                            nClust[2]=nClustAll[2]
                                        } else {
                                            i=which(annColAll$type==type3Flag)
                                            j=clustCAll$order
                                            arrayData=arrayData2=arrayDataAll[i,j]
                                            annCol=annColAll[i,]
                                            clustC=NA
                                            nClust[2]=NA
                                        }
                                    }
                                    if (geneBar=="clusterPr") {
                                        switch(distMethod,
                                               "pearson"={
                                               distMat=as.dist(1 - cor(t(arrayData2),method=distMethod,use="complete.obs"))
                                               },
                                               "spearman"={
                                               distMat=as.dist(1 - cor(t(arrayData2),method=distMethod,use="complete.obs"))
                                               },
                                               "kendall"={
                                               if (getCorFlag) {
                                                   corMat2=cor(t(arrayData2),method=distMethod,use="complete.obs")
                                                   save(corMat2,file=paste("corMat",fNameOut4,".RData",sep=""))
                                               } else {
                                                    load(file=paste(datadir,"corMat",fNameOut4,".RData",sep=""))
                                               }
                                               distMat=as.dist(1 - corMat2)
                                               },
                                               "euclidean"={
                                               distMat=dist(t(arrayData2), method=distMethod)
                                               }
                                               )
                                        clustR=hclust(distMat, method=linkMethod)
                                    } else {
                                        clustR=NA
                                        nClust[1]=NA
                                    }
                                }
                            }
                            
                            if (type3Flag==typeFlag) {
                                timeStamp=Sys.time()
                                print(format(timeStamp, "%x %X"))
                                set.seed(98765)
                                #distMat=dist(arrayData2[1:20,1:50])
                                res=predictionStrengthThis(distMat,minCPS,maxCPS,M=nPermPS,clustermethod=hclustCBI,classification=classifFlag,method="ward.D2",scaling=FALSE)
                                timeStamp=c(timeStamp,Sys.time())
                                print(format(timeStamp[2], "%x %X"))
                                print(diff(timeStamp))
                                save(res,file=paste("predictionStrength_",classifFlag,predStrFlag,"_",minCPS,"clustMin_",maxCPS,"clustMax_",nPermPS,"perms",fNameOut4,".RData",sep=""))
                            }
                        }
                    }
                }

            }
        }
    }
}

####################################################################
####################################################################
