 ## Set params:
## sepFClustFlag, getCorFlag, orderFlag, typeList, type2Flag, cohortList, datVerList

setwd("/Users/royr/UCSF/WallaceMarshall")

####################################################################
####################################################################
##load("results/tmp_w3830.RData"); datVer="w3830"
#load("results/tmp_w5969.RData")
load("data.RData")
#load("tmp_wt906mycRas289.RData")
datVer="w906"
datVerList=c("w906","w5969")
datVerList=c("w5969")
datVerList=c("w906")
datVerList=c("w906","m289")
datVerList=c("w906m289")
for (datVer in datVerList) {
    if (datVer=="w5969") {load("results/tmp_w5969.RData"); datVer="w5969"}
    
    exampleFlag=T
    exampleFlag=F
    
    sepFClustFlag=T ## NOT USED for typeFlag=""
    sepFClustFlag=F
    
    getCorFlag=F
    getCorFlag=T
    
    library(marray)
    #library(fpc)
    library(NbClust)
    #source(paste(dirSrc,"functions/heatmap.5.R",sep=""))
    #source(paste(dirSrc,"functions/heatmapAcgh.7.R",sep=""))
    #source(paste(dirSrc,"functions/heatmap.5.2.R",sep=""))
    #source(paste(dirSrc,"functions/heatmap.5.4.R",sep=""))
    #source(paste(dirSrc,"functions/heatmap.5.5.R",sep=""))
    #source(paste(dirSrc,"functions/heatmapAcgh.7.1.R",sep=""))
    #source(paste(dirSrc,"functions/heatmap.5.6.R",sep=""))
    #source(paste(dirSrc,"functions/heatmapAcgh.7.3.R",sep=""))
    source(paste(dirSrc,"functions/heatmap.5.7.R",sep=""))
    source(paste(dirSrc,"functions/heatmapAcgh.7.4.R",sep=""))
    
    cohort="_mycRas"
    cohort="_wt"
    cohort="_mycRasWt"
    cohortList=c("_mycRasWt","_wt")
    cohortList=c("_mycRas")
    cohortList=c("_wt")
    cohortList=c("_mycRasWt")

    subsetFFlag="_unselectedFeat"
    subsetFFlag="_selectedFeat"
    subsetFFlag=""
    subsetFList=c("_unselectedFeat","_selectedFeat")
    subsetFList=""

    predictionStrengthFlag=T
    predictionStrengthFlag=F
    classifFlag="averagedist"
    classifFlag="centroid"
    nbClustFlag=T
    nbClustFlag=F
    
    datadir="results/cor/"
    switch(datVer,
        "w3830"={
            datadir="results/wt3830/cor/"
        },
        "w906"={
            datadir="results/wt906/cor/"
            ann=annW2
        },
        "m289"={
            datadir="results/mycRas289/cor/"
            ann=annW2
        },
        "w906m289"={
            datadir="results/wt906mycRas289/cor/"
            ann=annW2
        }
    )
    datadir=""
    datadir="results/wt906mycRas289/"
    datadir="results/wt906mycRas289_noCorr/"
    
    typeFlag="cell"
    
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
    type2Flag=""; typeList=c("",sort(unique(ann$type))); centrFlag=""; scaleList=""
    typeList=sort(unique(ann$type))
    typeList=""
    typeList=c("",sort(unique(ann$type)))
    
    sampleBar=""
    sampleBar="cluster"

    distMethod="spearman"
    distMethod="kendall"
    distMethod="pearson"
    
    classDistFlag="euclidean"
    classDistFlag="kendall"
    classDistFlag=""
    
    linkMethod="ward.D2"

    centrFlag="_noCentering"
    centrFlag="_centered"
    centrFlag=""
    centrFlag="_centerByType"
    
    scaleList=c("","_noScaling")
    scaleList=c("_scaled")
    scaleList=c("_scaleByType")
    scaleList=c("")
    scaleFlag="_noScaling"
    scaleFlag="_scaled"
    scaleFlag="_scaleByType"
    scaleFlag=""
    
    outFormat="pdf"
    outFormat="png"

    colList=c("skyblue","blue","yellow","purple","red")
    colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
    colList2=c("skyblue","blue")
    colHM=c("red","blue","grey")

    limSdSam=c(20,270)
    limSdSam=c(20,300)
    limDist2classSam=c(1,5)

    subsetFlag=""
    subsetList=c("_zernike1Only")
    subsetList=c("_classSet1")
    subsetList=c("_classSet1MissClass")
    subsetList=c("_classSet1MissClassRpart")
    subsetList=c("_classSet1Set2MissClassRpart")
    subsetList=c("_classSet1","_classSet1Set2")
    subsetList=c("")
    
    
    if (any(c("_classSet1MissClassRpart","_classSet1Set2MissClassRpart")%in%subsetList)) annRpart=read.table(paste("rpart.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    
    
    orderFlag=rep("",2)
    datadir2=rep("",2)

    orderFlag[1]="_mycRasOrd"
    orderFlag[1]="_wt_classSet1Ord"
    orderFlag[1]=""
    orderFlag[1]="_wtOrd"
    datadir2[1]=""
    switch(datVer,
        "w3830"={
            datadir2[1]=paste("results/wt3830/heatmap/",sub("_","",sub("Ord","",orderFlag[1])),"/",sub("_","",orderFlag[1]),"/",sub("_","",type2Flag),"/",sep="")
        },
        "w906"={
            datadir2[1]=paste("results/wt906/heatmap/",sub("_","",sub("Ord","",orderFlag[1])),"/",sub("_","",orderFlag[1]),"/",sub("_","",type2Flag),"/",sep="")
        },
        "m289"={
            datadir2[1]=paste("results/mycRas289/heatmap/",sub("_","",sub("Ord","",orderFlag[1])),"/",sub("_","",orderFlag[1]),"/",sub("_","",type2Flag),"/",sep="")
        },
        "w906m289"={
            datadir2[1]=paste("results/wt906mycRas289/heatmap/",sub("_","",sub("Ord","",orderFlag[1])),"/",sub("_","",orderFlag[1]),"/",sub("_","",type2Flag),"/",sep="")
        }
    )

    orderFlag[2]="_kmeans2ClustOrd"
    datadir2[2]=paste("results/heatmap/",sub("_","",cohort),"/kmeans/",sep="")
    switch(datVer,
        "w3830"={
            datadir2[2]=paste("results/wt3830/heatmap/",sub("_","",cohort),"/kmeans/",sep="")
        },
        "w906"={
            datadir2[2]=paste("results/wt906/heatmap/",sub("_","",cohort),"/kmeans/",sep="")
        },
        "m289"={
            datadir2[2]=paste("results/mycRas289/heatmap/",sub("_","",cohort),"/kmeans/",sep="")
        },
        "w906m289"={
            datadir2[2]=paste("results/wt906mycRas289/heatmap/",sub("_","",cohort),"/kmeans/",sep="")
        }
    )
    orderFlag[2]=""
    datadir2[2]=""

    orderSamList=c("_hclust4ClustOrd","_kmeans2ClustOrd","_kmeans3ClustOrd","_kmeans4ClustOrd")
    orderSamList=c("_hclust4ClustOrd")
    orderSamList=c("_kmeans2ClustOrd","_kmeans3ClustOrd","_kmeans4ClustOrd")
    orderSamList=c("_kmeans3ClustOrd")
    orderSamList=c("")

    switch(datVer,
        "w5969"={
            cohortList=c("_wt","_mycRas","_mycRasWt")
            cohortList=c("_wt","_mycRas")
            cohortList=c("_mycRas")
        },
        "w906"={
            cohortList=c("_wt906")
        },
        "m289"={
            cohortList=c("_mycRas289")
        },
        "w906m289"={
            cohortList=c("_wt906mycRas289")
        }
    )
    subsetFList=""
    type2Flag=""; typeList=c("",sort(unique(ann$type)))
    type2Flag="_reducedBioFeatPC"; typeList=c("",sort(unique(ann$type)))
    type2Flag=""; typeList=c("",sort(unique(ann$type)))
    distMethod="pearson"
    distMethod="kendall"
    orderFlag[1]="_wtOrd"; datadir2[1]=paste(sub("_","",sub("Ord","",orderFlag[1])),"/",sub("_","",orderFlag[1]),"/",sub("_","",type2Flag),"/",sep="")
    orderFlag[1]="_mycRas289Ord"; datadir2[1]=""
    orderFlag[1]=""
    
    if (datVer=="m289") {
        orderFlag[1]=""
        orderFlag[1]="_wt906Ord"; datadir2[1]=""
    } else {
        orderFlag[1]="_mycRas289Ord"; datadir2[1]=""
        orderFlag[1]=""
    }
    
    for (orderSam in orderSamList) {
        orderFlag[2]=orderSam
        orderFlag[2]=""
        for (subsetFFlag in subsetFList) {
            
            for (cohort in cohortList) {
                #for (sepFClustFlag in c(F,T)) {
                if (F) {
                    if (sepFClustFlag) {
                        typeList=sort(unique(ann$type))
                    } else {
                        typeList=""
                    }
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
                        if (centrFlag=="_centered") {
                            j=1:nrow(cellW2)
                            cell=doCentering(rbind(cellW2,cellM2))
                            cell_rbf=doCentering(rbind(cell_rbfW2,cell_rbfM2))
                            cell_rbfm=doCentering(rbind(cell_rbfmW2,cell_rbfmM2))
                            if (scaleFlag=="_scaled") {
                                cell=doScaling(cell)[j,]
                                cell_rbf=doScaling(cell_rbf)[j,]
                                cell_rbfm=doScaling(cell_rbfm)[j,]
                            } else {
                                cell=cell[j,]
                                cell_rbf=cell_rbf[j,]
                                cell_rbfm=cell_rbfm[j,]
                            }
                        } else {
                            if (scaleFlag=="_scaled") {
                                j=1:nrow(cellW2)
                                cell=doScaling(rbind(cellW2,cellM2))[j,]
                                cell_rbf=doScaling(rbind(cell_rbfW2,cell_rbfM2))[j,]
                                cell_rbfm=doScaling(rbind(cell_rbfmW2,cell_rbfmM2))[j,]
                            }
                        }
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
                        if (centrFlag=="_centered") {
                            j=nrow(cellW2)+(1:nrow(cellM2))
                            cell=doCentering(rbind(cellW2,cellM2))
                            cell_rbf=doCentering(rbind(cell_rbfW2,cell_rbfM2))
                            cell_rbfm=doCentering(rbind(cell_rbfmW2,cell_rbfmM2))
                            if (scaleFlag=="_scaled") {
                                cell=doScaling(cell)[j,]
                                cell_rbf=doScaling(cell_rbf)[j,]
                                cell_rbfm=doScaling(cell_rbfm)[j,]
                            } else {
                                cell=cell[j,]
                                cell_rbf=cell_rbf[j,]
                                cell_rbfm=cell_rbfm[j,]
                            }
                        } else {
                            if (scaleFlag=="_scaled") {
                                j=nrow(cellW2)+(1:nrow(cellM2))
                                cell=doScaling(rbind(cellW2,cellM2))[j,]
                                cell_rbf=doScaling(rbind(cell_rbfW2,cell_rbfM2))[j,]
                                cell_rbfm=doScaling(rbind(cell_rbfmW2,cell_rbfmM2))[j,]
                            }
                        }
                    },
                    "_wt906mycRas289"={
                        cohortName="Wildtype 906 + Myc/Ras 289"
                        cell=rbind(cellW2,cellM2)
                        annCell=rbind(cbind(annCellW2,type=rep("wildtype",nrow(annCellW2)),stringsAsFactors=F),cbind(annCellM2,type=rep("myc/ras",nrow(annCellM2)),stringsAsFactors=F))
                        cell_rbf=NULL
                        cell_rbfm=NULL
                        sdSamPc=NULL
                        sdSamMn=NULL
                        sdFeatPc=NULL
                        sdFeatMn=NULL
                        if (centrFlag=="_centerByType") {
                            cell1=doCentering(cellW2)
                            cell2=doCentering(cellM2)
                        }
                        if (scaleFlag=="_scaleByType") {
                            cell1=doScaling(cell1)
                            cell2=doScaling(cell2)
                        }
                        cell=rbind(cell1,cell2)
                        rm(cell1,cell2)
                    }
                )
                if (cohort%in%c("_mycRas","_wt","_mycRasWt")) {
                    ann=annW
                    ann_rbf=ann_rbf1
                    ann_rbfm=ann_rbfm1
                } else {
                    ann=annW2
                    ann_rbf=ann_rbf2
                    ann_rbfm=ann_rbfm2
                    featId1=1:nrow(ann)
                }
                #featId1=1:nrow(ann)
                
                if (exampleFlag) {
                    i=order(apply(cell,2,sd,na.rm=T),decreasing=T)
                    j=order(apply(cell,1,sd,na.rm=T),decreasing=T)
                    cell=cell[j,i]
                    ann=ann[i,]
                    annCell=annCell[j,]
                    featId1=1:nrow(ann)
                    
                    grpUniq=unique(ann$type)
                    i=c()
                    for (gId in 1:length(grpUniq)) {
                        i=c(i,which(ann$type==grpUniq[gId])[1:20])
                    }
                    j=1:50
                    cell=cell[j,i]
                    ann=ann[i,]
                    annCell=annCell[j,]
                    featId1=1:nrow(ann)
                }
                
                if (F) {
                    arrayData=cell
                    centr=apply(arrayData,1,median,na.rm=T)
                    for (i in 1:ncol(arrayData)) {
                        arrayData[,i]=arrayData[,i]-centr[i]
                    }
                    scal=apply(arrayData,1,sd,na.rm=T)
                    for (i in 1:ncol(arrayData)) {
                        arrayData[,i]=arrayData[,i]/scal[i]
                    }
                    arrayData2=arrayData
                    for (j in 1:ncol(arrayData)) {
                        meanThis=mean(arrayData[,j],na.rm=T)
                        arrayData[is.na(arrayData[,j]),j]=meanThis
                    }
                    cell=arrayData
                }
                
                if (type2Flag=="_reducedBioFeatPC") {
                    sdFeat=sdFeatPc
                    sdSam=sdSamPc
                } else if (type2Flag=="_reducedBioFeatMean") {
                    sdFeat=sdFeatMn
                    sdSam=sdSamMn
                } else {
                    sdFeat=sdSam=rep(NA,nrow(cell))
                }
                
                for (typeFlag in typeList) {
                    if (typeFlag=="" & sepFClustFlag) {
                        cat("Cannot run typeFlag=='' & sepFClustFlag! Skipped .......\n")
                        next
                    }
                    for (scaleFlag in scaleList) {
                        for (subsetFlag in subsetList) {
                            #fNameOut=paste(cohort,orderFlag[1],orderFlag[2],ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,subsetFFlag,centrFlag,scaleFlag,"_",distMethod,sep="")
                            fNameOut=paste(cohort,subsetFlag,orderFlag[1],orderFlag[2],ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFFlag,centrFlag,scaleFlag,"_",distMethod,sep="")
                            header=paste(cohortName,typeFlag,sep="")
                            if (typeFlag=="") {
                                #header=paste(cohortName,": All samples clustered together",sep="")
                                header=paste(cohortName,": AllSamplesClusteredTogether",sep="")
                            } else {
                                #header=paste(cohortName,": Samples clustered within ",typeFlag,sep="")
                                header=paste(cohortName,": SamplesClusteredIn-",typeFlag,sep="")
                            }
                            nm=""
                            switch(subsetFlag,
                                "_classSet1"={nm=", set1"},
                                "_classSet1Set2"={nm=", set1set2"},
                                "_classSet1MissClass"={nm=", set1MissClassified"},
                                "_classSet1MissClassRpart"={nm=", set1MissClassified"},
                                "_classSet1Set2MissClassRpart"={nm=", set1set2MissClassified"}
                            )
                            if (nm!="") {
                                header=sub(cohortName,paste(cohortName,nm,sep=""),header)
                            }
                            if (classDistFlag!="") {
                                #header=paste(header,", classDist-",classDistFlag,sep="")
                                load("results/class.RData")
                                switch(classDistFlag,
                                    "euclidean"={classDistMat=classDistMatE},
                                    "pearson"={classDistMat=classDistMatP},
                                    "spearman"={classDistMat=classDistMatS},
                                    "kendall"={classDistMat=classDistMatK}
                                )
                            }
                            #annSam=data.frame(id=paste("sam",1:nrow(cell),sep=""),sd=sdSam,stringsAsFactors=F)
                            annSam=data.frame(id=rownames(cell),sd=sdSam,stringsAsFactors=F)
                            j=match(annSam$id,annCell$id); j1=which(!is.na(j)); j2=j[j1]
                            if (length(j1)!=0) {
                                if (cohort%in%c("_wt906mycRas289")) {
                                    annSam$type[j1]=annCell$type[j2]
                                }
                                #annSam$class=""
                                #annSam$class[j1]=annCell$class[j2]
                                annSam$set1Class=""
                                j2=which(annCell$set==1)
                                j=match(annSam$id,annCell$id[j2]); j1=which(!is.na(j)); j2=j2[j[j1]]
                                #annSam$set1Class[j1]=sub("class","",annCell$class[j2])
                                annSam$set1Class[j1]=sub("class","",paste("class",annCell$class[j2],sep=""))
                                annSam$set2Class=""
                                j2=which(annCell$set==2)
                                j=match(annSam$id,annCell$id[j2]); j1=which(!is.na(j)); j2=j2[j[j1]]
                                #annSam$set2Class[j1]=sub("class","",annCell$class[j2])
                                annSam$set2Class[j1]=sub("class","",paste("class",annCell$class[j2],sep=""))
                            }
                            if (classDistFlag!="") {
                                j=match(annSam$id,rownames(classDistMat)); j1=which(!is.na(j)); j2=j[j1]
                                if (any(!is.na(j))) {
                                    tmp=matrix(nrow=nrow(annSam),ncol=ncol(classDistMat),dimnames=list(annSam$id,paste("dist2",colnames(classDistMat),sep="")))
                                    tmp[j1,]=classDistMat[j2,]
                                    tmp[j1,]=t(apply(classDistMat[j2,],1,function(x) {
                                        y=order(order(x))
                                        y[is.na(x)]=NA
                                        y
                                    }))
                                    colnames(tmp)=paste("set1",capWords(colnames(tmp)),sep="")
                                    annSam=cbind(annSam,tmp)
                                }
                            }
                            
                            varFList="featRatio"
                            varFName=paste(c("featRatio")," ",sep="")
                            varFList=varFName=NULL
                            varList=c("sd")
                            varName=paste(c("stdDev")," ",sep="")
                            varList=varName=NULL
                            
                            if (cohort=="_mycRasWt") {
                                #annSam=cbind(annSam,cohort=c(rep("Myc/Ras",nrow(cellM)),rep("Wildtype",nrow(cellW))))
                                annSam$cohort="Myc/Ras"
                                annSam$cohort[substr(annSam$id,1,1)=="w"]="Wildtype"
                                #varList=c(varList,"cohort","set1Class")
                                #varName=c(varName,paste(c("cohort","set1Class")," ",sep=""))
                                varList=c(varList,"cohort")
                                varName=c(varName,paste(c("cohort")," ",sep=""))
                            } else if (cohort=="_wt") {
                                annSam$set1ClassPred=""
                                j=match(annSam$id,colnames(classPredMat)); j1=which(!is.na(j)); j2=j[j1]
                                annSam$set1ClassPred[j1]=sub("class","",classPredMat[paste(cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),sep=""),j2])
                                #varList=c(varList,"set1Class","set1ClassPred")
                                #varName=c(varName,paste(c("set1Class","set1ClassPred")," ",sep=""))
                                if (subsetFlag=="_classSet1Set2") {
                                    varList=c(varList,"set2Class")
                                    varName=c(varName,paste(c("set2Class")," ",sep=""))
                                }
                            } else if (cohort=="_wt906mycRas289") {
                                varList="type"
                                varName="type "
                            }
                            if (F) {
                                j=grep("Dist2class",names(annSam))
                                if (length(j)!=0) {
                                    #nm=names(annSam)[grep("dist2class",names(annSam))]
                                    nm=names(annSam)[j]
                                    varList=c(varList,nm)
                                    varName=c(varName,nm)
                                }
                            }
                            
                            if (type2Flag%in%c("_reducedFeatPC_PC","_reducedFeatPC","_reducedFeatMean","_noisyFeat")) {
                                limit1=c(-2,2)
                                nClust=c(NA,NA)
                            } else {
                                if (scaleFlag%in%c("","_scaled","_scaleByType")) {
                                    limit1=c(-2,2)
                                } else {
                                    limit1=c(-4*10^5,4*10^5)
                                    limit1=c(-200,200)
                                }
                                nClust=c(3,NA)
                                nClust=c(10,NA)
                                nClust=c(15,NA)
                                nClust=c(5,NA)
                            }
                            
                            
                            cat("\n\n=========== 61 ",fNameOut,"\nnClust[2] ",nClust[2],"\n",sep="")
                            
                            geneBar=""
                            geneBar="clusterPr"
                            if (type2Flag=="") {
                                #featId=featId1[ann$type[featId1]==typeFlag]
                                featId=featId1
                                
                                arrayData=t(cell[,featId])
                                annFeat=ann[featId,]
                                
                                if (subsetFlag=="_zernike1Only") {
                                    header=paste(header,", no zernike2...",sep="")
                                    i=grep("Zernike",annFeat$feature)
                                    i=i[which(!annFeat$feature[i]%in%c("CellZernike1","NucZernike1"))]
                                    if (length(i)==0) next
                                    arrayData=arrayData[-i,]
                                    annFeat=annFeat[-i,]
                                }
                                if (sepFClustFlag | typeFlag=="") {
                                    geneBar="clusterPr"
                                } else {
                                    geneBar=""
                                }
                                if (centrFlag=="") {
                                    i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
                                    arrayData[i2,]=t(cell[,i1])
                                    limit1=c(-2,2)
                                } else if (centrFlag%in%c("_centered","_centerByType")) {
                                    limit1=c(-2,2)
                                } else {
                                    limit1=c(-0.5,0.5)
                                }
                                nClust=c(ifelse(geneBar=="",NA,3),NA)
                                if (is.null(varFList)) {
                                    if (any(annFeat$type!="",na.rm=T)) {
                                        varFList=c("type")
                                        varFName=paste(c("type")," ",sep="")
                                    }
                                } else {
                                    varFList=c(varFList,c("type"))
                                    varFName=c(varFName,paste(c("type")," ",sep=""))
                                }
                            } else {
                                switch(type2Flag,
                                "_reducedFeatPC_PC"={
                                    header=paste(header,", based on prin comps after PCA based feature reduction",sep="")
                                    tbl=read.table(paste("results/pca/prinComp",cohort,"_",typeFlag,"_reducedFeatPCA.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                    arrayData=t(as.matrix(tbl[,-1]))
                                    annFeat=data.frame(feature=rownames(arrayData),cohort=cohort,type=typeFlag,stringsAsFactors=F)
                                    nClust=c(3,NA)
                                },
                                "_reducedFeatPC"={
                                    header=paste(header,", PC based feature reduction",sep="")
                                    arrayData=t(as.matrix(cell_rf[,which(ann_rf$type==typeFlag)]))
                                    if (centrFlag=="") {
                                        i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
                                        arrayData[i2,]=t(cell[,i1])
                                        limit1=c(-2,2)
                                    } else if (centrFlag%in%c("_centered","_centerByType")) {
                                        limit1=c(-2,2)
                                    } else {
                                        limit1=c(-0.5,0.5)
                                        if (typeFlag%in%c("mitochondria","nucleus")) {
                                            limit1=c(-0.0005,0.0005)
                                        }
                                    }
                                    annFeat=data.frame(feature=rownames(arrayData),cohort=cohort,type=typeFlag,stringsAsFactors=F)
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
                                    } else if (centrFlag%in%c("_centered","_centerByType")) {
                                        limit1=c(-2,2)
                                    } else {
                                        limit1=c(-0.5,0.5)
                                        if (typeFlag%in%c("mitochondria","nucleus")) {
                                            limit1=c(-0.0005,0.0005)
                                        }
                                    }
                                    annFeat=data.frame(feature=rownames(arrayData),cohort=cohort,type=typeFlag,stringsAsFactors=F)
                                    nClust=c(3,NA)
                                },
                                "_reducedFeatMean"={
                                    header=paste(header,", mean based feature reduction",sep="")
                                    arrayData=t(as.matrix(cell_rfm[,which(ann_rfm$type==typeFlag)]))
                                    if (centrFlag=="") {
                                        i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
                                        arrayData[i2,]=t(cell[,i1])
                                        limit1=c(-2,2)
                                    } else if (centrFlag%in%c("_centered","_centerByType")) {
                                        limit1=c(-2,2)
                                    } else {
                                        limit1=c(-0.5,0.5)
                                        if (typeFlag%in%c("mitochondria","nucleus")) {
                                            limit1=c(-0.0005,0.0005)
                                        }
                                    }
                                    annFeat=data.frame(feature=rownames(arrayData),cohort=cohort,type=typeFlag,stringsAsFactors=F)
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
                                    } else if (centrFlag%in%c("_centered","_centerByType")) {
                                        limit1=c(-2,2)
                                    } else {
                                        limit1=c(-0.5,0.5)
                                        if (typeFlag%in%c("mitochondria","nucleus")) {
                                            limit1=c(-0.0005,0.0005)
                                        }
                                    }
                                    annFeat=data.frame(feature=rownames(arrayData),cohort=cohort,type=typeFlag,stringsAsFactors=F)
                                    nClust=c(3,NA)
                                },
                                "_reducedBioFeatPC"={
                                    #header=paste(header,", prin comp based bio feature reduction",sep="")
                                    #header=paste(header,", PC based bio feature reduction",sep="")
                                    #header=paste(header,", features reduced by PCA",sep="")
                                    header=paste(header,", featReducedByPCA",sep="")
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
                                    } else if (centrFlag%in%c("_centered","_centerByType")) {
                                        limit1=c(-2,2)
                                    } else {
                                        limit1=c(-0.5,0.5)
                                    }
                                    annFeat=ann_rbf
                                    nClust=c(ifelse(geneBar=="",NA,3),NA)
                                    if (is.null(varFList)) {
                                        if (any(annFeat$type!="",na.rm=T)) {
                                            varFList=c("type")
                                            varFName=paste(c("type")," ",sep="")
                                        }
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
                                    } else if (centrFlag%in%c("_centered","_centerByType")) {
                                        limit1=c(-2,2)
                                    } else {
                                        limit1=c(-0.5,0.5)
                                    }
                                    annFeat=ann_rbfm
                                    nClust=c(ifelse(geneBar=="",NA,3),NA)
                                    if (is.null(varFList)) {
                                        if (any(annFeat$type!="",na.rm=T)) {
                                            varFList=c("type")
                                            varFName=paste(c("type")," ",sep="")
                                        }
                                    } else {
                                        varFList=c(varFList,c("type"))
                                        varFName=c(varFName,paste(c("type")," ",sep=""))
                                    }
                                },
                                "_noisyFeat"={
                                    header=paste("noisy features",sep="")
                                    i=which(substr(ann_rf$feature,1,nchar("cluster"))!="cluster")
                                    arrayData=t(as.matrix(cell_rf[,i]))
                                    annFeat=ann_rf[i,]
                                    if (any(annFeat$type!="",na.rm=T)) {
                                        varFList=c("type")
                                        varFName=paste(c("type")," ",sep="")
                                    }
                                }
                                )
                            }
                            if (orderFlag[1]!="") {
                                geneBar=""
                                nClust[1]=NA
                            }
                            if (orderFlag[2]!="") {
                                sampleBar=""
                                nClust[2]=NA
                            }
                            
                            annFeat$featRatio=annFeat$feature
                            i=match(annFeat$feature,ann$feature); i2=which(!is.na(i)); i1=i[i2]
                            annFeat$featRatio[i2][which(ann$complex[i1]==0)]=NA
                            annFeat$featRatio[is.na(i)]=NA
                            
                            phenAll=annSam
                            
                            if (subsetFFlag!="") {
                                k=sub("ClustOrd","",orderSam)
                                #k=as.integer(substr(k,nchar(k),nchar(k)))
                                k=as.integer(sub("_kmeans|_hclust","",k))
                                p=which(colnames(wssMat)==paste("_",k,sub("_","",sub(k,"",sub("ClustOrd","",orderSam))),sep=""))
                                i2=order(wssMat[,p],decreasing=T)
                                if (subsetFFlag=="_selectedFeat") {
                                    i2=i2[(wssCutoff[p]+1):nrow(wssMat)]
                                } else {
                                    i2=i2[1:wssCutoff[p]]
                                }
                                i=match(rownames(wssMat)[i2],annFeat$feature)
                                arrayData=arrayData[i,]
                                annFeat=annFeat[i,]
                            }
                            
                            #annFeatAll=annFeat
                            #nClustAll=nClust
                            
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
                            
                            if (length(grep("_class",subsetFlag))==1) {
                                #header=paste(header,", no zernike2...",sep="")
                                i=1:nrow(annFeat)
                                j=1:nrow(annSam)
                                if (subsetFlag=="_classSet1") {
                                    j=which(annSam$set1Class!="")
                                } else if (subsetFlag=="_classSet1Set2") {
                                    j=which(annSam$set1Class!="" | annSam$set2Class!="")
                                } else if (subsetFlag=="_classSet1MissClass") {
                                    j=which(annSam$set1Class!=annSam$set1ClassPred)
                                } else if (subsetFlag=="_classSet1MissClassRpart") {
                                    i=match(annFeat$feature,annRpart$feature[which(annRpart$type==typeFlag)]); i=which(!is.na(i))
                                    j=which(annSam$set1Class!=annSam$set1ClassPred)
                                    nClust[1]=NA
                                } else if (subsetFlag=="_classSet1Set2MissClassRpart") {
                                    i=match(annFeat$feature,annRpart$feature[which(annRpart$type==typeFlag)]); i=which(!is.na(i))
                                    j=which((annSam$set1Class!="" & annSam$set1Class!=annSam$set1ClassPred) | (annSam$set2Class!="" & annSam$set2Class!=annSam$set1ClassPred))
                                    nClust[1]=NA
                                }
                                arrayData=arrayData[i,j]
                                annSam=annSam[j,]
                                annFeat=annFeat[i,]
                            }
                            annFeatAll=annFeat
                            
                            varFListAll=varFList
                            varFNameAll=varFName
                            
                            cat("\n\n======================",fNameOut,"======================\n")
                            
                            
                            ## -------------------
                            if (typeFlag=="") arrayData2=arrayData else arrayData2=arrayData[which(annFeat$type==typeFlag),]
                            #save(arrayData2,file=paste("arrayData2",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,".RData",sep=""))
                            #save(arrayData2,file=paste("arrayData2",cohort,subsetFlag,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,".RData",sep=""))
                            
                            if (sampleBar=="cluster") {
                                fNameOut4=sub(orderFlag[2],"",sub(orderFlag[1],"",fNameOut))
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
                                        distMat=dist(t(arrayData2), method=distMethod)
                                    }
                                )
                                clustC=hclust(distMat, method=linkMethod)
                                
                                if (nbClustFlag) {
                                    cat("\n\n-----------nbClust\n")
                                    timeStamp=Sys.time()
                                    print(format(timeStamp, "%x %X"))
                                    set.seed(98765)
                                    res=NbClust(data=t(arrayData2),diss=distMat,distance=NULL,min.nc=2,max.nc=15,method="ward.D2",index="alllong")
                                    timeStamp=c(timeStamp,Sys.time())
                                    print(format(timeStamp[2], "%x %X"))
                                    print(diff(timeStamp))
                                    save(res,file=paste("nbClust",fNameOut4,".RData",sep=""))
                                }
                            } else {
                                if (orderFlag[2]!="") {
                                    #fNameOut4=sub(orderFlag[2],"",sub(orderFlag[1],"",fNameOut))
                                    fNameOut4=sub(orderFlag[1],"",fNameOut)
                                    load(file=paste("kmeans_",gsub("_kmeans|ClustOrd","",orderFlag[2]),"clusts_10000starts.RData",sep=""))
                                    j=match(colnames(arrayData2),names(fit$cluster))
                                    tbl=data.frame(id=colnames(arrayData2),clustId=paste("cluster",fit$cluster[j],sep=""),order=rep(NA,ncol(arrayData2)))
                                    grpUniq=sort(unique(fit$cluster))
                                    j=c()
                                    for (gId in 1:length(grpUniq)) {
                                        jj=match(names(fit$cluster)[which(fit$cluster==grpUniq[gId])],colnames(arrayData2))
                                        switch(distMethod,
                                            "pearson"={
                                                distMat=as.dist(1 - cor(arrayData2[,jj],method=distMethod,use="complete.obs"))
                                            },
                                            "spearman"={
                                                distMat=as.dist(1 - cor(arrayData2[,jj],method=distMethod,use="complete.obs"))
                                            },
                                            "kendall"={
                                                if (getCorFlag) {
                                                    corMatSam2=cor(arrayData2[,jj],method=distMethod,use="complete.obs")
                                                    #save(corMatSam2,file=paste("corMatSam",fNameOut4,".RData",sep=""))
                                                } else {
                                                }
                                                distMat=as.dist(1 - corMatSam2)
                                            },
                                            "euclidean"={
                                                distMat=dist(t(arrayData2[,jj]), method=distMethod)
                                            }
                                        )
                                        if (distMethod=="kendall" & !getCorFlag) {
                                            #clustInfo=read.table(paste(datadir2[2],fNameOut4,"/clusterInfoSample",paste(fNameOut4,"_",type4Flag,sep=""),".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                        } else {
                                            clustC=hclust(distMat, method=linkMethod)
                                            j=c(j,match(clustC$label[clustC$order],tbl$id))
                                        }
                                    }
                                    if (distMethod=="kendall" & !getCorFlag) {
                                        #clustInfo=read.table(paste(datadir2[2],fNameOut4,"/clusterInfoSample",paste(fNameOut4,"_",type4Flag,sep=""),".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                        clustInfo=read.table(paste(datadir2[2],"clusterInfoSample",paste(fNameOut4,sep=""),".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                        j=match(clustInfo$id,tbl$id)
                                    } else {
                                        if (distMethod=="kendall") {
                                            tbl2=tbl[j,]
                                            tbl2$order=1:nrow(tbl2)
                                            write.table(tbl2,paste("clusterInfoSample",fNameOut4,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                                        }
                                    }
                                    annSam$clust=tbl$clustId
                                    annSam=annSam[j,]
                                    arrayData=arrayData[,j]
                                    nm=paste("clust",gsub("_kmeans|ClustOrd","",orderFlag[2]),sep="")
                                    names(annSam)[match("clust",names(annSam))]=nm
                                    varList=c(varList,nm)
                                    varName=c(varName,paste(c("clust")," ",sep=""))
                                    
                                }
                                clustC=NA
                                nClust[2]=NA
                            }
                            #phenAll=annSam
                            arrayDataAll=arrayData
                            clustCAll=clustC
                            nClustAll=nClust
                            varListAll=varList
                            varNameAll=varName
                            
                            
                            if (predictionStrengthFlag) {
                                if (F) {
                                    x=arrayData2
                                    x=arrayData2[1:10,1:4]
                                    y=dist(x)
                                    timeStamp=Sys.time()
                                    print(format(timeStamp, "%x %X"))
                                    set.seed(98765)
                                    res=prediction.strength(,2,3,M=3,clustermethod=hclustCBI,classification=classifFlag,method="ward.D2",scaling=FALSE)
                                    timeStamp=c(timeStamp,Sys.time())
                                    print(format(timeStamp[2], "%x %X"))
                                    print(diff(timeStamp))
                                    res=hclustCBI(data=x,k=3,cut="number",method="ward.D2",scaling=F,noisecut=0)
                                    res=prediction.strength(x,2,3,M=3,clustermethod=hclustCBI,method="ward.D2",scaling=TRUE)
                                    res=prediction.strength(x,2,3,M=3,clustermethod=claraCBI)
                                }
                                
                                timeStamp=Sys.time()
                                print(format(timeStamp, "%x %X"))
                                set.seed(98765)
                                res=prediction.strength(distMat,2,15,M=10,clustermethod=hclustCBI,classification=classifFlag,method="ward.D2",scaling=FALSE)
                                timeStamp=c(timeStamp,Sys.time())
                                print(format(timeStamp[2], "%x %X"))
                                print(diff(timeStamp))
                                save(res,file=paste("predictionStrengthSam_",classifFlag,fNameOut4,".RData",sep=""))
                            }
                            
                            ## -------------------
                            #if (cohort!="_mycRasWt") {
                            if (sepFClustFlag) {
                                type3List=sort(unique(annFeatAll$type))
                            } else {
                                type3List=typeFlag
                            }
                            #par(mfrow=c(length(type3List),1))
                            for (type3Flag in type3List) {
                                #cat("===========",type3Flag,"1 =========\n\n")
                                if (sepFClustFlag & type2Flag%in%c("","_reducedBioFeatPC","_reducedBioFeatMean")) {
                                    fNameOut2=paste(fNameOut,"_",type3Flag,sep="")
                                } else {
                                    fNameOut2=fNameOut
                                }
                                fNameOut4=fNameOut2
                                if (typeFlag=="") i=1:nrow(annFeatAll) else i=which(annFeatAll$type==type3Flag)
                                arrayData=arrayData2=arrayDataAll[i,]
                                annFeat=annFeatAll[i,]
                                nClust=nClustAll
                                if (!is.na(nClust[1]) & nClust[1]>=nrow(annFeat)) nClust[1]=NA
                                if (!is.na(nClust[2]) & nClust[2]>=nrow(annSam)) nClust[2]=NA
                                #if (!sepFClustFlag & type2Flag%in%c("_reducedBioFeatPC","_reducedBioFeatMean")) {
                                if (!sepFClustFlag & typeFlag!="" & type2Flag%in%c("","_reducedBioFeatPC","_reducedBioFeatMean")) {
                                    arrayData=arrayDataAll
                                    annFeat=annFeatAll
                                    x=c()
                                    for (type4Flag in sort(unique(annFeatAll$type))) {
                                        fNameOut4=paste(fNameOut,"_",type4Flag,sep="")
                                        arrayData2=arrayData[which(annFeat$type==type4Flag),]
                                        if (orderFlag[1]%in%c("_mycRas289Ord","_wt906Ord","_mycRasOrd","_wtOrd","_wt_classSet1Ord")) {
                                            #fNameOut4=sub(cohort,"",sub(subsetFFlag,"",sub(subsetFlag,"",sub(orderFlag[2],"",sub(orderFlag[1],"",fNameOut2)))))
                                            fNameOut4=sub(cohort,"",sub(subsetFFlag,"",sub(subsetFlag,"",sub(orderFlag[2],"",sub(orderFlag[1],sub("Ord","",orderFlag[1]),fNameOut2)))))
                                            #fNameOut4=sub("spearman","kendall",sub(cohort,"",sub(subsetFFlag,"",sub(subsetFlag,"",sub(orderFlag[2],"",sub(orderFlag[1],"",fNameOut2))))))
                                            clustInfo=read.table(paste(datadir2[1],fNameOut4,"/clusterInfoFeature",paste(fNameOut4,"_",type4Flag,sep=""),".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                            #clustInfo=read.table(paste(datadir2[1],fNameOut4,"/clusterInfoFeature",paste(fNameOut4,sep=""),".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                            i=match(clustInfo$feature,annFeat$feature)
                                            i=i[which(!is.na(i))]
                                            x=c(x,annFeat$feature[i])
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
                                                    distMat=dist(arrayData2, method=distMethod)
                                                }
                                            )
                                            clustR=hclust(distMat, method=linkMethod)
                                            x=c(x,clustR$labels[clustR$order])
                                        }
                                    }
                                    i=match(x,annFeat$feature)
                                    i=i[which(!is.na(i))]
                                    arrayData=arrayData[i,]
                                    annFeat=annFeat[i,]
                                    clustR=NA
                                    nClust[1]=NA
                                } else {
                                    if (orderFlag[1]%in%c("_mycRas289Ord","_wt906Ord","_mycRasOrd","_wtOrd","_wt_classSet1Ord")) {
                                        #fNameOut4=sub(cohort,"",sub(subsetFFlag,"",sub(subsetFlag,"",sub(orderFlag[2],"",sub(orderFlag[1],"",fNameOut)))))
                                        fNameOut4=sub(cohort,"",sub(subsetFFlag,"",sub(subsetFlag,"",sub(orderFlag[2],"",sub(orderFlag[1],sub("Ord","",orderFlag[1]),fNameOut)))))
                                        clustInfo=read.table(paste(datadir2[1],fNameOut4,"/clusterInfoFeature",fNameOut4,ifelse(type3Flag=="","","_"),type3Flag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                        i=match(clustInfo$feature,annFeatAll$feature)
                                        i=i[which(!is.na(i))]
                                        arrayData=arrayDataAll[i,]
                                        annFeat=annFeatAll[i,]
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
                                        if (sepFClustFlag & type2Flag%in%c("","_reducedBioFeatPC","_reducedBioFeatMean")) {
                                            if (type3Flag==typeFlag) {
                                                clustC=clustCAll
                                                nClust[2]=nClustAll[2]
                                            } else {
                                                i=which(annFeatAll$type==type3Flag)
                                                j=clustCAll$order
                                                arrayData=arrayData2=arrayDataAll[i,j]
                                                annFeat=annFeatAll[i,]
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
                                #} ## for (type3Flag in type3List) {
                                
                                ## -------------------
                                if (nrow(arrayData)>10) {
                                    nameRow=rep("",nrow(arrayData))
                                    i=which(!is.na(annFeat$featRatio))
                                    nameRow[i]=annFeat$featRatio[i]
                                } else {
                                    nameRow=annFeat$feature
                                }
                                nameRow=rep("",nrow(arrayData))
                                if (is.null(varFList)) {
                                    colRow=NULL
                                } else {
                                    colRow=matrix(nrow=length(varFList),ncol=nrow(annFeat))
                                    for (varId in 1:length(varFList)) {
                                        if (varFList[varId]%in%c("sd")) {
                                            j=match(annFeat$feature,annFeatAll$feature)
                                            x=round(100*annFeatAll[,varFList[varId]])+1
                                            lim=range(x,na.rm=T)
                                            x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                                            grpUniq=lim[1]:lim[2]
                                            colRowUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                            colRow[varId,]=colRowUniq[x[j]]
                                        } else {
                                            x=as.character(annFeatAll[,varFList[varId]])
                                            x[x==""]=NA; x=as.integer(as.factor(x))
                                            grpUniq=sort(unique(x))
                                            x=x[match(annFeat$feature,annFeatAll$feature)]
                                            if (length(grpUniq)<=length(colList2)) {
                                                colRow[varId,]=colList2[x]
                                            } else if (length(grpUniq)<=length(colList)) {
                                                colRow[varId,]=colList[x]
                                            } else {
                                                colRow[varId,]=rainbow(length(grpUniq))[x]
                                            }
                                        }
                                    }
                                    rownames(colRow)=varFName
                                }
                                
                                ## -------------------
                                if (ncol(arrayData>10)) {
                                    nameCol=rep("",ncol(arrayData))
                                } else {
                                    nameCol=colnames(arrayData)
                                }
                                colCol=NULL
                                if (!is.null(varList) & (!sepFClustFlag | (sepFClustFlag & type3Flag==typeFlag & cohort=="_mycRasWt"))) {
                                    #if (!sepFClustFlag | (sepFClustFlag & type3Flag==typeFlag)) {
                                    #if (!sepFClustFlag) {
                                    colCol=matrix(nrow=length(varList),ncol=nrow(annSam))
                                    for (varId in 1:length(varList)) {
                                        if (varList[varId]%in%c("sd")) {
                                            #if (varList[varId]%in%c("sd")) {
                                            j=match(annSam$id,phenAll$id)
                                            x=round(100*phenAll[,varList[varId]])+1
                                            lim=range(x,na.rm=T)
                                            if (varList[varId]==c("sd")) lim=limSdSam
                                            if (length(grep("dist2class",varList[varId]))==1) lim=limDist2classSam
                                            x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                                            grpUniq=lim[1]:lim[2]
                                            colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                            colCol[varId,]=colColUniq[x[j]]
                                        } else if (length(grep("dist2class",varList[varId]))==1) {
                                            #if (varList[varId]%in%c("sd")) {
                                            j=match(annSam$id,phenAll$id)
                                            x=round(phenAll[,varList[varId]])
                                            lim=range(x,na.rm=T)
                                            #if (length(grep("dist2class",varList[varId]))==1) lim=limDist2classSam
                                            grpUniq=lim[1]:lim[2]
                                            colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                            colCol[varId,]=colColUniq[x[j]]
                                        } else {
                                            x=as.character(phenAll[,varList[varId]])
                                            x[x==""]=NA; x=as.integer(as.factor(x))
                                            grpUniq=sort(unique(x))
                                            x=x[match(annSam$id,phenAll$id)]
                                            if (substr(varList[varId],1,5)%in%c("clust") | length(grep("Class",varList[varId]))==1) {
                                                colCol[varId,]=colList[x]
                                            } else {
                                                if (length(grpUniq)<=length(colList2)) {
                                                    colCol[varId,]=colList2[x]
                                                } else if (length(grpUniq)<=length(colList)) {
                                                    colCol[varId,]=colList[x]
                                                } else {
                                                    colCol[varId,]=rainbow(length(grpUniq))[x]
                                                }
                                            }
                                        }
                                    }
                                    rownames(colCol)=varName
                                }
                                
                                ## -------------------
                                #print("summary(c(arrayData))")
                                #print(summary(c(arrayData)))
                                #print("quantile(abs(c(arrayData)),probs=seq(0,1,by=.1),na.rm=T)")
                                #print(quantile(abs(c(arrayData)),probs=seq(0,1,by=.1),na.rm=T))
                                main=NULL
                                main=header
                                
                                if (class(clustC)=="hclust") nClust[2]=ifelse(ncol(arrayData)>10,10,NA)
                                
                                #cat("===========",type3Flag,"3 =========\n\n")
                                subDir=""
                                subDir <- paste(fNameOut,sep="")			
                                if (!file.exists(subDir)){
                                    dir.create(file.path(subDir))
                                }
                                subDir=paste(subDir,"/",sep="")
                                if (outFormat=="png") {
                                    if (orderFlag[1]%in%c("_mycRas289Ord","_wt906Ord","_mycRasOrd","_wtOrd","_wt_classSet1Ord") & sepFClustFlag & type3Flag!=typeFlag) {
                                        margins=c(23,36)
                                    } else {
                                        margins=c(6,1)
                                        margins=c(10,20)
                                        margins=c(5,20)
                                    }
                                    png(paste(subDir,"heatmap",fNameOut2,".png",sep=""),width=480*2,height=480*2)
                                } else {
                                    margins=c(12,5)
                                    pdf(paste(subDir,"heatmap",fNameOut2,".pdf",sep=""))
                                }
                                if (class(clustC)=="hclust") nClust[2]=5

                                hcc=heatmap3(x=arrayData, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=colRow, labCol=nameCol, labRow=nameRow, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit2,cexCol=2, , high=colHM[1], low=colHM[2], mid=colHM[3])
                                dev.off()
                                
                                cat("\n\n=========== 101 ",fNameOut,"\nnClust[2] ",nClust[2],"\n",sep="")

                                ## -------------------
                                if (is.na(nClust[1])) {
                                    if (class(clustR)=="hclust") i=clustR$order else i=1:nrow(annFeat)
                                    tbl=cbind(annFeat[i,],clustId="cluster1",order=1:nrow(annFeat))
                                } else {
                                    if (F) {
                                        pdf(paste(subDir,"clusterFeatures",fNameOut2,".pdf",sep=""))
                                        plot(clustR,main=paste("Feature clusters with ",nClust[1]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustR,k=nClust[1])
                                        dev.off()
                                    }
                                    clustId=cutree(clustR,k=nClust[1])[clustR$order]
                                    k1=which(!duplicated(clustId))
                                    for (k in 1:length(k1)) {
                                        clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                                    }
                                    tbl=cbind(annFeat[clustR$order,],clustId,order=1:nrow(annFeat))
                                }
                                write.table(tbl, paste(subDir,"clusterInfoFeature",fNameOut2,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                                
                                if (!sepFClustFlag | (sepFClustFlag & type3Flag==typeFlag)) {
                                    if (is.na(nClust[2])) {
                                        if (class(clustC)=="hclust") j=clustR$order else j=1:nrow(annSam)
                                        tbl=cbind(annSam[j,],clustId="cluster1",order=1:nrow(annSam))
                                    } else {
                                        if (F) {
                                            pdf(paste(subDir,"clusterSamples",fNameOut,".pdf",sep=""))
                                            plot(clustC,main=paste("Sample clusters with ",nClust[2]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustC,k=nClust[2])
                                            dev.off()
                                        }
                                        clustId=cutree(clustC,k=nClust[2])[clustC$order]
                                        k1=which(!duplicated(clustId))
                                        for (k in 1:length(k1)) {
                                            clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                                        }
                                        tbl=cbind(annSam[clustC$order,],clustId,order=1:nrow(annSam))
                                        out=matrix(nrow=nrow(tbl),ncol=9)
                                        colnames(out)=paste("clustId_",2:10,sep="")
                                        for (kk in 1:ncol(out)) {
                                            clustId=cutree(clustC,k=as.integer(sub("clustId_","",colnames(out)[kk])))[clustC$order]
                                            k1=which(!duplicated(clustId))
                                            for (k in 1:length(k1)) {
                                                clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                                            }
                                            out[,kk]=clustId
                                        }
                                        tbl=cbind(tbl,out)
                                    }
                                    write.table(tbl, paste(subDir,"clusterInfoSample",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                                }
                            } ## for (type3Flag in type3List) {
                            
                            if (!is.null(colRow)) {
                                for (varId in 1:length(varFListAll)) {
                                    if (varFListAll[varId]%in%c("sd")) {
                                        width = 480; height = 140
                                    } else {
                                        width = 480; height = 480
                                    }
                                    if (outFormat=="png") {
                                        png(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],".png",sep=""),width=width,height=height)
                                    } else {
                                        pdf(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],".pdf",sep=""))
                                    }
                                    if (varFListAll[varId]%in%c("sd")) {
                                        x=round(100*annFeatAll[,varFListAll[varId]])+1
                                        lim=range(x,na.rm=T)
                                        grpUniq=lim[1]:lim[2]
                                        colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                        lim=(lim-1)/100
                                        #heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1)],median(1:length(colColUniq))))
                                        heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1,median(1:length(colColUniq)))]))
                                    } else {
                                        x=as.character(annFeatAll[,varFListAll[varId]]); x[x==""]=NA
                                        grpUniq=table(x)
                                        grpUniq=paste(names(grpUniq)," (",grpUniq,")",sep="")
                                        k=1:length(grpUniq)
                                        if (length(grpUniq)<=length(colList2)) {
                                            sampleColorLegend(tls=grpUniq[k],col=colList2,legendTitle=varFNameAll[varId])
                                        } else if (length(grpUniq)<=length(colList)) {
                                            cexThis=NULL
                                            if (length(grpUniq)>10) cexThis=1
                                            sampleColorLegend(tls=grpUniq[k],col=colList,legendTitle=varFNameAll[varId],cex=cexThis)
                                        } else {
                                            sampleColorLegend(tls=grpUniq[k],col=rainbow(length(grpUniq)),legendTitle=varFNameAll[varId])
                                        }
                                    }
                                    dev.off()
                                }
                            }
                            
                            if (!is.null(colCol)) {
                                for (varId in 1:length(varListAll)) {
                                    if (varListAll[varId]%in%c("sd")) {
                                        width=480; height=140
                                    } else {
                                        width=480; height=480
                                    }
                                    if (outFormat=="png") {
                                        png(paste("heatmapSampleColorBarLegend_",varListAll[varId],".png",sep=""),width=width,height=height)
                                    } else {
                                        pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],".pdf",sep=""))
                                    }
                                    if (varList[varId]%in%c("sd")) {
                                        #if (varListAll[varId]%in%c("sd")) {
                                        x=round(100*phenAll[,varListAll[varId]])+1
                                        lim=range(x,na.rm=T)
                                        if (varListAll[varId]==c("sd")) lim=limSdSam
                                        if (length(grep("dist2class",varList[varId]))==1) lim=100*limDist2classSam
                                        grpUniq=lim[1]:lim[2]
                                        colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                        lim=(lim-1)/100
                                        heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1,median(1:length(colColUniq)))]))
                                    } else if (length(grep("dist2class",varList[varId]))==1) {
                                        #if (varListAll[varId]%in%c("sd")) {
                                        x=round(phenAll[,varListAll[varId]])+1
                                        lim=range(x,na.rm=T)
                                        if (length(grep("dist2class",varList[varId]))==1) lim=limDist2classSam
                                        grpUniq=lim[1]:lim[2]
                                        colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                        lim=lim
                                        heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1,median(1:length(colColUniq)))]))
                                    } else {
                                        x=as.character(phenAll[,varListAll[varId]]); x[x==""]=NA
                                        grpUniq=table(x)
                                        grpUniq=paste(names(grpUniq)," (",grpUniq,")",sep="")
                                        k=1:length(grpUniq)
                                        if (substr(varListAll[varId],1,5)%in%c("clust") | length(grep("Class",varListAll[varId]))==1) {
                                            sampleColorLegend(tls=grpUniq[k],col=colList,legendTitle=varNameAll[varId])
                                        } else {
                                            if (length(grpUniq)<=length(colList2)) {
                                                sampleColorLegend(tls=grpUniq[k],col=colList2,legendTitle=varNameAll[varId])
                                            } else if (length(grpUniq)<=length(colList)) {
                                                sampleColorLegend(tls=grpUniq[k],col=colList,legendTitle=varNameAll[varId])
                                            } else {
                                                sampleColorLegend(tls=grpUniq[k],col=rainbow(length(grpUniq)),legendTitle=varNameAll[varId])
                                            }
                                        }
                                    }
                                    dev.off()
                                }
                            }
                        }
                        
                        if (outFormat=="png") {
                            png(paste("heatmapColorRange",fNameOut,".png",sep=""),width=480,height=140)
                        } else {
                            pdf(paste("heatmapColorRange",fNameOut,".pdf",sep=""))
                        }
                        heatmapColorBar=function(limit,cols=c("green","red","black")) {
                            try <- maPalette(high=cols[1], low=cols[2], mid=cols[3],k=15)
                            maColorBar(try, scale=limit,k=3)
                        }
                        heatmapColorBar(limit=limit1,cols=colHM)
                        dev.off()
                    }
                }
            }
        }
    }
}

####################################################################
####################################################################
## nbClust

colList=c("brown","red","orange","green","cyan","blue","magenta","purple","darkgreen")

png("nbClust.png",width=2*480, height=1*480)
par(mar=c(5, 4, 4, 2) + 0.1)
par(mar=c(7, 4, 4, 2) + 0.1)
datadir="/Users/royr/UCSF/WallaceMarshall/results/nbClust/"
datVerList=c("m289","w906","w906m289")
kk=1
grpUniq=c()
for (datVer in datVerList) {
    switch(datVer,
    "w5969"={
        cohortList=c("_wt","_mycRas","_mycRasWt")
        cohortList=c("_wt","_mycRas")
        cohortList=c("_mycRas")
    },
    "w906"={
        cohortList=c("_wt906")
    },
    "m289"={
        cohortList=c("_mycRas289")
    },
    "w906m289"={
        cohortList=c("_wt906mycRas289")
    }
    )
    typeList=c("",sort(unique(ann$type)))
    offset=0
    for (cohort in cohortList) {
        for (typeFlag in typeList) {
            cat("\n\n=========== ",cohort,", ",typeFlag,"\n",sep="")
            load(paste(datadir,"nbClust",cohort,ifelse(typeFlag=="","","_"),typeFlag,"_reducedBioFeatPC_kendall.RData",sep=""))
            x=apply(res$All.index,2,function(x) {
                nClust=(1:length(x))+1
                k=which(is.finite(x) & !is.na(x))
                y=nClust[k][which.max(x[k])]
                #print(y)
                if (length(y)==0) y=NA
                y
            })
            #print(table(x))
            #print(colnames(res$Best.nc)[which(res$Best.nc["Number_clusters",]==5)])
            #res[["All.index"]][,"Gap"]
            print(res$Best.nc["Number_clusters","SDbw"])
            if (kk==1) {
                plot((1:ncol(res$Best.nc))+offset, res$Best.nc["Number_clusters",],ylim=c(0,16),xlab="",ylab="Best number of clusters",xaxt="n",col=colList[kk])
                axis(side=1,at=1:ncol(res$Best.nc),labels=colnames(res$Best.nc),las=3)
            } else {
                points((1:ncol(res$Best.nc))+offset, res$Best.nc["Number_clusters",],xaxt="n",col=colList[kk])
            }
            grpUniq=c(grpUniq,paste(sub("_","",cohort),typeFlag))
            kk=kk+1
            offset=offset+0.05
        }
    }
}
dev.off()

png("nbClust_ColorBarLegend.png")
sampleColorLegend(tls=grpUniq,col=colList,legendTitle="Best number of clusters from")
dev.off()

####################################################################
####################################################################
## MDS plot

source(paste(dirSrc,"functions/heatmapAcgh.7.4.R",sep=""))

cohort="_wt906"
typeList=c("",sort(unique(ann$type)))
typeFlag="nucleus"
ordFlag="wtOrd"
distMethod="kendall"

nClust=5

switch(cohort,
    "_wt906"={
        cohortName="Wildtype 906"
    },
    "_mycRas289"={
        cohortName="Myc/Ras 289"
    }
)
colId=paste("clustId_",nClust,sep="")
for (typeFlag in typeList) {
    datadir="results/wt906mycRas289/"
    datadir="results/wt906mycRas289_noCorr/"
    load(paste(datadir,"corMatSam",cohort,ifelse(typeFlag=="","","_"),typeFlag,"_",distMethod,".RData",sep=""))
    datadir=paste(datadir,"heatmap/",sub("_","",cohort),"/",ordFlag,"/",cohort,ifelse(typeFlag=="","","_"),typeFlag,"_",distMethod,"/",sep="")
    clustInfo=read.table(paste(datadir,"clusterInfoSample",cohort,ifelse(typeFlag=="","","_"),typeFlag,"_",distMethod,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    clustInfo=clustInfo[match(colnames(corMatSam2),clustInfo$id),]
    header=paste(cohortName,", ",ifelse(typeFlag=="","all features",typeFlag),sep="")
    colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
    clustInfo$col=rep("black",ncol(corMatSam2))
    clustInfo$col=colList[match(clustInfo[,colId],paste("cluster",1:length(colList),sep=""))]
    dat=as.dist(1-abs(corMatSam2))
    mds1=cmdscale(dat, k = 2)
    clustInfo=clustInfo[match(rownames(mds1),clustInfo$id),]

    png(paste("mdsPlot_corMatSam",cohort,typeFlag,"_",distMethod,".png",sep=""))
    plot(mds1[,1], mds1[,2], type = "p", xlab = "MDS1", ylab = "MDS2", axes = T, main = paste(header,": MDS plot of samples",sep=""),col=clustInfo$col)
    #plot(mds1[,1], mds1[,2], type = "n", xlab = "", ylab = "", axes = FALSE, main = "cmdscale (stats)")
    #text(mds1[,1], mds1[,2], labels(dat), cex = 0.9, xpd = TRUE)
    dev.off()
}

grpUniq=paste("cluster",1:nClust,sep="")
png(paste("mdsPlotSampleColorBarLegend_heatmapClusters.png",sep=""))
sampleColorLegend(tls=grpUniq,col=colList[1:length(grpUniq)],legendTitle="Heatmap clusters")
dev.off()

####################################################################
####################################################################
## Compare wt & myc/ras proportions in clusters

cohort="_wt906mycRas289"

distMethod="kendall"
typeFlag=""

datadir=paste("results/wt906mycRas289_noCorr/heatmap/centerByType/",cohort,ifelse(typeFlag=="","","_"),typeFlag,"_centerByType_scaleByType_",distMethod,"/",sep="")
clustInfo=read.table(paste(datadir,"clusterInfoSample",cohort,ifelse(typeFlag=="","","_"),typeFlag,"_centerByType_scaleByType_",distMethod,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
x=as.matrix(table(clustInfo$clustId_5,clustInfo$type))
cbind(round(x[,1]/sum(x[,1]),2),round(x[,2]/sum(x[,2]),2))

datadir=paste("results/wt906mycRas289_noCorr/heatmap/centerTogether/",cohort,ifelse(typeFlag=="","","_"),typeFlag,"_",distMethod,"/",sep="")
clustInfo=read.table(paste(datadir,"clusterInfoSample",cohort,ifelse(typeFlag=="","","_"),typeFlag,"_",distMethod,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
x=as.matrix(table(clustInfo$clustId_5,clustInfo$type))
cbind(round(x[,1]/sum(x[,1]),2),round(x[,2]/sum(x[,2]),2))

####################################################################
####################################################################
