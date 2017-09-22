library("Rtsne")


## -------------------------------
cohort="_mycRas289"
cohort="_wt906"

type2Flag=""
type2Flag="_reducedBioFeatPC"

colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
colListSam=c("green","cyan","red","blue","purple","magenta","purple","darkgreen","brown","red","orange","yellow")
centrFlag=""
scaleFlag=""

for (cohort in c("_mycRas289","_wt906")) {
    switch(cohort,
        "_wt906"={
            cohortName="WT 906"
            ordFlag="_mycRas289Ord"
        },
        "_mycRas289"={
            cohortName="Myc/Ras 289"
            ordFlag=""
        }
    )
    if (type2Flag=="_reducedBioFeatPC") {
        corFlag=F
        switch(cohort,
        "_wt906"={
            cell=cell_rbfW2
            ann=ann_rbf2
            annCell=annCellW2
        },
        "_mycRas289"={
            cell=cell_rbfM2
            ann=ann_rbf2
            annCell=annCellM2
        }
        )
    } else {
        corFlag=T
        switch(cohort,
            "_wt906"={
                cell=cellW2
                ann=annW2
                annCell=annCellW2
            },
            "_mycRas289"={
                cell=cellM2
                ann=annM2
                annCell=annCellM2
            }
        )
    }

    typeList=c("",sort(unique(ann$type)))

    corThres=.8
    sdThres=.05
    sdThres=NA
    datTypeList=c("_data","_corKend")
    datTypeList=c("_data")
    
    x=apply(cell,2,function(x) {mean(!is.na(x))})
    cell=apply(cell,2,function(x) {
        y=x
        y0=mean(y,na.rm=T)
        y[is.na(y)]=y0
        y
    })
    featId=1:10; samId=1:20
    featId=1:ncol(cell); samId=1:nrow(cell)
    if (corFlag) {
        datadir=""
        corInfo=read.table(paste(datadir,"corrFeature",cohort,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        featId=which(!ann$feature%in%corInfo$feature2[which(corInfo$corKend>=corThres)]); samId=1:nrow(cell)
    }
    cell=cell[samId,featId]
    ann=ann[featId,]
    annCell=annCell[samId,]
    arrayData=cell
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
    colVec=rep("black",nrow(ann))
    grpUniq=sort(unique(ann$type))
    for (gId in 1:length(grpUniq)) {
        colVec[which(ann$type==grpUniq[gId])]=colList[gId]
    }
    datadir=paste("results/",sub("_","",cohort),"/",sep="")
    for (datType in datTypeList) {
        #for (dat2Type in c("_sample","_feature")) {
        for (dat2Type in c("_sample")) {
            for (typeFlag in typeList) {
                fNameOut=paste(cohort,orderFlag,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),sep="")
                datadirG=paste(cohort,ordFlag,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,"_kendall/",sep="")
                fNameG=paste("clusterInfoSample",cohort,ordFlag,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,"_kendall.txt",sep="")
                
                clustInfo=read.table(paste(datadirG,fNameG,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                j=match(annCell$id,clustInfo$id); j1=which(!is.na(j)); j2=j[j1]
                annCell$clustId=""
                annCell$clustId[j1]=clustInfo$clustId_5[j2]
                
                type3List=c("",sort(unique(annFeatAll$type)))
                for (type3Flag in type3List) {
                    #fName=paste(cohort,ifelse(typeFlag=="","","_"),typeFlag,dat2Type,datType,sep="")
                    fName=paste(fNameOut,ifelse(type3Flag=="","",paste("_",type3Flag,sep="")),dat2Type,datType,sep="")
                    header=paste(cohortName,ifelse(type3Flag=="","",", "),type3Flag,": ",sep="")
                    if (corFlag) {
                        i=match(colnames(arrayData),ann$feature[!ann$feature%in%corInfo$feature2[which(abs(corInfo$corKend)>=corThres)]]); i=which(!is.na(i))
                    } else {
                        i=1:nrow(ann)
                    }
                    if (type3Flag=="") {
                        #i=1:nrow(ann)
                    } else {
                        i=i[which(ann$type[i]==type3Flag)]
                    }
                    if (dat2Type=="_sample") {
                        header=paste(header,"Cell",sep="")
                        ttl="Cell"
                        col2Vec=rep("black",nrow(ann))
                    } else {
                        header=paste(header,"Feature",sep="")
                        ttl="Feature"
                        col2Vec=colVec
                    }
                    switch(datType,
                        "_data"={
                            header=paste(header,", data based",sep="")
                            if (dat2Type=="_sample") {
                                dat=arrayData[,i]
                                col2Vec=colListSam[as.integer(sub("cluster","",annCell$clustId))]
                            } else {
                                dat=t(arrayData[,i])
                                col2Vec=rep("black",nrow(dat))
                                i=match(rownames(dat),ann$feature)
                                col2Vec=colVec[i]
                            }
                        },
                        "_corKend"={
                            header=paste(header,", Kendall's correlation based",sep="")
                            if (dat2Type=="_sample") {
                                load(paste(datadir,"corMatSam",cohort,"_kendall.RData",sep=""))
                                dat=corMatSam2[samId,samId]
                            } else {
                                load(paste(datadir,"corMat",cohort,"_kendall.RData",sep=""))
                                dat=corMat2[featId,featId][i,i]
                                i=match(ann$feature,rownames(dat)); i=i[which(!is.na(i))]
                                dat=dat[i,i]
                                col2Vec=rep("black",nrow(dat))
                                i=match(rownames(dat),ann$feature)
                                col2Vec=colVec[i]
                            }
                        }
                    )
                    header=paste(header,"\n(sample color from hclust clusters based on ",ifelse(typeFlag=="","all",typeFlag)," features)",sep="")
                    if (!is.na(sdThres)) {
                        i=apply(dat,1,function(x) {sd(x,na.rm=T)>sdThres})
                        j=apply(dat,2,function(x) {sd(x,na.rm=T)>sdThres})
                        dat=dat[i,j]
                        col2Vec=col2Vec[i]
                    }
                    subDir=""
                    subDir=paste("tsne",fNameOut,sep="")
                    if (!file.exists(subDir)){
                        dir.create(file.path(subDir))
                    }
                    subDir=paste(subDir,"/",sep="")
                    ttl=paste(ttl,": Dim ",1:2,sep="")
                    set.seed(42) # Set a seed if you want reproducible results
                    perplexFlag=30
                    fit=try(Rtsne(dat,perplexity=perplexFlag))
                    if (class(fit)!="try-error") {
                        png(paste(subDir,"plotTsne_perplex",perplexFlag,fName,".png",sep=""))
                        plot(fit$Y,col=col2Vec,main=paste(header,"\nperplexity ",perplexFlag,sep=""),xlab=ttl[1],ylab=ttl[2],pch=20)
                        dev.off()
                    }
                    perplexFlag=5
                    fit=try(Rtsne(dat,perplexity=perplexFlag))
                    if (class(fit)!="try-error") {
                        png(paste(subDir,"plotTsne_perplex",perplexFlag,fName,".png",sep=""))
                        plot(fit$Y,col=col2Vec,main=paste(header,"\nperplexity ",perplexFlag,sep=""),xlab=ttl[1],ylab=ttl[2],pch=20)
                        dev.off()
                    }
                }
                png(paste("tsneSampleColorLegend_heatmapCluster",fNameOut,".png",sep=""))
                x=as.character(annCell[,"clustId"]); x[x==""]=NA
                grpUniq=table(x)
                grpUniq=paste(names(grpUniq)," (",grpUniq,")",sep="")
                k=1:length(grpUniq)
                sampleColorLegend(tls=grpUniq[k],col=colListSam,legendTitle=paste(cohortName,": Heatmap clusters",sep=""),cex=1.5)
                dev.off()
            }
        }
    }
}

## -------------------------------
## NOT USED

if (F) {
    data(iris)
    ris_unique <- unique(iris) # Remove duplicates
    iris_matrix <- as.matrix(iris_unique[,1:4])
    fit <- Rtsne(iris_matrix) # Run TSNE
    plot(fit$Y,col=iris_unique$Species)

    # Using a dist object
    fit <- Rtsne(dist(iris_matrix))
    plot(fit$Y,col=iris_unique$Species)

    # Use a given initialization of the locations of the points
    tsne_part1 <- Rtsne(iris_unique[,1:4], theta=0.0, pca=FALSE,max_iter=350)
    tsne_part2 <- Rtsne(iris_unique[,1:4], theta=0.0, pca=FALSE, max_iter=150,Y_init=tsne_part1$Y)
}

## -------------------------------
