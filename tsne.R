library("Rtsne")


## -------------------------------
cohort="_wt906"
cohort="_mycRas105"

colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
centrFlag=""
scaleFlag=""

switch(cohort,
    "_wt906"={
        cell=cellW2
        ann=annW2
        annCell=annCellW2
    },
    "_mycRas105"={
        cell=cellM2
        ann=annM2
        annCell=annCellM2
    }
)
typeList=c("",sort(unique(ann$type)))

corThres=.8
sdThres=.05

x=apply(cell,2,function(x) {mean(!is.na(x))})
cell=apply(cell,2,function(x) {
    y=x
    y0=mean(y,na.rm=T)
    y[is.na(y)]=y0
    y
})
datadir=""
corInfo=read.table(paste(datadir,"corrFeature",cohort,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
featId=1:10; samId=1:20
featId=1:ncol(cell); samId=1:nrow(cell)
featId=which(!ann$feature%in%corInfo$feature2[which(corInfo$corKend>=corThres)]); samId=1:nrow(cell)
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
for (datType in c("_data","_corKend")) {
    for (dat2Type in c("_sample","_feature")) {
        for (typeFlag in typeList) {
            fName=paste(ifelse(typeFlag=="","","_"),typeFlag,dat2Type,datType,sep="")
            header=paste("WT 906",ifelse(typeFlag=="","",", "),typeFlag,": ",sep="")
            i=match(colnames(arrayData),ann$feature[!ann$feature%in%corInfo$feature2[which(abs(corInfo$corKend)>=corThres)]]); i=which(!is.na(i))
            if (typeFlag=="") {
                #i=1:nrow(ann)
            } else {
                i=i[which(ann$type[i]==typeFlag)]
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
            i=apply(dat,1,function(x) {sd(x,na.rm=T)>sdThres})
            j=apply(dat,2,function(x) {sd(x,na.rm=T)>sdThres})
            dat=dat[i,j]
            col2Vec=col2Vec[i]
            ttl=paste(ttl,": Dim ",1:2,sep="")
            set.seed(42) # Set a seed if you want reproducible results
            perplexFlag=30
            fit=try(Rtsne(dat,perplexity=perplexFlag))
            if (class(fit)!="try-error") {
                png(paste("plotTsne_perplex",perplexFlag,fName,".png",sep=""))
                plot(fit$Y,col=col2Vec,main=paste(header,"\nperplexity ",perplexFlag,sep=""),xlab=ttl[1],ylab=ttl[2],pch=20)
                dev.off()
            }
            perplexFlag=5
            fit=try(Rtsne(dat,perplexity=perplexFlag))
            if (class(fit)!="try-error") {
                png(paste("plotTsne_perplex",perplexFlag,fName,".png",sep=""))
                plot(fit$Y,col=col2Vec,main=paste(header,"\nperplexity ",perplexFlag,sep=""),xlab=ttl[1],ylab=ttl[2],pch=20)
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
