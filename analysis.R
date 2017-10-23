## Run sections 1 and 2

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
capWords=function(s, strict=FALSE) {
    cap=function(s) paste(toupper(substring(s, 1, 1)),
    {s=substring(s, 2); if(strict) tolower(s) else s},
    sep="", collapse=" " )
    sapply(strsplit(s, split=" "), cap, USE.NAMES=!is.null(names(s)))
}

####################################################################
####################################################################
## Section 1
## ----------------------------------------------
datadir="docs/"
cellW2=read.table(paste(datadir,"20170707_WT_RowAB_CONCATENATED.csv",sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T)
annW2=read.table(paste(datadir,"20170729_cellFeatures_categorized.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
for (k in 1:ncol(annW2)) if (is.character(annW2[,k])) annW2[,k]=gsub("'","",annW2[,k])
out=as.data.frame(t(sapply(annW2[which(annW2[,11]!=""),11],function(x) {
    y=strsplit(x," = ")[[1]]
    gsub(" |:","",y)
},USE.NAMES=F)),stringsAsFactors=F)
names(out)=c("featId","featUniq")
annW2=annW2[,c(1,7)]
names(annW2)=c("feature","featId")
annW2$featUniq=out$featUniq[match(annW2$featId,out$featId)]
rm(out)
annW2$type=""
if (F) {
    annW2=read.table(paste(datadir,"cellFeatures.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
    tmpC=rep("",ncol(cellW2))
    annW2=data.frame(feature=annW2[,1],type=tmpC,stringsAsFactors=F)
}
annW2$type[which(substr(tolower(annW2$feature),1,nchar("cell"))=="cell")]="cell"
annW2$type[which(substr(tolower(annW2$feature),1,nchar("mito"))=="mito")]="mitochondria"
annW2$type[which(substr(tolower(annW2$feature),1,nchar("nuc"))=="nuc")]="nucleus"
annW2$type[which(substr(tolower(annW2$feature),1,nchar("CelVarEnergy"))==tolower("CelVarEnergy"))]="cell"
annCellW2=data.frame(id=paste("w2_",1:nrow(cellW2),sep=""),stringsAsFactors=F)
colnames(cellW2)=annW2$feature
rownames(cellW2)=annCellW2$id
annCellW2$junk=0
cellW2=as.matrix(cellW2)

## ----------------------------------------------
datadir="docs/"
cellM2=read.table(paste(datadir,"MYCRAS_20150522_RowA1-6_CONCATENATED.csv",sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T)
tbl=read.table(paste(datadir,"MYCRAS_20150522_RowA7-18_CONCATENATED.csv",sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T)
cellM2=rbind(cellM2,tbl)
rm(tbl)
annM2=read.table(paste(datadir,"20170729_cellFeatures_categorized.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
for (k in 1:ncol(annM2)) if (is.character(annM2[,k])) annM2[,k]=gsub("'","",annM2[,k])
out=as.data.frame(t(sapply(annM2[which(annM2[,11]!=""),11],function(x) {
    y=strsplit(x," = ")[[1]]
    gsub(" |:","",y)
},USE.NAMES=F)),stringsAsFactors=F)
names(out)=c("featId","featUniq")
annM2=annM2[,c(1,7)]
names(annM2)=c("feature","featId")
annM2$featUniq=out$featUniq[match(annM2$featId,out$featId)]
rm(out)
annM2$type=""
if (F) {
    annM2=read.table(paste(datadir,"cellFeatures.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
    tmpC=rep("",ncol(cellM2))
    annM2=data.frame(feature=annM2[,1],type=tmpC,stringsAsFactors=F)
}
annM2$type[which(substr(tolower(annM2$feature),1,nchar("cell"))=="cell")]="cell"
annM2$type[which(substr(tolower(annM2$feature),1,nchar("mito"))=="mito")]="mitochondria"
annM2$type[which(substr(tolower(annM2$feature),1,nchar("nuc"))=="nuc")]="nucleus"
annM2$type[which(substr(tolower(annM2$feature),1,nchar("CelVarEnergy"))==tolower("CelVarEnergy"))]="cell"
annCellM2=data.frame(id=paste("m2_",1:nrow(cellM2),sep=""),stringsAsFactors=F)
colnames(cellM2)=annM2$feature
rownames(cellM2)=annCellM2$id
annCellM2$junk=0
cellM2=as.matrix(cellM2)

## -------------------
load("results/tmp_w5969.RData")
annW=ann
ann_rbf1=ann_rbf
ann_rbfm1=ann_rbfm

####################################################################
####################################################################
## Section

## -------------------
## Generate tmp_w5969.RData

datadir="docs/"
cellW=read.table(paste(datadir,"WT_Plate1WellA1-K16_MeasurementsNUMBERSONLY.csv",sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T)

annCell=read.table(paste(datadir,"20160222_CellClassIdentification.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
annCell=data.frame(id=paste("w",annCell$cellNumber,sep=""),class=annCell$classManuallyCurated,stringsAsFactors=F)
tbl=read.table(paste(datadir,"20160523_CellClassIdentification_PartTwo.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
names(tbl)[match(c("Cell.Number..in.A1.K16.dataset.","Class..Manually.Curated."),names(tbl))]=c("id","class")
tbl$id=paste("w",tbl$id,sep="")
tbl$class=tbl$class
annCell$set=1
tbl$set=2
annCell=rbind(annCell,tbl)
tmp=rep(NA,nrow(cellW))
annCellW=data.frame(id=paste("w",1:nrow(cellW),sep=""),class=tmp,set=tmp,stringsAsFactors=F)
j=match(annCellW$id,annCell$id); j1=which(!is.na(j)); j2=j[j1]
annCellW[j1,]=annCell[j2,]

tbl=read.table(paste(datadir,"WT_Plate1WellK17-P24_MeasurementsNUMBERSONLY.csv",sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T)
tmp=annCellW[1:nrow(tbl),]
for (k in 1:ncol(tmp)) {
    tmp[,k]=NA
}
tmp$id=paste("w",(1:nrow(tmp))+nrow(annCellW),sep="")
tmp$set=3
cellW=rbind(cellW,tbl)
annCellW=rbind(annCellW,tmp)

datadir="docs/myc_ras_data/"
cellM=read.table(paste(datadir,"MYCRASPlate1WellRowsA-C_MeasurementsNUMBERSONLY.csv",sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T)

ann=read.table(paste(datadir,"FeaturesList.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=1)
names(ann)=c("feature","colId")
ann$colId=as.integer(ann$colId)
ann=ann[!is.na(ann$colId),]
grep("CellZernike",ann$feature)
k=36:130
ann$feature[k]=sub("zervike","CellZernike",ann$feature[k])
grep("NucZernike",ann$feature)
k=201:297
ann$feature[k]=sub("zervike","NucZernike",ann$feature[k])
ann$type=""
ann$type[which(substr(ann$feature,1,4)=="Cell")]="cell"
ann$type[which(substr(ann$feature,1,3)=="Nuc")]="nucleus"
ann$type[which(substr(ann$feature,1,4)=="Mito")]="mitochondria"

## ----------------------
cell=cellM
cell=as.matrix(cell)
colnames(cell)[ann$colId]=ann$feature
rownames(cell)=paste("m",1:nrow(cell),sep="")
cellM=cell

cell=cellW
cell=as.matrix(cell)
colnames(cell)[ann$colId]=ann$feature
rownames(cell)=paste("w",1:nrow(cell),sep="")
cellW=cell

ann=ann[,which(names(ann)!="colId")]

## ----------------------

for (cohort in c("_mycRas","_wt","_mycRasWt")) {
    switch(cohort,
    "_mycRas"={
        cell=cellM
    },
    "_wt"={
        cell=cellW
    },
    "_mycRasWt"={
        cell=rbind(cellM,cellW)
        #rownames(cell)=c(paste("m",rownames(cellM),sep=""),paste("w",rownames(cellW),sep=""))
    }
    )
    ann2=t(sapply(c("CellAreaByMitoSumArea",
    "CellAreaByNucArea",
    "MitoSumAreaByNucArea",
    "NucAreaByNucNumber",
    "MitoFusedCountByMitoFragmentedCount",
    "MitoFusedSumAreaByMitoFragmentedSumArea",
    "CellContrastByMitoContrast",
    "CellCorrelationByMitoCorrelation",
    "CellEnergyByMitoEnergy",
    "CellHomogeneityByMitoHomogeneity",
    "MitoFusedSumAreaByMitoSumArea",
    "MitoFragmentedSumAreaByMitoSumArea"),function(x) {
        y=strsplit(x,"By")[[1]]
        c(x,y)
    },USE.NAMES=F))
    colnames(ann2)=c("newVar","numVar","dinVar")
    tmp=matrix(nrow=nrow(cell),ncol=nrow(ann2))
    colnames(tmp)=ann2[,"newVar"]
	for (k in 1:nrow(ann2)) {
		j1=which(ann$feature==ann2[k,"numVar"])
		j2=which(ann$feature==ann2[k,"dinVar"])
		tmp[,k]=cell[,j1]/cell[,j2]
		if (any(cell[,j2]==0,na.rm=T)) {
			i=which(cell[,j2]==0)
			tmp[i,k]=NA
			y=cell[,j2]
			y1=cell[,j1]/cell[,j2]
			y2=cell[,j1]/(cell[,j2]+min(cell[which(cell[,j2]>0),j2]))
			png(paste("tmp",cohort,"_",ann$feature[j2],".png",sep=""))
			y=cell[,j1]/cell[,j2]
			i=order(y)
			plot(y)
			y=cell[,j1]/(cell[,j2]+min(cell[which(cell[,j2]>0),j2]))
			points(y[i],col="red")
			dev.off()
			
			if (F) {
			png(paste("tmp",cohort,"_",ann$feature[j2],".png",sep=""))
			plot(y1,y2)
			points(y1[i],y2[i],col="red")
			dev.off()
			}
			
		}
	}

	if (F) {
		tmp=cbind(cell[,"CellArea"]/cell[,"MitoSumArea"],
		cell[,"CellArea"]/cell[,"NucArea"],
		cell[,"MitoSumArea"]/cell[,"NucArea"],
		cell[,"NucArea"]/cell[,"NucNumber"],
		cell[,"MitoFusedCount"]/cell[,"MitoFragmentedCount"],
		cell[,"MitoFusedSumArea"]/cell[,"MitoFragmentedSumArea"],
		cell[,"CellContrast"]/cell[,"MitoContrast"],
		cell[,"CellCorrelation"]/cell[,"MitoCorrelation"],
		cell[,"CellEnergy"]/cell[,"MitoEnergy"],
		cell[,"CellHomogeneity"]/cell[,"MitoHomogeneity"],
		cell[,"MitoFusedSumArea"]/cell[,"MitoSumArea"],
		cell[,"MitoFragmentedSumArea"]/cell[,"MitoSumArea"])
		colnames(tmp)=c("CellAreaByMitoSumArea",
		"CellAreaByNucArea",
		"MitoSumAreaByNucArea",
		"NucAreaByNucNumber",
		"MitoFusedCountByMitoFragmentedCount",
		"MitoFusedSumAreaByMitoFragmentedSumArea",
		"CellContrastByMitoContrast",
		"CellCorrelationByMitoCorrelation",
		"CellEnergyByMitoEnergy",
		"CellHomogeneityByMitoHomogeneity",
		"MitoFusedSumAreaByMitoSumArea",
		"MitoFragmentedSumAreaByMitoSumArea")
	}

	cell=cbind(cell,tmp)
	switch(cohort,
		   "_mycRas"={
		   cellM=cell
		   },
		   "_wt"={
		   cellW=cell
           },
           "_mycRasWt"={
            cellMW=cell
		   }
		   )
	
}

tmp2=ann[match(ann2[,"numVar"],ann$feature),]
tmp2$feature=colnames(tmp)
ann=rbind(ann,tmp2)
ann$complex=0
ann$complex[which(ann$feature%in%ann2[,"newVar"])]=1


## -------------------

x=ann$feature[grep("2",ann$feature)]
x=x[which(substr(x,nchar(x),nchar(x))==2)]
x=x[which(is.na(as.integer(substr(x,nchar(x)-1,nchar(x)-1))))]
x=substr(x,1,nchar(x)-1)
ann$featUniq=ann$feature
for (k in 1:length(x)) {
	cat("\n\n------------------",x[k],"\n")
	i=which(substr(ann$feature,1,nchar(x[k]))==x[k])
	ann$featUniq[i]=x[k]
	y=max(diff(i))
	if (y!=1) cat("Should be 1: ",y,"\n",sep="")
}

annW=ann

####################################################################
####################################################################
## Pair-wise correlation of samples and features

cohort="_wt906"
cohort="_mycRas289"

distMethod="pearson"
distMethod="kendall"

switch(cohort,
    "_mycRas"={
        cohortName="Myc/Ras"
        cell=cellM
    },
    "_wt"={
        cohortName="Wildtype"
        cell=cellW
    },
    "_mycRasWt"={
        cohortName="Myc/Ras & Wildtype"
        cell=cellMW
    },
    "_wt906"={
        cohortName="Wildtype 906"
        cell=cellW2
    },
    "_mycRas289"={
        cohortName="Myc/Ras 289"
        cell=cellM2
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
}
#cell=cell[1:10,]

if (F) {
	out=t(apply(cell,2,function(x) {
		y=c(summary(x))
		res=vector(mode="numeric",length=7)
		names(res)=c( "Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","NA's")
		res[match(names(y),names(res))]=y
		res
	}))
	sumInfo=cbind(ann,out)
	write.table(sumInfo,file=paste("summaryFeature",cohort,".txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)
}

if (F) {
    timeStamp=Sys.time()
	corMatSam=cor(t(cell),use="complete.obs",method=distMethod)
    timeStamp=c(timeStamp,Sys.time())
    print(diff(timeStamp))
	save(corMatSam,file="corMatSam.RData")
	
    timeStamp=Sys.time()
	corMat=cor(cell,use="complete.obs",method=distMethod)
    timeStamp=c(timeStamp,Sys.time())
    print(diff(timeStamp))
	save(corMat,file="corMat.RData")

	summary(abs(c(corMat[lower.tri(corMat)])))
	#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
	#0.00000 0.01226 0.02827 0.07226 0.06898 1.00000 

	table(abs(c(corMat[lower.tri(corMat)]))==1)
    #corMat[lower.tri(corMat)][abs(c(corMat[lower.tri(corMat)]))==1]
    #corMat[abs(corMat)==1]
    
    png(paste("pairwiseCorFeature_kendall",cohort,".png",sep=""))
    plot(sort(abs(c(corMat[lower.tri(corMat)]))),main=cohortName)
    dev.off()

	n=ncol(cell)
	tmp=vector(mode="numeric",length=n*(n-1)/2)
	tmpC=vector(mode="character",length=n*(n-1)/2)
	corInfo=data.frame(feature1=tmpC,featUniq1=tmpC,feature2=tmpC,featUniq2=tmpC,corKend=tmp,stringsAsFactors=F)
	k=1
	for (k1 in 1:(nrow(corMat)-1)) {
		for (k2 in (k1+1):nrow(corMat)) {
			corInfo$feature1[k]=rownames(corMat)[k1]
			corInfo$feature2[k]=rownames(corMat)[k2]
            corInfo$featUniq1[k]=ann$featUniq[match(rownames(corMat)[k1],ann$feature)]
            corInfo$featUniq2[k]=ann$featUniq[match(rownames(corMat)[k2],ann$feature)]
            corInfo$corKend[k]=corMat[k1,k2]
			k=k+1
		}
	}
	write.table(corInfo,file=paste("corrFeature",cohort,".txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)
}

datadir="results/"
datadir=""
sumInfo=read.table(paste(datadir,"summaryFeature",cohort,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
load(file=paste(datadir,"corMat.RData",sep=""))
load(file=paste(datadir,"corMatSam.RData",sep=""))
corInfo=read.table(paste(datadir,"corrFeature",cohort,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

corInfo[abs(corInfo$corKend)==1,]
#feature1              feature2 corKend
#56059          NucArea      NucEquivDiameter       1
#56077          NucArea           NucZernike1       1
#56761 NucEquivDiameter           NucZernike1       1
#80741      MitoMaxArea      MitoFusedMaxArea       1
#80810      MitoMinArea MitoFragmentedMinArea       1

table(highCorr=abs(corInfo$corKend)>=.8,sameFeat=corInfo$featUniq1==corInfo$featUniq2)
"
            sameFeat
highCorr  FALSE  TRUE
    FALSE 18222  2427
    TRUE     26   235
"
xlim=c(0,1); ylim=c(0,4000)
png(paste("corKend_featureGroup",cohort,".png",sep=""),width=2*480,height=480)
par(mfrow=c(1,2))
hist(abs(corInfo$corKend[which(corInfo$featUniq1==corInfo$featUniq2)]),xlim=xlim,ylim=ylim,main=paste(cohortName,": Features from same group",sep=""),xlab="Absolute Kendall's correlation coefficient")
hist(abs(corInfo$corKend[which(corInfo$featUniq1!=corInfo$featUniq2)]),xlim=xlim,ylim=ylim,main=paste(cohortName,": Features from different groups",sep=""),xlab="Absolute Kendall's correlation coefficient")
dev.off()

## -------------------

featId1=which(!ann$feature%in%c("NucArea","NucEquivDiameter","MitoFusedMaxArea","MitoFragmentedMinArea"))

quantile(abs(c(cell[,featId1])),probs=seq(0,1,by=.1),na.rm=T)

## -------------------
## Correlation plots

subsetFlag="_zernike1Only"
for (typeFlag in sort(unique(ann$type))) {
	header=typeFlag	
	featId=featId1[ann$type[featId1]==typeFlag]
	annFeat=ann[featId,]
	x=corInfo[which(corInfo$feature1%in%annFeat$feature & corInfo$feature2%in%annFeat$feature),]
	if (subsetFlag=="_zernike1Only") {
		header=paste(header,", no zernike2...",sep="")
		i=grep("Zernike",annFeat$feature)
		i=i[which(!annFeat$feature[i]%in%c("CellZernike1","NucZernike1"))]
		if (length(i)!=0) {
			annFeat=annFeat[-i,]
			x=corInfo[which(corInfo$feature1%in%annFeat$feature & corInfo$feature2%in%annFeat$feature),]
		}
	}
	subDir <- paste(typeFlag,subsetFlag,sep="")			
	if (!file.exists(subDir)){
		dir.create(file.path(subDir))
	}
	subDir=paste(subDir,"/",sep="")
	png(paste(subDir,"corr_feature",cohort,"_",typeFlag,subsetFlag,"_%03d.png",sep=""), width = 6*480, height = 3*480)
#	par(mar=c(5, 4, 4, 2) + 0.1)
	par(mar=c(5, 5, 4, 1) + 0.1)
	par(mfrow=c(3,6))
	par(cex.main=3,cex.lab=3,cex.axis=1.5)
	for (k in order(abs(x$corKend),decreasing=T)[1:min(nrow(x),180)]) {
		j1=which(ann$feature==x$feature1[k])
		j2=which(ann$feature==x$feature2[k])
		plot(cell[,j1],cell[,j2],main=paste("Kendall corr: ",signif(x$corKend[k],2),sep=""),xlab=ann$feature[j1],ylab=ann$feature[j2])
	}
	dev.off()
}

grpUniq=unique(sapply(ann$feature[grep("2",ann$feature)],function(x) {
	y=strsplit(x,"2")[[1]][1]
	z=strsplit(y,"")[[1]]
	for (i in length(z):1) {
		if (is.na(as.integer(z[i]))) {
			break
		} else {
			y=substr(y,1,i-1)
		}
	}
	y
},USE.NAMES=F))
for (typeFlag in sort(unique(ann$type))) {
	for (gId in 1:length(grpUniq)) {
		header=paste(typeFlag,", ",grpUniq[gId],"s only",sep="")
		featId=featId1[ann$type[featId1]==typeFlag]
		annFeat=ann[featId,]
		
		i=grep(grpUniq[gId],annFeat$feature)
		if (length(i)==0) next
		annFeat=annFeat[i,]
		x=corInfo[which(corInfo$feature1%in%annFeat$feature & corInfo$feature2%in%annFeat$feature),]

		subDir <- paste(typeFlag,"_",grpUniq[gId],sep="")			
		if (!file.exists(subDir)){
			dir.create(file.path(subDir))
		}
		subDir=paste(subDir,"/",sep="")
		png(paste(subDir,"corr_feature",cohort,"_",typeFlag,"_",grpUniq[gId],"_%03d.png",sep=""), width = 6*480, height = 3*480)
		par(mar=c(5, 5, 4, 1) + 0.1)
		par(mfrow=c(3,6))
		par(cex.main=3,cex.lab=3,cex.axis=1.5)
		for (k in order(abs(x$corKend),decreasing=T)[1:min(nrow(x),180)]) {
			j1=which(ann$feature==x$feature1[k])
			j2=which(ann$feature==x$feature2[k])
			plot(cell[,j1],cell[,j2],main=paste("Kendall corr: ",signif(x$corKend[k],2),sep=""),xlab=ann$feature[j1],ylab=ann$feature[j2])
		}
		dev.off()
	}
}


####################################################################
####################################################################
## NOT USED. Check
## Heatmap

i=which(ann$feature=="CellSumIntensity")
summary(cell[,i])
summary(cell[,i]-median(cell[,i],na.rm=T))
i=which(ann$feature=="CellMeanExteriorChordLength")
summary(cell[,i])
summary(cell[,i]-median(cell[,i],na.rm=T))


i=which(annFeat$feature=="CellSumIntensity")
summary(arrayData[i,])

## -------------------
"
> table(phen2$dist2class1)

1    2    3    4    5
778 1351 1242  213  246

> table(phen2$dist2class2)

1    2    3    4    5
632 1161 1566  258  213

> table(phen2$dist2class3)

1    2    3    4    5
1858 1304  532   92   44

> table(phen2$dist2class4)

1    2    3    4    5
313   13  252 2604  648

> table(phen2$dist2class5)

1    2    3    4    5
249    1  238  663 2679

> nrow(phen2)
[1] 3830
"

## -------------------
## Which features contribute most to the sample clustering
## For each feature
## * Something like t-test
## * Deviation from the centroids of the sample clusters - WSS
## Then rank the features
## For myc/ras+wt, all features, do this for 4 and 7 sample clusters

library(coin)

cohort="_mycRasWt"
typeflag=""
nClust=4

datadir1=paste("results/heatmap/",sub("_","",cohort),"/",sep="")

cohortName="Myc/Ras & Wildtype"
cell=cellMW
cell_rbf=cell_rbfMW

arrayData=t(cell_rbf)
i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
arrayData[i2,]=t(cell[,i1])
annFeat=ann_rbf

centr=apply(arrayData,1,median,na.rm=T)
for (i in 1:nrow(arrayData)) {
    arrayData[i,]=arrayData[i,]-centr[i]
}

scal=apply(arrayData,1,sd,na.rm=T)
for (i in 1:nrow(arrayData)) {
    arrayData[i,]=arrayData[i,]/scal[i]
}

## ----
distribFlag="approximate"
distribFlag="asymptotic"
pvalMat=matrix(nrow=nrow(arrayData),ncol=4,dimnames=list(rownames(arrayData),paste("_",c(paste("4hclust",sep=""),paste(2:4,"kmeans",sep="")),sep="")))
for (gId in 1) {
    nClust=as.integer(substr(colnames(pvalMat)[gId],2,2))
    if (length(grep("hclust",colnames(pvalMat)[gId]))==1) {
        datadir=paste(datadir1,"wtOrd/reducedBioFeatPC/_mycRasWt_wtOrd_reducedBioFeatPC_kendall/",sep="")
        clustInfo=read.table(paste(datadir,"clusterInfoSample",cohort,typeflag,"_wtOrd_reducedBioFeatPC_kendall_",nClust,"clust.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    } else {
        datadir=paste(datadir1,"kmeans/",sep="")
        clustInfo=read.table(paste(datadir,"clusterInfoSample",cohort,typeflag,"_kmeans",nClust,"ClustOrd_reducedBioFeatPC_kendall.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    }
    j=match(colnames(arrayData),clustInfo$id); j1=which(!is.na(j)); j2=j[j1]
    x=as.factor(clustInfo$clustId[j2])
    pvalMat[,gId]=apply(arrayData[,j1],1,function(y,x) {pvalue(kruskal_test(y~x,distribution=distribFlag))}, x=x)
}
tbl=cbind(feature=rownames(pvalMat),as.data.frame(pvalMat))
rownames(tbl)=NULL
write.table(tbl, paste("pvalFeature_forClusterSample",cohort,typeflag,"_reducedBioFeatPC_kendall.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)


pvalCutoff=rep(NA,ncol(pvalMat))
lim=range(c(pvalMat),na.rm=T)
png(paste("pvalFeature_forClusterSample",cohort,typeflag,"_reducedBioFeatPC_kendall.png",sep=""))
par(mfrow=c(2,2))
for (gId in 1) {
    lim=c(0,.5)
    y=sort(pvalMat[,gId],decreasing=F)
    plot(1:nrow(pvalMat), y, ylim=lim, type="b", main=paste(sub("_","",colnames(pvalMat)[gId]),sep=""), xlab="Feature (ordered by p-value)",ylab="P-value")
    y2=y
    while (T) {
        i=which.min(diff(y2))
        if (i<.75*length(y)) {
            break
        }
        y2=y[1:i]
    }
    abline(h=y[i],lty="dotted")
    pvalCutoff[gId]=i
}
dev.off()

cbind(ord=1:length(y),y,diff=c(diff(y),1))
which.min(diff(y[1:i]))

## ----
wssMat=matrix(0,nrow=nrow(arrayData),ncol=4,dimnames=list(rownames(arrayData),paste("_",c(paste("4hclust",sep=""),paste(2:4,"kmeans",sep="")),sep="")))
for (gId in 1:ncol(wssMat)) {
    nClust=as.integer(substr(colnames(wssMat)[gId],2,2))
    if (length(grep("hclust",colnames(wssMat)[gId]))==1) {
        datadir=paste(datadir1,"wtOrd/reducedBioFeatPC/_mycRasWt_wtOrd_reducedBioFeatPC_kendall/",sep="")
        clustInfo=read.table(paste(datadir,"clusterInfoSample",cohort,typeflag,"_wtOrd_reducedBioFeatPC_kendall_",nClust,"clust.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    } else {
        datadir=paste(datadir1,"kmeans/",sep="")
        clustInfo=read.table(paste(datadir,"clusterInfoSample",cohort,typeflag,"_kmeans",nClust,"ClustOrd_reducedBioFeatPC_kendall.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        #clustInfo$clustId=paste("cluster",clustInfo$clustId,sep="")
        #write.table(clustInfo, paste("clusterInfoSample",cohort,typeflag,"_kmeans",nClust,"ClustOrd_reducedBioFeatPC_kendall.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
        
    }
    
    centr=matrix(nrow=nrow(arrayData),ncol=nClust,dimnames=list(rownames(arrayData),paste("cluster",1:nClust,sep="")))
    for (k in 1:nClust) {
        j2=which(clustInfo$clustId==paste("cluster",k,sep=""))
        j1=match(clustInfo$id[j2],colnames(arrayData))
        centr[,k]=apply(arrayData[,j1],1,mean,na.rm=T)
        for (i in 1:nrow(arrayData)) {
            wssMat[i,gId]=wssMat[i,gId]+sum((arrayData[i,j1]-centr[i,k])^2,na.rm=T)
        }
    }
}
tbl=cbind(feature=rownames(wssMat),as.data.frame(wssMat))
rownames(tbl)=NULL
write.table(tbl, paste("wssFeature_forClusterSample",cohort,typeflag,"_reducedBioFeatPC_kendall.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

wssCutoff=rep(NA,ncol(wssMat))
lim=range(c(wssMat),na.rm=T)
png(paste("wssFeature_forClusterSample",cohort,typeflag,"_reducedBioFeatPC_kendall.png",sep=""))
par(mfrow=c(2,2))
for (gId in 1:ncol(wssMat)) {
    y=sort(wssMat[,gId],decreasing=T)
    plot(1:nrow(wssMat), y, ylim=lim, type="b", main=paste(sub("_","",colnames(wssMat)[gId]),sep=""), xlab="Feature (ordered by WSS)",ylab="Within groups sum of squares")
    y2=y
    while (T) {
        i=which.min(diff(y2))
        if (i<.75*length(y)) {
            break
        }
        y2=y[1:i]
    }
    abline(h=y[i],lty="dotted")
    wssCutoff[gId]=i
}
dev.off()

cbind(ord=1:length(y),y,diff=c(diff(y),1))
which.min(diff(y[1:i]))


if (F) {
wssCutoff=rep(NA,ncol(wssMat))
lim=range(c(wssMat),na.rm=T)
png(paste("wssFeature_forClusterSample",cohort,typeflag,"_reducedBioFeatPC_kendall.png",sep=""))
par(mfrow=c(2,2))
for (gId in 1:ncol(wssMat)) {
    y=sort(wssMat[,gId],decreasing=F)
    plot(1:nrow(wssMat), y, ylim=lim, type="b", main=paste(sub("_","",colnames(wssMat)[gId]),sep=""), xlab="Feature (ordered by WSS)",ylab="Within groups sum of squares")
    y2=y
    i=0
    while (T) {
        i=((i+1):length(y))[which.max(diff(y2))]
        if (i<.75*length(y)) {
            break
        }
        y2=y[(i+1):length(y)]
    }
    abline(h=y[i],lty="dotted")
    wssCutoff[gId]=i
}
dev.off()

cbind(ord=1:length(y),y,diff=c(diff(y),1))
which.max(diff(y[1:i]))
}

## -------------------
library(fpc)

classifFlag="averagedist"
classifFlag="centroid"

datTypeFlag="sample"
datTypeFlag="feature"

#for (datTypeFlag in c("sample","feature")) {
for (datTypeFlag in c("sample")) {
    datadir2=paste("predictionStrength/",classifFlag,"/",sep="")
    fileList=dir(datadir2)
    if (datTypeFlag=="feature") {
        fId=grep("predictionStrength_",fileList)
    } else {
        fId=grep("predictionStrengthSam_",fileList)
    }
    fileList=fileList[fId]
    fileList=fileList[grep("_mycRasWt",fileList)]
    fileList=fileList[grep("_2clustMin",fileList)]
    cat("\n\n===============================================\n")
    cat("============ ",datTypeFlag," ================\n")
    cat("===============================================\n")

    for (fId in 1:length(fileList)) {
        fNameOut=sub(".RData","",fileList[fId],fixed=T)
        x=c(sapply(fNameOut,function(x){
            y=strsplit(x,"_")[[1]]
            if (y[3]!="mean") y=c(y[1:2],"min",y[3:length(y)])
            y
        },USE.NAMES=F))
        names(x)=c("predictionStrengthSam","centroid","aggregation","2clustMin","15clustMax","20perms","cohort","reducedBioFeatPC","distMethod")
        x["cohort"]=paste("_",x["cohort"],sep="")
        
        cohort=x["cohort"]
        distMethod=x["distMethod"]
        
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
        
        #heading=paste(cohortName,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),subsetFlag,centrFlag,scaleFlag,", ",distMethod,", sd ",sub("_sd","",stdDev),sep="")
        heading=paste(cohortName,", ",distMethod,", ",x["aggregation"],sep="")
        cat("\n\n==========================================\n")
        #cat("\n\n============ ",heading," ================\n")
        print(fileList[fId])
        load(paste(datadir2,fileList[fId],sep=""))
        cat("Optimum k: ",res$optimalk,"\n")
        
        x=unlist(lapply(res$predcorr,function(x) {
            if (is.null(x)) y=NA else y=mean(x,na.rm=T)
            y
        }))
        names(x)=paste(1:(length(x)),"clust",sep="")
        cat("Mean prediction strength\n")
        print(x)
        
        lim=c(0,1)
        png(paste(fNameOut,".png",sep=""))
        plot((1:length(res$predcorr)), x, ylim=lim, type="b", main=paste(heading,sep=""), xlab="Number of clusters",ylab="Prediction strength")
        abline(h=0.8,lty="dotted")
        dev.off()
    }
}
"
===============================================
============  sample  ================
===============================================


============  mycRas  ================
[1] "predictionStrengthSam_centroid_mycRas_cell_reducedBioFeatPC_kendall.RData"
Optimum k:  1


============  mycRas  ================
[1] "predictionStrengthSam_centroid_mycRas_mitochondria_reducedBioFeatPC_kendall.RData"
Optimum k:  1


============  mycRas  ================
[1] "predictionStrengthSam_centroid_mycRas_nucleus_reducedBioFeatPC_kendall.RData"
Optimum k:  1


============  mycRas  ================
[1] "predictionStrengthSam_centroid_mycRas_reducedBioFeatPC_kendall.RData"
Optimum k:  1


============  wt  ================
[1] "predictionStrengthSam_centroid_wt_cell_reducedBioFeatPC_kendall.RData"
Optimum k:  1


============  wt  ================
[1] "predictionStrengthSam_centroid_wt_reducedBioFeatPC_kendall.RData"
Optimum k:  1


===============================================
============  feature  ================
===============================================


============  mycRas, cell  ================
[1] "predictionStrength_centroid_mycRas_cell_reducedBioFeatPC_kendall_cell.RData"
Optimum k:  2


============  mycRas, mitochondria  ================
[1] "predictionStrength_centroid_mycRas_mitochondria_reducedBioFeatPC_kendall_mitochondria.RData"
Optimum k:  1


============  mycRas, nucleus  ================
[1] "predictionStrength_centroid_mycRas_nucleus_reducedBioFeatPC_kendall_nucleus.RData"
Optimum k:  2


============  mycRas  ================
[1] "predictionStrength_centroid_mycRas_reducedBioFeatPC_kendall.RData"
Optimum k:  1


============  wt, cell  ================
[1] "predictionStrength_centroid_wt_cell_reducedBioFeatPC_kendall_cell.RData"
Optimum k:  2


============  wt  ================
[1] "predictionStrength_centroid_wt_reducedBioFeatPC_kendall.RData"
Optimum k:  2
"

####################################################################
####################################################################
## Section 2
## -------------------
## Reduce features by PCA
## kindFlag="biology" - based on biology
## kindFlag="cluster" - based on clustering
##     Separate out the feature clusters with good profiles from those with bad profiles
##     Use 15 clusters

cohort="_mycRas"
cohort="_wt"
cohort="_mycRasWt"

for (cohort in c("_wt906","_mycRas289")) {
    #for (cohort in c("_mycRas","_wt","_mycRasWt","_wt906","_mycRas289")) {
    kindFlag="cluster"
    kindFlag="biology"

    switch(cohort,
        "_mycRas"={
           cohortName="Myc/Ras"
           cell=cellM
        },
        "_wt"={
           cohortName="Wildtype"
           cell=cellW
        },
        "_mycRasWt"={
           cohortName="Myc/Ras & Wildtype"
           cell=rbind(cellM,cellW)
           #rownames(cell)=c(paste("m",rownames(cellM),sep=""),paste("w",rownames(cellW),sep=""))
        },
        "_wt906"={
           cohortName="Wildtype 906"
           cell=cellW2
        },
        "_mycRas289"={
            cohortName="Myc/Ras 289"
            cell=cellM2
        }
    )
    if (cohort%in%c("_mycRas","_wt","_mycRasWt")) {
        ann=annW
        #ann_rbf=ann_rbf1
        #ann_rbfm=ann_rbfm1
    } else {
        ann=annW2
        #ann_rbf=ann_rbf2
        #ann_rbfm=ann_rbfm2
    }
    
    if (kindFlag=="cluster") {
        good=list(c(1,6,8,9,11),c(1,2,3,4,5,6,8,9,10,14,15),c(1,2,3,10,14,15),c(1,4,5,6,7,8,9,10,11,14),c(1,2,3,6,8,13))
        names(good)=c("cell","mitochondria","nucleus","cell_zernike1Only","nucleus_zernike1Only")
        grpUniq=names(good)

        n=15
        datadir=paste("results/heatmap/initial/",n,"trees/clusterInfo/",sep="")
        fileList=dir(datadir)
        fileList=fileList[grep(cohort,fileList)]
    } else {
        grpUniq=unique(ann$type)
    }

    arrayData=cell
    for (j in 1:ncol(arrayData)) {
        meanThis=mean(arrayData[,j],na.rm=T)
        arrayData[is.na(arrayData[,j]),j]=meanThis
    }
    arrayData1=arrayData

    arrayData=cell
    centr=apply(arrayData,2,median,na.rm=T)
    for (i in 1:ncol(arrayData)) {
        arrayData[,i]=arrayData[,i]-centr[i]
    }
    scal=apply(arrayData,2,sd,na.rm=T)
    for (i in 1:ncol(arrayData)) {
        arrayData[,i]=arrayData[,i]/scal[i]
    }
    arrayData2=arrayData
    for (j in 1:ncol(arrayData)) {
        meanThis=mean(arrayData[,j],na.rm=T)
        arrayData[is.na(arrayData[,j]),j]=meanThis
    }
    arrayData2_1=arrayData

    cell_rbf=ann_rbf=NULL
    cell_rbfm=ann_rbfm=NULL

    cell_rf=ann_rf=NULL
    cell_rfm=ann_rfm=NULL
    cell_nzrf=ann_nzrf=NULL
    cell_nzrfm=ann_nzrfm=NULL

    for (gId in 1:length(grpUniq)) {
        cat("\n\n===============",grpUniq[gId],"===========\n")
        typeThis=sub("_zernike1Only","",grpUniq[gId])
        if (kindFlag=="cluster") {
            #clustInfo=read.table(paste(datadir,fileList[gId],sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            clustInfo=read.table(paste(datadir,"clusterInfoFeature",cohort,"_",grpUniq[gId],"_kendall.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            clustInfo$good=0
            k=which(clustInfo$clustId%in%paste("cluster",good[[grpUniq[gId]]],sep=""))
            clustInfo$good[k]=1
        } else {
            clustInfo=ann[which(ann$type==grpUniq[gId]),]
            clustInfo$clustId=clustInfo$featUniq
            if (cohort%in%c("_wt906","_mycRas289")) {
                clustInfo$good=1
            } else {
                clustInfo$good=0
                ## CHECK !!! Not sure why unique features with no replicates are excluded
                clustInfo$good[which(clustInfo$clustId%in%clustInfo$clustId[duplicated(clustInfo$clustId)])]=1
            }
        }
        
        clustUniq=unique(clustInfo$clustId[which(clustInfo$good==1)])
        cellMetaMean=cellMeta=matrix(nrow=nrow(arrayData2),ncol=length(clustUniq),dimnames=list(rownames(arrayData2),clustUniq))
        for (cId in 1:length(clustUniq)) {
            cat("\n\n---------------",clustUniq[cId],"\n")
            j=which(colnames(arrayData)%in%clustInfo$feature[which(clustInfo$clustId==clustUniq[cId])])
            #if (length(j)==0) {
            #    cellMetaMean[,clustUniq[cId]]=cellMeta[,clustUniq[cId]]=NA
            #} else if (length(j)==1) {
            if (length(j)==1) {
                cellMetaMean[,clustUniq[cId]]=cellMeta[,clustUniq[cId]]=arrayData2[,j]
            } else {
                j=j[apply(arrayData2_1[,j],2,function(x) mean(!is.na(x)))]
                cellMetaMean[,clustUniq[cId]]=apply(arrayData2[,j],1,mean,na.rm=T)
                #fit=prcomp(arrayData[,j], center=T, scale=T)
                fit=prcomp(arrayData2_1[,j], center=F, scale=F)
                cellMeta[,clustUniq[cId]]=fit$x[,1]
                png(paste("plot_meanVpca1_cellMeta","_",clustUniq[cId],cohort,"_",grpUniq[gId],".png",sep=""))
                plot(apply(arrayData2_1[,j],1,mean,na.rm=T),fit$x[,1],main=paste(grpUniq[gId],": ",clustUniq[cId],sep=""),xlab="mean",ylab="pca1")
                dev.off()
                png(paste("screeplot_cellMeta","_",clustUniq[cId],cohort,"_",grpUniq[gId],".png",sep=""))
                plot(fit,main=paste(grpUniq[gId],": ",clustUniq[cId],sep=""))
                dev.off()
            }
            print(ann$feature[j])
        }
        
        colId=c("feature","type")
        j=which(ann$feature%in%clustInfo$feature[which(clustInfo$good==0)])
        if (length(grep("_zernike1Only",grpUniq[gId]))==1) {
            cell_nzrf=cbind(cell_nzrf,arrayData2_1[,j],cellMeta)
            ann_nzrf=rbind(ann_nzrf,ann[j,colId],data.frame(feature=colnames(cellMeta),type=typeThis,stringsAsFactors=F))
            cell_nzrfm=cbind(cell_nzrfm,arrayData2_1[,j],cellMetaMean)
            ann_nzrfm=rbind(ann_nzrfm,ann[j,colId],data.frame(feature=colnames(cellMetaMean),type=typeThis,stringsAsFactors=F))
            
            tbl=cbind(ann_nzrf,t(cell_nzrf))
            names(tbl)=c(names(ann_nzrf),paste("sam",1:nrow(cell_nzrf),sep=""))
            write.table(tbl, paste("cell_zernike1Only_reducedFeatByPCA",cohort,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            tbl=cbind(ann_nzrfm,t(cell_nzrfm))
            names(tbl)=c(names(ann_nzrfm),paste("sam",1:nrow(cell_nzrfm),sep=""))
            write.table(tbl, paste("cell_zernike1Only_reducedFeatByMean",cohort,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
        } else {
            if (kindFlag=="cluster") {
                #cell_rf=cbind(cell_rf,cell[,j],cellMeta)
                cell_rf=cbind(cell_rf,arrayData2_1[,j],cellMeta)
                ann_rf=rbind(ann_rf,ann[j,colId],data.frame(feature=colnames(cellMeta),type=typeThis,stringsAsFactors=F))
                #cell_rfm=cbind(cell_rfm,arrayData2[,j],cellMetaMean)
                cell_rfm=cbind(cell_rfm,arrayData2_1[,j],cellMetaMean)
                ann_rfm=rbind(ann_rfm,ann[j,colId],data.frame(feature=colnames(cellMetaMean),type=typeThis,stringsAsFactors=F))
                
                tbl=cbind(ann_rf,t(cell_rf))
                names(tbl)=c(names(ann_rf),paste("sam",1:nrow(cell_rf),sep=""))
                write.table(tbl, paste("cell_reducedFeatByPCA",cohort,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                tbl=cbind(ann_rfm,t(cell_rfm))
                names(tbl)=c(names(ann_rfm),paste("sam",1:nrow(cell_rfm),sep=""))
                write.table(tbl, paste("cell_reducedFeatByMean",cohort,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            } else {
                cell_rbf=cbind(cell_rbf,arrayData2_1[,j],cellMeta)
                ann_rbf=rbind(ann_rbf,ann[j,colId],data.frame(feature=colnames(cellMeta),type=typeThis,stringsAsFactors=F))
                cell_rbfm=cbind(cell_rbfm,arrayData2_1[,j],cellMetaMean)
                ann_rbfm=rbind(ann_rbfm,ann[j,colId],data.frame(feature=colnames(cellMetaMean),type=typeThis,stringsAsFactors=F))
                
                tbl=cbind(ann_rbf,t(cell_rbf))
                names(tbl)=c(names(ann_rbf),paste("sam",1:nrow(cell_rbf),sep=""))
                write.table(tbl, paste("cell_reducedBioFeatByPCA",cohort,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                tbl=cbind(ann_rbfm,t(cell_rbfm))
                names(tbl)=c(names(ann_rbfm),paste("sam",1:nrow(cell_rbfm),sep=""))
                write.table(tbl, paste("cell_reducedBioFeatByMean",cohort,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            }
        }
    }

    sdSamPc=apply(cell_rbf[,which(substr(colnames(cell_rbf),1,nchar("cluster"))!="cluster")],1,sd,na.rm=T)
    sdSamMn=apply(cell_rbfm[,which(substr(colnames(cell_rbfm),1,nchar("cluster"))!="cluster")],1,sd,na.rm=T)
    #sdSam=apply(arrayData2[,which(colnames(arrayData2)%in%colnames(cell_rfm))],1,sd,na.rm=T)
    #sdSam=apply(cell[,which(colnames(cell)%in%colnames(cell_rfm))],1,sd,na.rm=T)

    sdFeatPc=apply(cell_rbf[,which(substr(colnames(cell_rbf),1,nchar("cluster"))!="cluster")],2,sd,na.rm=T)
    sdFeatMn=apply(cell_rbfm[,which(substr(colnames(cell_rbfm),1,nchar("cluster"))!="cluster")],2,sd,na.rm=T)
    meanFeatPc=apply(cell_rbf[,which(substr(colnames(cell_rbf),1,nchar("cluster"))!="cluster")],2,mean,na.rm=T)
    meanFeatMn=apply(cell_rbfm[,which(substr(colnames(cell_rbfm),1,nchar("cluster"))!="cluster")],2,mean,na.rm=T)
    quarFeatPc=t(apply(cell_rbf[,which(substr(colnames(cell_rbf),1,nchar("cluster"))!="cluster")],2,quantile,probs=c(.25,.75),na.rm=T))
    quarFeatMn=t(apply(cell_rbfm[,which(substr(colnames(cell_rbfm),1,nchar("cluster"))!="cluster")],2,quantile,probs=c(.25,.75),na.rm=T))

    if (cohort%in%c("_mycRas","_wt","_mycRasWt")) {
        ann_rbf1=ann_rbf
        ann_rbfm1=ann_rbfm
    } else {
        ann_rbf2=ann_rbf
        ann_rbfm2=ann_rbfm
    }
    switch(cohort,
        "_mycRas"={
           sdSamPcM=sdSamPc
           sdSamMnM=sdSamMn
           sdFeatPcM=sdFeatPc
           sdFeatMnM=sdFeatMn
           meanFeatPcM=meanFeatPc
           meanFeatMnM=meanFeatMn
           quarFeatPcM=quarFeatPc
           quarFeatMnM=quarFeatMn
           cell_rbfM=cell_rbf
           cell_rbfmM=cell_rbfm
        },
        "_wt"={
            sdSamPcW=sdSamPc
            sdSamMnW=sdSamMn
            sdFeatPcW=sdFeatPc
            sdFeatMnW=sdFeatMn
            meanFeatPcW=meanFeatPc
            meanFeatMnW=meanFeatMn
            quarFeatPcW=quarFeatPc
            quarFeatMnW=quarFeatMn
            cell_rbfW=cell_rbf
            cell_rbfmW=cell_rbfm
        },
        "_mycRasWt"={
            sdSamPcMW=sdSamPc
            sdSamMnMW=sdSamMn
            sdFeatPcMW=sdFeatPc
            sdFeatMnMW=sdFeatMn
            meanFeatPcMW=meanFeatPc
            meanFeatMnMW=meanFeatMn
            quarFeatPcMW=quarFeatPc
            quarFeatMnMW=quarFeatMn
            cell_rbfMW=cell_rbf
            cell_rbfmMW=cell_rbfm
        },
        "_wt906"={
            sdSamPcW2=sdSamPc
            sdSamMnW2=sdSamMn
            sdFeatPcW2=sdFeatPc
            sdFeatMnW2=sdFeatMn
            meanFeatPcW2=meanFeatPc
            meanFeatMnW2=meanFeatMn
            quarFeatPcW2=quarFeatPc
            quarFeatMnW2=quarFeatMn
            cell_rbfW2=cell_rbf
            cell_rbfmW2=cell_rbfm
        },
        "_mycRas289"={
            sdSamPcM2=sdSamPc
            sdSamMnM2=sdSamMn
            sdFeatPcM2=sdFeatPc
            sdFeatMnM2=sdFeatMn
            meanFeatPcM2=meanFeatPc
            meanFeatMnM2=meanFeatMn
            quarFeatPcM2=quarFeatPc
            quarFeatMnM2=quarFeatMn
            cell_rbfM2=cell_rbf
            cell_rbfmM2=cell_rbfm
        }
    )
}

save.image("tmp.RData")

####################################################################
####################################################################
## Section

for (thisFlag in c("_pc","_mean")) {
    switch(thisFlag,
    "_pc"={
        meanFeat1=meanFeatPcM
        meanFeat2=meanFeatPcW
        sdFeat1=sdFeatPcM
        sdFeat2=sdFeatPcW
        quarFeat1=quarFeatPcM
        quarFeat2=quarFeatPcW
        header="Features reduced by PCA"
    },
    "_mean"={
        meanFeat1=meanFeatMnM
        meanFeat2=meanFeatMnW
        sdFeat1=sdFeatMnM
        sdFeat2=sdFeatMnW
        quarFeat1=quarFeatMnM
        quarFeat2=quarFeatMnW
        header="Features reduced by averaging"
    }
    )
    for (metricFlag in c("sd","cv","mean","diffMeanVar","qcd")) {
        switch(metricFlag,
        "sd"={
            x1=sdFeat1; x2=sdFeat2; nm="Feature Standard Deviation"
            x1=sdFeat1^2; x2=sdFeat2^2; nm="Feature Variance"
            x1=log(sdFeat1); x2=log(sdFeat2); nm="log(Feature Standard Deviation)"
            x1=log2(sdFeat1); x2=log2(sdFeat2); nm="log2(Feature Standard Deviation)"
        },
        "cv"={
            x1=sdFeat1/meanFeat1; x2=sdFeat2/meanFeat2; nm="Feature Coefficient of Variation"
            x1=log2(sdFeat1/abs(meanFeat1)); x2=log2(sdFeat2/abs(meanFeat2)); nm="log2(Feature Coefficient of Variation)"
        },
        "mean"={
            x1=meanFeat1; x2=meanFeat2; nm="Feature Mean"
            if (any(x1<=0 | x2<=0)) {
                x1=x1-min(c(x1))+1
                x2=x2-min(c(x2))+1
            }
            x1=log2(x1); x2=log2(x2); nm="log2(Feature Mean)"
        },
        "diffMeanVar"={
            x1=meanFeat2-meanFeat1; x2=sdFeat2^2-sdFeat1^2; nm="Feature: Wildtype Mean - Myc/Ras Mean"; nm="Feature: Wildtype Variance - Myc/Ras Variance"
            if (any(x1<=0 | x2<=0)) {
                x1=x1-min(c(x1))+1
                x2=x2-min(c(x2))+1
            }
            x1=log2(x1); x2=log2(x2); nm=c("Feature: log2(Wildtype Mean - Myc/Ras Mean)","Feature: log2(Wildtype Variance - Myc/Ras Variance)")
        },
        "qcd"={
            x1=apply(quarFeat1,1,diff)/apply(quarFeat1,1,sum); x2=apply(quarFeat2,1,diff)/apply(quarFeat2,1,sum); nm="Feature Quartile Coefficient of Dispersion"
            if (any(x1<=0 | x2<=0)) {
                x1=x1-min(c(x1))+1
                x2=x2-min(c(x2))+1
            }
            x1=log2(x1); x2=log2(x2); nm="log2(Feature Quartile Coefficient of Dispersion)"
        }
        )
        if (length(nm)==1) nm=paste(c("Myc/Ras","Wildtype"),": ",nm,sep="")
        header2=paste(header,"\nProportion higher in wildtype:",round(mean(x2>x1),2))
        png(paste(metricFlag,"Feat_mycRasVsWt",thisFlag,".png",sep=""))
        lim=range(c(x1,x2),na.rm=T)
        plot(x1,x2,xlim=lim,ylim=lim,main=header2,xlab=nm[1],ylab=nm[2])
        abline(c(0,1))
        dev.off()
    }
}

"
summary(sdFeatPcW-sdFeatPcM)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#-0.673800 -0.000002  0.000001  0.071610  0.000087  9.086000
mean(sdFeatPcW>sdFeatPcM)
#[1] 0.6595745
"

for (thisFlag in c("_mycRas","_wt")) {
    switch(thisFlag,
    "_mycRas"={
        sdFeat1=sdFeatPcM
        sdFeat2=sdFeatMnM
        header="Myc/Ras"
    },
    "_wt"={
        sdFeat1=sdFeatPcW
        sdFeat2=sdFeatMnW
        header="Wildtype"
    }
    )
    x1=log2(sdFeat1); x2=log2(sdFeat2); nm="log2(Feature Standard Deviation)"
    png(paste("sdFeat_pcVsMeanBasedFeatureReduction",thisFlag,".png",sep=""))
    lim=range(c(x1,x2),na.rm=T)
    plot(x1,x2,xlim=lim,ylim=lim,main=header,xlab=paste("Features reduced by PCA: ",nm,sep=""),ylab=paste("Features reduced by averaging: ",nm,sep=""))
    abline(c(0,1))
    dev.off()
}


## -------------------
## Test for association of clusters defined by the 3 feature types

cohort="_mycRas289"; orderFlag=""
cohort="_wt906"; orderFlag="_mycRas289Ord"

datadir=""

typeList=sort(unique(ann$type))
nClustList=2:10

switch(cohort,
    "_wt906"={
        cell=cellW2
        annCell=annCellW2
    },
    "_mycRas289"={
        cell=cellM2
        annCell=annCellM2
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
}

type2Flag=""
type2Flag="_reducedBioFeatPC"

out=matrix(nrow=nrow(annCell),ncol=length(nClustList)*length(typeList),dimnames=list(annCell$id,paste("clustId_",nClustList,"_",rep(typeList,each=length(nClustList)),sep="")))
for (typeFlag in typeList) {
    clustInfo=read.table(paste(datadir,cohort,orderFlag,"_",typeFlag,type2Flag,"_kendall/clusterInfoSample",cohort,orderFlag,"_",typeFlag,type2Flag,"_kendall.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    for (nClust in nClustList) {
        out[match(clustInfo$id,rownames(out)),paste("clustId_",nClust,"_",typeFlag,sep="")]=clustInfo[,paste("clustId_",nClust,sep="")]
    }
}

n=length(typeList)*(length(typeList)-1)/2
pvalMat=matrix(nrow=length(nClustList),ncol=n)
for (cId in 1:length(nClustList)) {
    nClust=nClustList[cId]
    nm=rep("",n)
    k=1
    for (k1 in 1:(length(typeList)-1)) {
        for (k2 in (k1+1):length(typeList)) {
            x=table(out[,paste("clustId_",nClust,"_",typeList[k1],sep="")],out[,paste("clustId_",nClust,"_",typeList[k2],sep="")])
            #print(chisq.test(x))
            ## Fisher's exact test giving error
            pvalMat[cId,k]=chisq.test(x)$p.value
            nm[k]=paste(typeList[k1],"_",typeList[k2],sep="")
            k=k+1
        }
    }
}
colnames(pvalMat)=nm
tbl=cbind(id=paste("nClust_",nClustList,sep=""),as.data.frame(pvalMat))
rownames(tbl)=NULL
write.table(tbl, paste("pvalFeature_forClusterSample",cohort,typeflag,"_reducedBioFeatPC_kendall.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

out=cbind(pvalMat,rep(1,nrow(pvalMat)))
ttl=rep("",nrow(out)*ncol(out))
ttl[seq(2,length(ttl),by=length(typeList)+1)]=paste("nClust ",nClustList,sep="")
png(paste("pvalFeature_forClusterSample",cohort,typeflag,"_reducedBioFeatPC_kendall.png",sep=""))
barplot(-log10(c(t(out))),ylab="-log10(p-value)",names.arg=ttl,col=colList[1:(length(typeList)+1)],las=3)
dev.off()

## Optimum concordance among the 3 feature types: 5 sample clusters

## ----------------------------------------

####################################################################
####################################################################
## PCA

cohort="_mycRas"
cohort="_wt"

cohort="_wt906"
cohort="_mycRas289"

featureFlag=""
featureFlag="_reducedFeatMean"
featureFlag="_reducedFeatPCA"
featureFlag="_reducedBioFeatMean"
featureFlag="_reducedBioFeatPCA"

subsetList=c("","_zernike1Only")
subsetList=c("")

plotCutoffFlag=T
plotCutoffFlag=F

for (cohort in c("_wt906","_mycRas289")) {
    #cohortName=ifelse(cohort=="_mycRas","Myc/Ras","Wildtype")
    switch(cohort,
        "_mycRas"={
           cohortName="Myc/Ras"
           cell=cellM
        },
        "_wt"={
           cohortName="Wildtype"
           cell=cellW
        },
        "_mycRasWt"={
            cohortName="Myc/Ras & Wildtype"
            cell=cellMW
        },
        "_wt906"={
            cohortName="Wildtype 906"
            cell=cellW2
        },
        "_mycRas289"={
            cohortName="Myc/Ras 289"
            cell=cellM2
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
    }

    for (typeFlag in sort(unique(ann$type))) {
        for (subsetFlag in subsetList) {
            fNameOut=paste(cohort,"_",typeFlag,subsetFlag,featureFlag,sep="")
            header=paste(cohortName,": ",typeFlag,sep="")
            cat("\n\n===========",header,"\n")

            if (featureFlag=="_reducedFeatPCA") {
                featId=which(ann_rf$type==typeFlag)			
                arrayData=cell_rf[,featId]
                annFeat=ann_rf[featId,]
            } else if (featureFlag=="_reducedFeatMean") {
                featId=which(ann_rfm$type==typeFlag)
                arrayData=cell_rfm[,featId]
                annFeat=ann_rfm[featId,]
            } else if (featureFlag=="_reducedBioFeatPCA") {
                featId=which(ann_rbf$type==typeFlag)			
                arrayData=cell_rbf[,featId]
                annFeat=ann_rbf[featId,]			
            } else if (featureFlag=="_reducedBioFeatMean") {
                featId=which(ann_rbfm$type==typeFlag)			
                arrayData=cell_rbfm[,featId]
                annFeat=ann_rbfm[featId,]						
            } else {
                featId=featId1[ann$type[featId1]==typeFlag]		
                arrayData=cell[,featId]
                annFeat=ann[featId,]
                corMat2=corMat[featId,featId]
                corMatSam2=corMatSam
            }

            if (subsetFlag=="_zernike1Only") {
                break
                
                i=grep("Zernike",annFeat$feature)
                i=i[which(!annFeat$feature[i]%in%c("CellZernike1","NucZernike1"))]
                if (length(i)==0) next
                arrayData=arrayData[,-i]
                annFeat=annFeat[-i,]
                corMat2=corMat2[-i,-i]
                header=paste(header,", no zernike2...",sep="")
            }

            rownames(arrayData)=paste("sam",1:nrow(arrayData),sep="")
            
            for (j in 1:ncol(arrayData)) {
                meanThis=mean(arrayData[,j],na.rm=T)
                arrayData[is.na(arrayData[,j]),j]=meanThis
            }

            fit=prcomp(arrayData, center=T, scale=T)
            print("length(fit$sdev)")
            print(length(fit$sdev))
            print("(1:length(fit$sdev))[order(c(diff(fit$sdev),999))]")
            print((1:length(fit$sdev))[order(c(diff(fit$sdev),999))])
            
            cutoff=3
            
            if (featureFlag=="_reducedFeatPCA") {
                switch(typeFlag,
                       "cell"={
                       fit1=fit
                       cutoff=12
                       },
                       "mitochondria"={
                       fit2=fit
                       cutoff=16
                       },
                       "nucleus"={
                       fit3=fit
                       cutoff=16
                       }
                       )
            }
            
            if (featureFlag=="_reducedFeatMean") {
                ## CHECK THIS !!!
                switch(typeFlag,
                       "cell"={
                       fit1=fit
                       cutoff=18
                       },
                       "mitochondria"={
                       fit2=fit
                       cutoff=9
                       },
                       "nucleus"={
                       fit3=fit
                       cutoff=18
                       }
                       )
            }
            #fit$rotation[,"PC1"]
            
            png(paste("screePlot_pca",fNameOut,".png",sep=""))
            screeplot(fit,npcs=length(fit$sdev),main=paste("PCA screeplot: ",header,sep=""))
            if (plotCutoffFlag) abline(v=cutoff+1,col="red",lty="dotted")
            dev.off()
            
            if (F) {
            for (k in c(10,20,30)) {
                png(paste("screePlot_pca_",k,"comp",fNameOut,".png",sep=""))
                screeplot(fit,npcs=k,main=paste("PCA screeplot (",k," comp): ",header,sep=""))
                dev.off()
            }
            }
            
            png(paste("biPlot_pca",fNameOut,"_%1d.png",sep=""))
            #par(mfcol=c(2,2))
            biplot(fit,choices=c(1,2),main=paste("PCA biplot: ",header,sep=""))
            biplot(fit,choices=c(1,3),main=paste("PCA biplot: ",header,sep=""))
            biplot(fit,choices=c(2,3),main=paste("PCA biplot: ",header,sep=""))
            dev.off()
            
            png(paste("rotationPlot_pca",fNameOut,".png",sep=""))
            par(mfcol=c(2,2))
            plot(fit$rotation[,"PC1"],fit$rotation[,"PC2"],main=paste("PCA: ",header,sep=""),xlab="Rotation: PC1",ylab="Rotation: PC2")
            plot(fit$rotation[,"PC1"],fit$rotation[,"PC3"],main=paste("PCA: ",header,sep=""),xlab="Rotation: PC1",ylab="Rotation: PC3")
            plot(fit$rotation[,"PC2"],fit$rotation[,"PC3"],main=paste("PCA: ",header,sep=""),xlab="Rotation: PC2",ylab="Rotation: PC3")
            dev.off()
            
            png(paste("scorePlot_pca",fNameOut,".png",sep=""))
            par(mfcol=c(2,2))
            plot(fit$x[,"PC1"],fit$x[,"PC2"],main=paste("PCA: ",header,sep=""),xlab="Score: PC1",ylab="Score: PC2")
            plot(fit$x[,"PC1"],fit$x[,"PC3"],main=paste("PCA: ",header,sep=""),xlab="Score: PC1",ylab="Score: PC3")
            plot(fit$x[,"PC2"],fit$x[,"PC3"],main=paste("PCA: ",header,sep=""),xlab="Score: PC2",ylab="Score: PC3")
            dev.off()
            
            tbl=cbind(feature=rownames(fit$x),fit$x[,1:cutoff])
            write.table(tbl, paste("prinComp",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
        }
    }
}

## -------------------
## NOT USED

x1=apply(arrayData,2,function(x) x-median(x,na.rm=T))

png("pca_0.png")
biplot(fit)
dev.off()

fit0=prcomp(-arrayData, center=T, scale=T)
png("pca_0.png")
biplot(fit0,xlabs=rep("*",nrow(arrayData)),ylabs=rep("*",ncol(arrayData)))
dev.off()

fitF=prcomp(t(x1), center=F, scale=F)
png("pca_f_biplot.png")
biplot(fitF,xlabs=rep("*",ncol(arrayData)),ylabs=rep("*",nrow(arrayData)))
dev.off()

fit1=princomp(x1, cor = TRUE)
png("pca_1_biplot.png")
biplot(fit1,xlabs=rep("*",nrow(arrayData)),ylabs=rep("*",ncol(arrayData)))
dev.off()
png("pca_1_score.png")
lim=NULL
plot(fit1$scores[,1],fit1$scores[,2],xlim=lim,ylim=lim)
abline(c(0,1),lty="dotted")
dev.off()
png("pca_1_eigenvector.png")
lim=NULL
plot(fit1$loadings[,1],fit1$loadings[,2],xlim=lim,ylim=lim)
abline(c(0,1),lty="dotted")
dev.off()

library(FactoMineR)

fit2 <- PCA(x1,scale.unit = TRUE, graph = F, ncp =6)
pc <- fit2$var$coord[,c(1:6)]
pc <- fit2$ind$coord[,c(1:6)]
png("pca_2.png")
lim=c(-2000,1000)
lim=NULL
plot(fit2,xlim=lim,ylim=lim)
abline(c(0,1),lty="dotted")
dev.off()
png("pca_2_2.png")
lim=c(-2000,1000)
lim=NULL
plot(fit2$ind$coord[,1],fit2$ind$coord[,2],xlim=lim,ylim=lim)
abline(c(0,1),lty="dotted")
dev.off()
png("pca_2_3.png")
lim=c(-2000,1000)
lim=NULL
plot(fit2$svd$U[,1],fit2$svd$U[,2],xlim=lim,ylim=lim)
abline(c(0,1),lty="dotted")
dev.off()

plot(fit$x[1,],pc[,1])

fit$rotation[1:3,1:3]
fit1$loadings[1:3,1:3]
k=6
plot(fit$rotation[,k],-fit1$loadings[,k]); abline(c(0,1))
plot(fit1$loadings[,k],fit2$svd$V[,k]); abline(c(0,1))
plot(fit$x[,k],-fit2$ind$coord[,k]); abline(c(0,1))
plot(fit1$scores[,k],fit2$ind$coord[,k]); abline(c(0,1))


####################################################################
####################################################################
## k-means clustering

## ---------------------------------------------
typeList=c("",sort(unique(ann$type)))
typeList=""

centrFlag="_noCentering"
centrFlag=""

scaleList=c("_noScaling","")
scaleList=""

cohort="_mycRasWt"
cohort="_mycRas"
cohort="_wt"
cohortList=c("_mycRasWt","_mycRas","_wt")
cohortList=c("_mycRasWt")

for (scaleFlag in scaleList) {
    out1=out2=NULL
    for (cohort in cohortList) {
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
        }
        )
        for (typeFlag in typeList) {
            if (typeFlag=="") {
                header=paste(cohortName,": All samples clustered together",sep="")
            } else {
                header=paste(cohortName,": Samples clustered within ",typeFlag,sep="")
            }
            arrayData=t(cell_rbf)
            annFeat=ann_rbf
            if (typeFlag=="") arrayData2=arrayData else arrayData2=arrayData[which(annFeat$type==typeFlag),]
            if (centrFlag=="") {
                centr=apply(arrayData2,1,median,na.rm=T)
                for (i in 1:nrow(arrayData2)) {
                    arrayData2[i,]=arrayData2[i,]-centr[i]
                }
            }
            if (scaleFlag=="") {
                scal=apply(arrayData2,1,sd,na.rm=T)
                for (i in 1:nrow(arrayData2)) {
                    arrayData2[i,]=arrayData2[i,]/scal[i]
                }
            }
            for (nClust in 2:15) {
                fit=kmeans(arrayData2, nClust)
                out=c("feature",cohort,typeFlag,nClust,fit$tot.withinss,fit$betweenss)
                out1=rbind(out1,out)
                fit=kmeans(t(arrayData2), nClust)
                out=c("sample",cohort,typeFlag,nClust,fit$tot.withinss,fit$betweenss)
                out2=rbind(out2,out)
                if (F) {
                    png(paste("kmeans_",nClust,"cluster_samples_",typeFlag,cohort,centrFlag,scaleFlag,".png",sep=""))
                    plot(arrayData2, col = fit$cluster)
                    points(fit$centers, col = 1:2, pch = 8, cex = 2)
                    plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
                    dev.off()
                }
            }
        }
    }
    colnames(out1)=colnames(out2)=c("feature","cohort","type","nClust","totWSS","BSS")
    kmeanF=as.data.frame(out1,stringsAsFactors=T)
    kmeanS=as.data.frame(out2,stringsAsFactors=T)
    for (k in c("nClust","totWSS","BSS")) {
        kmeanF[,k]=as.numeric(as.character(kmeanF[,k]))
        kmeanS[,k]=as.numeric(as.character(kmeanS[,k]))
    }
    rownames(kmeanF)=rownames(kmeanS)=NULL
    kmeanF$BSSbyWSS=kmeanF$BSS/kmeanF$totWSS
    kmeanS$BSSbyWSS=kmeanS$BSS/kmeanS$totWSS
    kmeanF[which(kmeanF$BSSbyWSS==max(kmeanF$BSSbyWSS)),]
    kmeanS[which(kmeanS$BSSbyWSS==max(kmeanS$BSSbyWSS)),]

    for (featureFlag in c("feature","sample")) {
        switch(featureFlag,
        "feature"={
            tbl=kmeanF
        },
        "sample"={
            tbl=kmeanS
        }
        )
        fName=featureFlag
        heading=featureFlag
        lim=range(tbl$totWSS)
        grp=paste(tbl$cohort,tbl$type,sep="/")
        grpUniq=unique(grp)
        for (gId in 1:length(grpUniq)) {
            k=which(grp==grpUniq[gId])
            png(paste("screeplot_kmeans_",fName,"_",sub("/","_",grpUniq[gId]),centrFlag,scaleFlag,".png",sep=""))
            plot(tbl$nClust[k], tbl$totWSS[k], ylim=lim, type="b", main=paste("K-means clustering of ",heading,"s: ",grpUniq[gId],sep=""), xlab="Number of Clusters",ylab="Within groups sum of squares")
            dev.off()
        }
    }
}

#save.image("tmp.RData")

## ---------------------------------------------
typeList=c("",sort(unique(ann$type)))
typeList=""

centrFlag=""

scaleList=c("_noScaling","")
scaleList=""

cohort="_mycRas"
cohort="_wt"
cohort="_mycRasWt"
cohortList=c("_mycRasWt","_mycRas","_wt")
cohortList=c("_mycRasWt")

for (scaleFlag in scaleList) {
    out1=out2=NULL
    for (cohort in cohortList) {
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
        }
        )
        for (typeFlag in typeList) {
            if (typeFlag=="") {
                header=paste(cohortName,": All samples clustered together",sep="")
            } else {
                header=paste(cohortName,": Samples clustered within ",typeFlag,sep="")
            }
            arrayData=t(cell_rbf)
            annFeat=ann_rbf
            if (typeFlag=="") arrayData2=arrayData else arrayData2=arrayData[which(annFeat$type==typeFlag),]
            if (centrFlag=="") {
                centr=apply(arrayData2,1,median,na.rm=T)
                for (i in 1:nrow(arrayData2)) {
                    arrayData2[i,]=arrayData2[i,]-centr[i]
                }
            }
            if (scaleFlag=="") {
                scal=apply(arrayData2,1,sd,na.rm=T)
                for (i in 1:nrow(arrayData2)) {
                    arrayData2[i,]=arrayData2[i,]/scal[i]
                }
            }

            nClustList=2:15
            nPerm=10
            nPerm=1000
            nPerm=10000

            ## ------------------------------
            ## Takes 1.04 min for 10 perms for k=2:15
            if (F) {
                out=matrix(nrow=length(nClustList),ncol=nPerm)
                rownames(out)=nClustList
                cat("No. of clusters: ",nrow(out),"\n")
                timeStamp=Sys.time()
                print(format(timeStamp, "%x %X"))
                for (k in 1:length(nClustList)) {
                    cat("Cluster: ",nClustList[k],"\n")
                    res=sapply(1:nPerm,function(x) {
                        fit=kmeans(t(arrayData2),nClustList[k])
                        fit$tot.withinss
                    },USE.NAMES=F)
                    out[k,]=res
                }
                timeStamp=c(timeStamp,Sys.time())
                print(format(timeStamp[2], "%x %X"))
                print(diff(timeStamp))
                save(out,file=paste("kmeans_",nPerm,"perms.RData",sep=""))
                res=apply(out,1,min,na.rm=T)
                sort(res,decreasing=T)
                png(paste("screeplot_kmeans",cohort,"_",nPerm,"perms.png",sep=""))
                plot((1:length(res))+1, res, ylim=lim, type="b", main=paste("K-means clustering: ",cohortName,sep=""), xlab="Number of Clusters",ylab="Within groups sum of squares")
                dev.off()
            }
            
            ## ------------------------------
            timeStamp=Sys.time()
            print(format(timeStamp, "%x %X"))
            for (k in 1:length(nClustList)) {
                cat("Cluster: ",nClustList[k],"\n")
                fit=kmeans(t(arrayData2),nClustList[k],nstart=nPerm)
                save(fit,file=paste("kmeans_",nClustList[k],"clusts_",nPerm,"starts.RData",sep=""))
            }
            timeStamp=c(timeStamp,Sys.time())
            print(format(timeStamp[2], "%x %X"))
            print(diff(timeStamp))
            
            res=rep(NA,length(nClustList))
            names(res)=nClustList
            for (k in 1:length(nClustList)) {
                load(file=paste("kmeans_",nClustList[k],"clusts_",nPerm,"starts.RData",sep=""))
                res[k]=fit$tot.withinss

            }
            sort(res,decreasing=T)
            png(paste("screeplot_kmeans",cohort,"_",nPerm,"starts.png",sep=""))
            plot((1:length(res))+1, res, ylim=lim, type="b", main=paste("K-means clustering: ",cohortName,sep=""), xlab="Number of Clusters",ylab="Within groups sum of squares")
            dev.off()

        }
    }
    colnames(out1)=colnames(out2)=c("feature","cohort","type","nClust","totWSS","BSS")
    kmeanF=as.data.frame(out1,stringsAsFactors=T)
    kmeanS=as.data.frame(out2,stringsAsFactors=T)
    for (k in c("nClust","totWSS","BSS")) {
        kmeanF[,k]=as.numeric(as.character(kmeanF[,k]))
        kmeanS[,k]=as.numeric(as.character(kmeanS[,k]))
    }
    rownames(kmeanF)=rownames(kmeanS)=NULL
    kmeanF$BSSbyWSS=kmeanF$BSS/kmeanF$totWSS
    kmeanS$BSSbyWSS=kmeanS$BSS/kmeanS$totWSS
    kmeanF[which(kmeanF$BSSbyWSS==max(kmeanF$BSSbyWSS)),]
    kmeanS[which(kmeanS$BSSbyWSS==max(kmeanS$BSSbyWSS)),]
    
    for (featureFlag in c("feature","sample")) {
        switch(featureFlag,
        "feature"={
            tbl=kmeanF
        },
        "sample"={
            tbl=kmeanS
        }
        )
        fName=featureFlag
        heading=featureFlag
        lim=range(tbl$totWSS)
        grp=paste(tbl$cohort,tbl$type,sep="/")
        grpUniq=unique(grp)
        for (gId in 1:length(grpUniq)) {
            k=which(grp==grpUniq[gId])
            png(paste("screeplot_kmeans_",fName,"_",sub("/","_",grpUniq[gId]),centrFlag,scaleFlag,".png",sep=""))
            plot(tbl$nClust[k], tbl$totWSS[k], ylim=lim, type="b", main=paste("K-means clustering of ",heading,"s: ",grpUniq[gId],sep=""), xlab="Number of Clusters",ylab="Within groups sum of squares")
            dev.off()
        }
    }
}
load("/Users/royr/UCSF/WallaceMarshall/kmeans_2clusts_10000starts.RData")
fit2=fit
load("/Users/royr/UCSF/WallaceMarshall/kmeans_3clusts_10000starts.RData")
fit3=fit
load("/Users/royr/UCSF/WallaceMarshall/kmeans_4clusts_10000starts.RData")
fit4=fit
table(names(fit2$cluster)==names(fit3$cluster))
table(names(fit2$cluster)==names(fit4$cluster))
table(clust4=fit4$cluster,clust3=fit3$cluster,clust2=fit2$cluster)
"
, , clust2 = 1

        clust3
clust4    1    2    3
1    0    2    0
2    0  835    0
3    0    0    0
4    0    0 2226

, , clust2 = 2

        clust3
clust4    1    2    3
1    3    0    1
2    0    0    0
3  985    0    0
4    1    0  824
"

## ---------------------------------------------
## NOT USED
## For testing
# Dataset for Clustering
n = 150
g = 6
set.seed(g)
d <- data.frame(x = unlist(lapply(1:g, function(i) rnorm(n/g, runif(1)*i^2))),
y = unlist(lapply(1:g, function(i) rnorm(n/g, runif(1)*i^2))))
mydata<-d
#Plot 3X2 plots
attach(mtcars)
par(mfrow=c(3,2))

#Plot the original dataset
plot(mydata$x,mydata$y,main="Original Dataset")

#Scree plot to deterine the number of clusters
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) {
    wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
}
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

# Ward Hierarchical Clustering
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=5, border="red")

#Silhouette analysis for determining the number of clusters
library(fpc)
asw <- numeric(20)
for (k in 2:20)
asw[[k]] <- pam(mydata, k) $ silinfo $ avg.width
k.best <- which.max(asw)

cat("silhouette-optimal number of clusters:", k.best, "\n")
plot(pam(d, k.best))

# K-Means Cluster Analysis
fit <- kmeans(mydata,k.best)
mydata
# get cluster means
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, clusterid=fit$cluster)
plot(mydata$x,mydata$y, col = fit$cluster, main="K-means Clustering results")

Hope it helps!!
shareimprove this answer

answered Aug 5 '13 at 11:07
Udeep Shakya
27924

add a comment
up vote
5
down vote


For k-mean you might want to have a look at gap statistic

http://blog.echen.me/2011/03/19/counting-clusters/
shareimprove this answer

answe

## ---------------------------------------------
require(graphics)

# a 2-dimensional example
x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
colnames(x) <- c("x", "y")
x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2).
matrix(rnorm(100, mean = 2, sd = 0.3), ncol = 2))
colnames(x) <- c("x", "y","z")
(cl <- kmeans(x, 2))
plot(x, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex = 2)



####################################################################
####################################################################

if (F) {

fileList=dir(datadir2[2],pattern="clusterInfoSample")
print(dim(phen2))
for (fId in 1:length(fileList)) {
    print(fileList[fId])
    clustInfo=read.table(paste(datadir2[2],fileList[fId],sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    j=match(clustInfo$id,phen2$id); j2=which(!is.na(j)); j1=j[j2]
    print(table(is.na(j)))
    print(dim(clustInfo))
    clustInfo=cbind(clustInfo[j2,],class=phen2$class[j1])
    clustInfo=clustInfo[,c("id","class","clustId","order")]
    write.table(clustInfo,file=paste(fileList[fId],sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)
}

}

####################################################################
####################################################################

library(rpart)

cohort="_wt"

cvFlag="_2foldCV"
cvFlag=""

nPerm=2

typeflag=""
typeList=""
typeList=c("",sort(unique(ann$type)))

rpartCtrlFlag="_cp0.05"
rpartCtrlFlag="_minbct10"
rpartCtrlList=c("","_minbct10","_cp0.05")
rpartCtrlFlag=""
rpartCtrlList=c("_minbct5","_minbct10","_minbct20")
rpartCtrlList=c("")

write.table(paste("feature","type",sep="\t"),file=paste("rpart.txt",sep=""),append=F,col.names=F,row.names=F,sep="\t",quote=F)

classPredMat=NULL; nm=c()
x=paste("_wt_",typeList,sep=""); x[x=="_wt_"]="_wt"
classPredMat=matrix("",nrow=length(typeList),ncol=nrow(cellW),dimnames=list(x,rownames(cellW)))
for (rpartCtrlFlag in rpartCtrlList) {
    for (typeFlag in typeList) {
        fNameOut1=paste(cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),sep="")
        if (typeFlag=="") {
            header1=paste(cohortName,", all Features",sep="")
        } else {
            header1=paste(cohortName,", ",typeFlag," features",sep="")
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
            }
        )

        arrayData=t(cell_rbf)
        i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
        arrayData[i2,]=t(cell[,i1])
        annFeat=ann_rbf

        centr=apply(arrayData,1,median,na.rm=T)
        for (i in 1:nrow(arrayData)) {
            arrayData[i,]=arrayData[i,]-centr[i]
        }

        scal=apply(arrayData,1,sd,na.rm=T)
        for (i in 1:nrow(arrayData)) {
            arrayData[i,]=arrayData[i,]/scal[i]
        }

        ## ---------------

        #j=match(colnames(arrayData),annCell$id); j1=which(!is.na(j)); j2=j[j1]
        j2=which(annCell$set==1)
        j=match(colnames(arrayData),annCell$id[j2]); j1=which(!is.na(j)); j2=j2[j[j1]]
        if (typeFlag=="") arrayData2=arrayData else arrayData2=arrayData[which(annFeat$type==typeFlag),]
        
        for (iPerm in 1:ifelse(cvFlag=="",1,nPerm)) {
            if (cvFlag=="") {
                fNameOut=fNameOut1
                header=header1
                samId=1:length(j1)
                kFold=1
                #j=1:length(j1)
            } else {
                samId=sample(1:length(j1),size=length(j1),replace=F)
                kFold=as.integer(gsub("_|foldCV","",cvFlag))
                #j=sample(1:length(j1),size=round(length(j1)/2),replace=F)
            }
            for (thisFold in 1:kFold) {
                if (cvFlag=="") {
                    j=samId
                } else {
                    fNameOut=paste(fNameOut1,cvFlag,"_",thisFold,"_",iPerm,"iter",sep="")
                    header=paste(header1,", ",thisFold," k, ",iPerm," iteration",sep="")
                    jMat=matrix(rep(round(seq(1,length(samId)+1,length=kFold+1)),2),ncol=2,byrow=F)
                    jMat[1:(nrow(jMat)-1),2]=jMat[2:nrow(jMat),1]-1
                    jMat=jMat[-nrow(jMat),]
                    if (!is.matrix(jMat)) jMat=t(as.matrix(jMat))
                    j=samId[-(jMat[thisFold,1]:jMat[thisFold,2])]
                }
                #dat=cbind(class=annCell$class[j2],as.data.frame(t(arrayData[1:10,j1])))
                #dat=cbind(class=annCell$class[j2[j]],as.data.frame(t(arrayData2[,j1[j]])),stringsAsFactors=F)
                dat=cbind(class=paste("class",annCell$class[j2[j]],sep=""),as.data.frame(t(arrayData2[,j1[j]])),stringsAsFactors=F)
                dat2=dat

                if (F) {
                    fit=rpart(class ~ ., method="class",data = dat)
                    #png(paste("rpart",fNameOut,".png",sep=""))
                    png(paste("rpart",fNameOut,".png",sep=""),width=1.5*480,height=480)
                    #par(mfrow = c(1,2), xpd = NA) # otherwise on some devices the text is clipped
                    par(xpd = NA) # otherwise on some devices the text is clipped
                    plot(fit,main=header)
                    text(fit, use.n = TRUE)
                    dev.off()
                }

                if (length(grep("_minbct",rpartCtrlFlag))==1) {
                    fit=rpart(class ~ ., control=rpart.control(minbucket=as.integer(sub("_minbct","",rpartCtrlFlag))), data = dat)
                } else if (length(grep("_cp",rpartCtrlFlag))==1) {
                    fit=rpart(class ~ ., control=rpart.control(cp=as.numeric(sub("_cp","",rpartCtrlFlag))), data = dat)
                } else {
                    fit=rpart(class ~ ., data = dat)
                }
                
                x=unique(as.character(fit$frame$var)); x=x[x!="<leaf>"]
                write.table(cbind(feature=x,type=typeFlag),file=paste("rpart.txt",sep=""),append=T,col.names=F,row.names=F,sep="\t",quote=F)
                
                fit1=fit
                png(paste("rpart",fNameOut,rpartCtrlFlag,".png",sep=""),width=2.5*480,height=480)
                par(xpd = NA) # otherwise on some devices the text is clipped
                plot(fit)
                text(fit, use.n = TRUE,cex=.9)
                dev.off()
                
                datadir="rpart/"
                datadir=""
                sink(paste(datadir,"rpart",fNameOut,rpartCtrlFlag,".txt",sep=""),append=F)
                print(fit1)
                sink()

                ## ----------------------------------
                ## Get the prediction criteria

                datadir="rpart/"
                datadir=""

                fit=fit1; tbl=read.table(paste(datadir,"rpart",fNameOut,rpartCtrlFlag,".txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
                x=sapply(tbl[5:nrow(tbl),1],function(x) {
                    offset=2
                    z=strsplit(strsplit(x,")")[[1]][1],"")[[1]]
                    y=gsub(" +"," ",sub(" +","",sub(" +$", "", x)))
                    for (gId in 1:length(grpUniq)) {
                        y=sub(grpUniq[gId],sub(" ","",grpUniq[gId]),y)
                    }
                    y=strsplit(y," ")[[1]]
                },USE.NAMES=F)


                xx=tbl[5:nrow(tbl),1]
                pp=c()
                ppLeaf=c()
                grpLeaf=c()
                kk=1:16
                kk=1:length(xx)
                #for (k in 1:length(xx)) {
                nodeFlag=c()
                nodePrevFlag=0
                for (k in kk) {
                    x=xx[k]
                    #pId=pPrevId+1
                    offset=2
                    z=strsplit(strsplit(x,")")[[1]][1],"")[[1]]
                    pId=length(grep(" ",z))-offset
                    z=nchar(strsplit(x,")")[[1]][1])+1
                    nodeThisFlag=z
                    
                    y=gsub(" +"," ",sub(" +","",sub(" +$", "", x)))
                    y=sub("> ",">",sub("< ","<",y))
                    for (gId in 1:length(grpUniq)) {
                        y=sub(grpUniq[gId],sub(" ","",grpUniq[gId]),y)
                    }
                    y=strsplit(y," ")[[1]]
                    k2=grep("-",y[2])
                    if (length(k2)==1) {
                        y[2]=paste(sub("-","(-",y[2]),")",sep="")
                    }
                    #ppThis=paste(y[2:3],collapse="")
                    ppThis=paste("dat$",y[2],sep="")
                    
                    if (nodePrevFlag<nodeThisFlag) {
                        nodeFlag=c(nodeFlag,z)
                        pp=c(pp,ppThis)
                    } else {
                        for (k2 in length(nodeFlag):0) {
                            if (k2==0) {
                                pp=c()
                                nodeFlag=c()
                                break
                            }
                            if (nodeFlag[k2]<nodeThisFlag) {
                                pp=pp[1:k2]
                                nodeFlag=nodeFlag[1:k2]
                                break
                            }
                        }
                        nodeFlag=c(nodeFlag,z)
                        pp=c(pp,ppThis)
                    }
                    if (y[length(y)]=="*") {
                        ppLeaf=c(ppLeaf,paste(pp,collapse=" & "))
                        grpLeaf=c(grpLeaf,y[5])
                    }
                    nodePrevFlag=nodeThisFlag
                }

                ## ----------------------------------
                ## Check observed vs. predicted class

                if (T) {
                    cat("\n\n=================== ",paste("rpart",fNameOut,rpartCtrlFlag,"\n\n",sep="")," ==================",sep="")
                    cat("\n\n=================== Set1 (training set) ==================\n\n",sep="")
                    dat=dat2
                    classObs=sub(" ","",dat$class)
                    classObs=dat$class
                    if (class(dat$class)%in%c("integer","numeric")) {
                        classPred=rep(NA,nrow(dat))
                    } else {
                        classPred=rep("",nrow(dat))
                    }
                    for (k in 1:length(ppLeaf)) {
                        a=paste("classPred[which(",ppLeaf[k],")]=grpLeaf[k]",sep="")
                        eval(parse(text=a))
                    }
                    if (class(dat$class)%in%c("integer","numeric")) {
                        classPred=as.numeric(classPred)
                        cat("\nCor: ",round(cor(classObs,classPred,use="complete.obs"),2),"\n")
                    } else {
                        print(table(classObs,classPred))
                        cat("\nMisclassification rate: ",round(mean(classObs!=classPred),2),sep="")
                    }
                    classPred2=classPred
                    if (cvFlag=="") {
                        #out=rep("",nrow(cell))
                        #out[match(rownames(dat),rownames(cell))]=classPred
                        #classPredMat=rbind(classPredMat,out)
                        classPredMat[which(rownames(classPredMat)==paste(fNameOut,rpartCtrlFlag,sep="")),match(rownames(dat),colnames(classPredMat))]=classPred
                        nm=c(nm,paste(fNameOut,rpartCtrlFlag,sep=""))
                    }
                    
                    ## -----------------------
                    if (2%in%annCell$set | cvFlag!="") {
                        if (2%in%annCell$set) {
                            j=which(annCell$set==2)
                            #dat=cbind(class=annCell$class[j],as.data.frame(t(arrayData2[,match(annCell$id[j],colnames(arrayData2))])),stringsAsFactors=F)
                            dat=cbind(class=paste("class",annCell$class[j],sep=""),as.data.frame(t(arrayData2[,match(annCell$id[j],colnames(arrayData2))])),stringsAsFactors=F)
                        } else {
                            j=which(!annCell$id[j2]%in%rownames(dat2))
                            #dat=cbind(class=annCell$class[j2[j]],as.data.frame(t(arrayData2[,j1[j]])),stringsAsFactors=F)
                            dat=cbind(class=paste("class",annCell$class[j2[j]],sep=""),as.data.frame(t(arrayData2[,j1[j]])),stringsAsFactors=F)
                        }
                        dat1=dat
                        
                        cat("\n\n=================== Set2 (validation set) ==================\n\n",sep="")
                        dat=dat1
                        classObs=sub(" ","",dat$class)
                        classObs=dat$class
                        if (class(dat$class)%in%c("integer","numeric")) {
                            classPred=rep(NA,nrow(dat))
                        } else {
                            classPred=rep("",nrow(dat))
                        }
                        for (k in 1:length(ppLeaf)) {
                            a=paste("classPred[which(",ppLeaf[k],")]=grpLeaf[k]",sep="")
                            eval(parse(text=a))
                        }
                        if (class(dat$class)%in%c("integer","numeric")) {
                            classPred=as.numeric(classPred)
                            cat("\nCor: ",round(cor(classObs,classPred,use="complete.obs"),2),"\n")
                        } else {
                            print(table(classObs,classPred))
                            cat("\nMisclassification rate: ",round(mean(classObs!=classPred),2),sep="")
                        }
                        classPred1=classPred
                        #out=rep("",nrow(cell))
                        #out[match(rownames(dat),rownames(cell))]=classPred
                        #classPredMat=rbind(classPredMat,out)
                        classPredMat[which(rownames(classPredMat)==paste(fNameOut,rpartCtrlFlag,sep="")),match(rownames(dat),colnames(classPredMat))]=classPred
                        nm=c(nm,paste(fNameOut,rpartCtrlFlag,sep=""))
                        #classPredMat[nrow(classPredMat),match(rownames(dat),rownames(cell))]=classPred
                    }
                }
            }
        }
    }
}
#colnames(classPredMat)=rownames(cell)
#rownames(classPredMat)=nm

table(classPredMat[1,]=="")

dat=dat2
#dat=cbind(class=annCell$class[j2],as.data.frame(t(arrayData2[,j1])),stringsAsFactors=F)
dat=cbind(class=paste("class",annCell$class[j2],sep=""),as.data.frame(t(arrayData2[,j1])),stringsAsFactors=F)
x=apply(classPredMat[,match(rownames(dat),colnames(classPredMat))],1,function(x) {
    mean(x!=dat$class,na.rm=T)
})
round(x,2)
round(x[1:4],2)

tbl=t(classPredMat)
colnames(tbl)=paste("classPred",colnames(tbl),sep="")
tbl[tbl==""]=NA
classObs=rep("",nrow(tbl))
classObs[classObs==""]=NA
j=match(rownames(tbl),annCell$id); j1=which(!is.na(j)); j2=j[j1]
#classObs[j1]=annCell$class[j2]
classObs[j1]=paste("class",annCell$class[j2],sep="")
tbl=cbind(id=rownames(tbl),classObs,as.data.frame(tbl))
write.table(tbl,file=paste("predictedClass_rpart_wt.txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)



"
round(x,2)
_wt_minbct5          _wt_cell_minbct5  _wt_mitochondria_minbct5
0.05                      0.05                      0.12
_wt_nucleus_minbct5              _wt_minbct10         _wt_cell_minbct10
0.21                      0.05                      0.05
_wt_mitochondria_minbct10      _wt_nucleus_minbct10              _wt_minbct20
0.30                      0.33                      0.28
_wt_cell_minbct20 _wt_mitochondria_minbct20      _wt_nucleus_minbct20
0.28                      0.35                      0.45


round(x[1:4],2)
_wt         _wt_cell _wt_mitochondria      _wt_nucleus
0.05             0.05             0.14             0.26
"
## all features, cells perform best
## min bucket, cp do not make a difference for these 2 subsets

####################################################################
####################################################################
## Distance of each sample to class centroid

cohort="_wt"
cohort="_mycRasWt"
typeFlag=""

fNameOut1=paste(cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),sep="")
if (typeFlag=="") {
    header1=paste(cohortName,", all Features",sep="")
} else {
    header1=paste(cohortName,", ",typeFlag," features",sep="")
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
}
)

arrayData=t(cell_rbf)
i=match(colnames(cell),rownames(arrayData)); i1=which(!is.na(i)); i2=i[i1]
arrayData[i2,]=t(cell[,i1])
annFeat=ann_rbf

centr=apply(arrayData,1,median,na.rm=T)
for (i in 1:nrow(arrayData)) {
    arrayData[i,]=arrayData[i,]-centr[i]
}

scal=apply(arrayData,1,sd,na.rm=T)
for (i in 1:nrow(arrayData)) {
    arrayData[i,]=arrayData[i,]/scal[i]
}

## ---------------

j=match(colnames(arrayData),annCell$id); j1=which(!is.na(j)); j2=j[j1]
if (typeFlag=="") arrayData2=arrayData else arrayData2=arrayData[which(annFeat$type==typeFlag),]

getDistance=function(x,y,method="euclidean") {
    if (method%in%c("pearson","spearman","kendall")) {
        z=cor(x,y,use="complete.obs",method=method)
    } else {
        switch(method,
        "euclidean"={
            z=sqrt(sum((x-y)^2,na.rm=T))
        }
        )
    }
    z
}

tmpC=rep("",ncol(arrayData))
dat=data.frame(id=colnames(arrayData),class=tmpC,stringsAsFactors=F)
#dat$class[j1]=annCell$class[j2]
#grpUniq=sort(unique(annCell$class))
dat$class[j1]=paste("class",annCell$class[j2],sep="")
grpUniq=sort(unique(paste("class",annCell$class,sep="")))
classDistMatE=classDistMatP=classDistMatS=classDistMatK=matrix(nrow=nrow(dat),ncol=length(grpUniq),dimnames=list(dat$id,grpUniq))
jj=which(dat$class=="")
jj=1:nrow(dat)
for (gId in 1:length(grpUniq)) {
    j=which(dat$class==grpUniq[gId])
    centr=apply(arrayData2[,j],1,mean,na.rm=T)
    x=apply(arrayData2[,jj],2,getDistance,method="euclidean",y=centr)
    classDistMatE[jj,gId]=x
    x=apply(arrayData2[,jj],2,getDistance,method="pearson",y=centr)
    classDistMatP[jj,gId]=1-abs(x)
    centr=apply(arrayData2[,j],1,median,na.rm=T)
    x=apply(arrayData2[,jj],2,getDistance,method="spearman",y=centr)
    classDistMatS[jj,gId]=1-abs(x)
    x=apply(arrayData2[,jj],2,getDistance,method="kendall",y=centr)
    classDistMatK[jj,gId]=1-abs(x)
}


i=2
centr=apply(arrayData2[,j],1,median,na.rm=T)
x2=apply(arrayData2[,jj],2,getDistance,method="spearman",y=centr)
x=arrayData2[,jj][,i]; y=centr
z1=cor(x,y,use="complete.obs",method="spearman")
z2=sqrt(sum((x-y)^2,na.rm=T))
c(z1,z2)

plot(classDistMatE[,1], abs(classDistMatP[,1]))

j=match(dat$id,rownames(classDistMatE)); j1=which(!is.na(j)); j2=j[j1]
classDistMat=classDistMatE
tmp=matrix(nrow=nrow(dat),ncol=ncol(classDistMat),dimnames=list(dat$id,paste("dist2",colnames(classDistMat),sep="")))
tmp[j1,]=classDistMat[j2,]
tmp[j1,]=t(apply(classDistMat[j2,],1,function(x) {
    y=order(order(x))
    y[is.na(x)]=NA
    y
}))
tmpE=tmp
classDistMat=1-abs(classDistMatP)
tmp=matrix(nrow=nrow(dat),ncol=ncol(classDistMat),dimnames=list(dat$id,paste("dist2",colnames(classDistMat),sep="")))
tmp[j1,]=classDistMat[j2,]
tmp[j1,]=t(apply(classDistMat[j2,],1,function(x) {
    y=order(order(x))
    y[is.na(x)]=NA
    y
}))
tmpP=tmp
classDistMat=1-abs(classDistMatS)
tmp=matrix(nrow=nrow(dat),ncol=ncol(classDistMat),dimnames=list(dat$id,paste("dist2",colnames(classDistMat),sep="")))
tmp[j1,]=classDistMat[j2,]
tmp[j1,]=t(apply(classDistMat[j2,],1,function(x) {
    y=order(order(x))
    y[is.na(x)]=NA
    y
}))
tmpS=tmp
classDistMat=1-abs(classDistMatK)
tmp=matrix(nrow=nrow(dat),ncol=ncol(classDistMat),dimnames=list(dat$id,paste("dist2",colnames(classDistMat),sep="")))
tmp[j1,]=classDistMat[j2,]
tmp[j1,]=t(apply(classDistMat[j2,],1,function(x) {
    y=order(order(x))
    y[is.na(x)]=NA
    y
}))
tmpK=tmp
for (grpThis in sort(unique(dat$class))) {
    cat("\n\n-------- ",grpThis,"\n",sep="")
    j=which(dat$class==grpThis)
    k=which(colnames(tmpE)==paste("dist2",grpThis,sep=""))
    print(table(tmpE[j,k]))
    print(table(tmpP[j,k]))
    print(table(tmpS[j,k]))
    print(table(tmpK[j,k]))
}

save(classPredMat,classDistMatE,classDistMatP,classDistMatS,classDistMatK,file="class.RData")



####################################################################
####################################################################
