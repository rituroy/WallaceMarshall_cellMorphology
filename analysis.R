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

## ----------------------------------------------
datadir="docs/"
cellW2=read.table(paste(datadir,"20170707_WT_RowAB_CONCATENATED.csv",sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T)
tmpC=rep("",ncol(cellW2))
annW2=data.frame(feature=paste("feat_",1:ncol(cellW2),sep=""),type=tmpC,stringsAsFactors=F)
annCellW2=data.frame(id=paste("w2_",1:nrow(cellW2),sep=""),stringsAsFactors=F)

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

####################################################################
####################################################################

cell=cellM

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
	corMatSam=cor(t(cell),use="complete.obs",method="kendall")
	save(corMatSam,file="corMatSam.RData")
	
	corMat=cor(cell,use="complete.obs",method="kendall")
	save(corMat,file="corMatSam.RData")

	summary(abs(c(corMat[lower.tri(corMat)])))
	#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	#0.00000 0.01226 0.02827 0.07226 0.06898 1.00000 

	table(abs(c(corMat[lower.tri(corMat)]))==1)
	corMat[lower.tri(corMat)][abs(c(corMat[lower.tri(corMat)]))==1]
	corMat[abs(corMat)==1]

	n=ncol(cell)
	tmp=vector(mode="numeric",length=n*(n-1)/2)
	tmpC=vector(mode="character",length=n*(n-1)/2)
	corInfo=data.frame(feature1=tmpC,feature2=tmpC,corKend=tmp,stringsAsFactors=F)
	k=1
	for (k1 in 1:(nrow(corMat)-1)) {
		for (k2 in (k1+1):nrow(corMat)) {
			corInfo$feature1[k]=rownames(corMat)[k1]
			corInfo$feature2[k]=rownames(corMat)[k2]
			corInfo$corKend[k]=corMat[k1,k2]
			k=k+1
		}
	}
	write.table(corInfo,file=paste("corrFeature",cohort,".txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)
}

datadir="results/"
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
## NOT USED

#load("tmp.RData")
#load("tmp_w3830.RData"); datVer="w3830"
load("tmp_w5969.RData"); datVer="w5969"


annRpart=read.table(paste("rpart.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

library(marray)
#library(fpc)
#source(paste(dirSrc,"functions/heatmap.5.R",sep=""))
#source(paste(dirSrc,"functions/heatmapAcgh.7.R",sep=""))
#source(paste(dirSrc,"functions/heatmap.5.2.R",sep=""))
source(paste(dirSrc,"functions/heatmap.5.4.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.1.R",sep=""))

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

datadir=""
if (datVer=="w3830") {
    datadir="results/w3830/cor/"
} else {
    datadir="results/cor/"
}

typeFlag="cell"

sepFClustFlag=T ## NOT USED for typeFlag=""
sepFClustFlag=F

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
typeList=sort(unique(ann$type))
typeList=c("",sort(unique(ann$type)))
typeList=""

sampleBar=""
sampleBar="cluster"

distMethod="pearson"
distMethod="spearman"
distMethod="kendall"

classDistFlag="euclidean"
classDistFlag="kendall"

linkMethod="ward.D2"

centrFlag="_noCentering"
centrFlag=""

scaleList=c("","_noScaling")
scaleList=c("")
scaleFlag="_noScaling"
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

getCorFlag=T
getCorFlag=F

orderFlag=rep("",2)
datadir2=rep("",2)

orderFlag[1]="_mycRasOrd"
orderFlag[1]="_wt_classSet1Ord"
orderFlag[1]=""
orderFlag[1]="_wtOrd"
datadir2[1]=""
if (datVer=="w3830") {
    datadir2[1]=paste("results/w3830/heatmap/",sub("_","",sub("Ord","",orderFlag[1])),"/",sub("_","",orderFlag[1]),"/",sub("_","",type2Flag),"/",sep="")
} else {
    datadir2[1]=paste("results/heatmap/",sub("_","",sub("Ord","",orderFlag[1])),"/",sub("_","",orderFlag[1]),"/",sub("_","",type2Flag),"/",sep="")
}

orderFlag[2]="_kmeans2ClustOrd"
if (datVer=="w3830") {
    datadir2[2]=paste("results/w3830/heatmap/",sub("_","",cohort),"/kmeans/",sep="")
} else {
    datadir2[2]=paste("results/heatmap/",sub("_","",cohort),"/kmeans/",sep="")
}
orderFlag[2]=""
datadir2[2]=""

orderSamList=c("_hclust4ClustOrd","_kmeans2ClustOrd","_kmeans3ClustOrd","_kmeans4ClustOrd")
orderSamList=c("_hclust4ClustOrd")
orderSamList=c("_kmeans2ClustOrd","_kmeans3ClustOrd","_kmeans4ClustOrd")
orderSamList=c("_kmeans3ClustOrd")
orderSamList=c("")

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
            #header=paste(header,", classDist-",classDistFlag,sep="")
            load("results/class.RData")
            switch(classDistFlag,
            "euclidean"={classDistMat=classDistMatE},
            "pearson"={classDistMat=classDistMatP},
            "spearman"={classDistMat=classDistMatS},
            "kendall"={classDistMat=classDistMatK}
            )
            #phen2=data.frame(id=paste("sam",1:nrow(cell),sep=""),sd=sdSam,stringsAsFactors=F)
            phen2=data.frame(id=rownames(cell),sd=sdSam,stringsAsFactors=F)
            j=match(phen2$id,annCell$id); j1=which(!is.na(j)); j2=j[j1]
            if (length(j1)!=0) {
                #phen2$class=""
                #phen2$class[j1]=annCell$class[j2]
                phen2$set1Class=""
                j2=which(annCell$set==1)
                j=match(phen2$id,annCell$id[j2]); j1=which(!is.na(j)); j2=j2[j[j1]]
                #phen2$set1Class[j1]=sub("class","",annCell$class[j2])
                phen2$set1Class[j1]=sub("class","",paste("class",annCell$class[j2],sep=""))
                phen2$set2Class=""
                j2=which(annCell$set==2)
                j=match(phen2$id,annCell$id[j2]); j1=which(!is.na(j)); j2=j2[j[j1]]
                #phen2$set2Class[j1]=sub("class","",annCell$class[j2])
                phen2$set2Class[j1]=sub("class","",paste("class",annCell$class[j2],sep=""))
            }
            j=match(phen2$id,rownames(classDistMat)); j1=which(!is.na(j)); j2=j[j1]
            tmp=matrix(nrow=nrow(phen2),ncol=ncol(classDistMat),dimnames=list(phen2$id,paste("dist2",colnames(classDistMat),sep="")))
            tmp[j1,]=classDistMat[j2,]
            tmp[j1,]=t(apply(classDistMat[j2,],1,function(x) {
                y=order(order(x))
                y[is.na(x)]=NA
                y
            }))
            colnames(tmp)=paste("set1",capWords(colnames(tmp)),sep="")
            phen2=cbind(phen2,tmp)
			
			varFList="featRatio"
			varFName=paste(c("featRatio")," ",sep="")
			varFList=varFName=NULL
			varList=c("sd")
			varName=paste(c("stdDev")," ",sep="")
            
            if (cohort=="_mycRasWt") {
                #phen2=cbind(phen2,cohort=c(rep("Myc/Ras",nrow(cellM)),rep("Wildtype",nrow(cellW))))
                phen2$cohort="Myc/Ras"
                phen2$cohort[substr(phen2$id,1,1)=="w"]="Wildtype"
                varList=c(varList,"cohort","set1Class")
                varName=c(varName,paste(c("cohort","set1Class")," ",sep=""))
            } else if (cohort=="_wt") {
                phen2$set1ClassPred=""
                j=match(phen2$id,colnames(classPredMat)); j1=which(!is.na(j)); j2=j[j1]
                phen2$set1ClassPred[j1]=sub("class","",classPredMat[paste(cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),sep=""),j2])
                varList=c(varList,"set1Class","set1ClassPred")
                varName=c(varName,paste(c("set1Class","set1ClassPred")," ",sep=""))
                if (subsetFlag=="_classSet1Set2") {
                    varList=c(varList,"set2Class")
                    varName=c(varName,paste(c("set2Class")," ",sep=""))
                }
            }
            nm=names(phen2)[grep("dist2class",names(phen2))]
            varList=c(varList,nm)
            varName=c(varName,nm)
			
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
#				featId=featId1[ann$type[featId1]==typeFlag]
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
					   } else {
							limit1=c(-0.5,0.5)
					   }
					   annFeat=ann_rbf
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
					   annFeat=ann_rbfm
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
					   annFeat=ann_rf[i,]
					   varFList=c("type")
					   varFName=paste(c("type")," ",sep="")
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
            
            phenAll=phen2

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
                j=1:nrow(phen2)
                if (subsetFlag=="_classSet1") {
                    j=which(phen2$set1Class!="")
                } else if (subsetFlag=="_classSet1Set2") {
                    j=which(phen2$set1Class!="" | phen2$set2Class!="")
                } else if (subsetFlag=="_classSet1MissClass") {
                    j=which(phen2$set1Class!=phen2$set1ClassPred)
                } else if (subsetFlag=="_classSet1MissClassRpart") {
                    i=match(annFeat$feature,annRpart$feature[which(annRpart$type==typeFlag)]); i=which(!is.na(i))
                    j=which(phen2$set1Class!=phen2$set1ClassPred)
                    nClust[1]=NA
                } else if (subsetFlag=="_classSet1Set2MissClassRpart") {
                    i=match(annFeat$feature,annRpart$feature[which(annRpart$type==typeFlag)]); i=which(!is.na(i))
                    j=which((phen2$set1Class!="" & phen2$set1Class!=phen2$set1ClassPred) | (phen2$set2Class!="" & phen2$set2Class!=phen2$set1ClassPred))
                    nClust[1]=NA
                }
                arrayData=arrayData[i,j]
                phen2=phen2[j,]
                annFeat=annFeat[i,]
            }
            annFeatAll=annFeat

			varFListAll=varFList
			varFNameAll=varFName

			cat("\n\n======================",fNameOut,"======================\n")


			## -------------------
			if (typeFlag=="") arrayData2=arrayData else arrayData2=arrayData[which(annFeat$type==typeFlag),]
            #save(arrayData2,file=paste("arrayData2",cohort,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,subsetFlag,".RData",sep=""))
            save(arrayData2,file=paste("arrayData2",cohort,subsetFlag,ifelse(typeFlag=="","",paste("_",typeFlag,sep="")),type2Flag,".RData",sep=""))
            
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
                    phen2$clust=tbl$clustId
                    phen2=phen2[j,]
                    arrayData=arrayData[,j]
                    nm=paste("clust",gsub("_kmeans|ClustOrd","",orderFlag[2]),sep="")
                    names(phen2)[match("clust",names(phen2))]=nm
                    varList=c(varList,nm)
                    varName=c(varName,paste(c("clust")," ",sep=""))

                }
				clustC=NA
				nClust[2]=NA
			}
            #phenAll=phen2
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
    #			par(mfrow=c(length(type3List),1))
                for (type3Flag in type3List) {
                    cat("===========",type3Flag,"1 =========\n\n")
                    if (sepFClustFlag & type2Flag%in%c("_reducedBioFeatPC","_reducedBioFeatMean")) {
                        fNameOut2=paste(fNameOut,"_",type3Flag,sep="")
                    } else {
                        fNameOut2=fNameOut
                    }
                    fNameOut4=fNameOut2
                    if (typeFlag=="") i=1:nrow(annFeatAll) else i=which(annFeatAll$type==type3Flag)
                    arrayData=arrayData2=arrayDataAll[i,]
                    annFeat=annFeatAll[i,]
    #				if (!sepFClustFlag & type2Flag%in%c("_reducedBioFeatPC","_reducedBioFeatMean")) {
                    if (!sepFClustFlag & typeFlag!="" & type2Flag%in%c("","_reducedBioFeatPC","_reducedBioFeatMean")) {
                        arrayData=arrayDataAll
                        annFeat=annFeatAll
                        x=c()
                        for (type4Flag in sort(unique(annFeatAll$type))) {
                            fNameOut4=paste(fNameOut,"_",type4Flag,sep="")
                            arrayData2=arrayData[which(annFeat$type==type4Flag),]
                            if (orderFlag[1]%in%c("_mycRasOrd","_wtOrd","_wt_classSet1Ord")) {
                                #fNameOut4=sub(cohort,"",sub(subsetFFlag,"",sub(subsetFlag,"",sub(orderFlag[2],"",sub(orderFlag[1],"",fNameOut2)))))
                                fNameOut4=sub(cohort,"",sub(subsetFFlag,"",sub(subsetFlag,"",sub(orderFlag[2],"",sub(orderFlag[1],sub("Ord","",orderFlag[1]),fNameOut2)))))
                                #fNameOut4=sub("spearman","kendall",sub(cohort,"",sub(subsetFFlag,"",sub(subsetFlag,"",sub(orderFlag[2],"",sub(orderFlag[1],"",fNameOut2))))))
                                clustInfo=read.table(paste(datadir2[1],fNameOut4,"/clusterInfoFeature",paste(fNameOut4,"_",type4Flag,sep=""),".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
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
                        if (orderFlag[1]%in%c("_mycRasOrd","_wtOrd","_wt_classSet1Ord")) {
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
                            if (sepFClustFlag & type2Flag%in%c("_reducedBioFeatPC","_reducedBioFeatMean")) {
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
                    
                    if (predictionStrengthFlag) {
                        if (type3Flag==typeFlag) {
                            timeStamp=Sys.time()
                            print(format(timeStamp, "%x %X"))
                            set.seed(98765)
                            res=prediction.strength(distMat,2,15,M=10,clustermethod=hclustCBI,classification=classifFlag,method="ward.D2",scaling=FALSE)
                            timeStamp=c(timeStamp,Sys.time())
                            print(format(timeStamp[2], "%x %X"))
                            print(diff(timeStamp))
                            save(res,file=paste("predictionStrength_",classifFlag,fNameOut4,".RData",sep=""))
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
                if (!sepFClustFlag | (sepFClustFlag & type3Flag==typeFlag & cohort=="_mycRasWt")) {
                #if (!sepFClustFlag | (sepFClustFlag & type3Flag==typeFlag)) {
                #if (!sepFClustFlag) {
					colCol=matrix(nrow=length(varList),ncol=nrow(phen2))
					for (varId in 1:length(varList)) {
                        if (varList[varId]%in%c("sd")) {
                            #if (varList[varId]%in%c("sd")) {
							j=match(phen2$id,phenAll$id)
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
                            j=match(phen2$id,phenAll$id)
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
							x=x[match(phen2$id,phenAll$id)]
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
				print("summary(c(arrayData))")
				print(summary(c(arrayData)))
				print("quantile(abs(c(arrayData)),probs=seq(0,1,by=.1),na.rm=T)")
				print(quantile(abs(c(arrayData)),probs=seq(0,1,by=.1),na.rm=T))
				main=NULL
				main=header
                
                if (class(clustC)=="hclust") nClust[2]=ifelse(ncol(arrayData)>10,10,NA)
				
                cat("===========",type3Flag,"3 =========\n\n")
				subDir=""
				subDir <- paste(fNameOut,sep="")			
				if (!file.exists(subDir)){
					dir.create(file.path(subDir))
				}
				subDir=paste(subDir,"/",sep="")
				if (outFormat=="png") {
                    if (orderFlag[1]%in%c("_mycRasOrd","_wtOrd","_wt_classSet1Ord") & sepFClustFlag & type3Flag!=typeFlag) {
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
				hcc=heatmap3(x=arrayData, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=colRow, labCol=nameCol, labRow=nameRow, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit2,cexCol=2, , high=colHM[1], low=colHM[2], mid=colHM[3])
				dev.off()

				## -------------------
				if (is.na(nClust[1])) {
                    tbl=cbind(annFeat,clustId="cluster1",order=1:nrow(annFeat))
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
				
				if (is.na(nClust[2])) {
                    tbl=cbind(phen2,clustId="cluster1",order=1:nrow(phen2))
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
					
					tbl=cbind(phen2[clustC$order,],clustId,order=1:nrow(phen2))
				}
                write.table(tbl, paste(subDir,"clusterInfoSample",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
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
						width = 480; height = 140
					} else {
						width = 480; height = 480
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

## -------------------
## Separate out the feature clusters with good profiles from those with bad profiles
## Use 15 clusters

cohort="_mycRas"
cohort="_wt"

for (cohort in c("_mycRas","_wt","_mycRasWt")) {
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
           }
           )

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
        #	clustInfo=read.table(paste(datadir,fileList[gId],sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            clustInfo=read.table(paste(datadir,"clusterInfoFeature",cohort,"_",grpUniq[gId],"_kendall.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            clustInfo$good=0
            k=which(clustInfo$clustId%in%paste("cluster",good[[grpUniq[gId]]],sep=""))
            clustInfo$good[k]=1
        } else {
            clustInfo=ann[which(ann$type==grpUniq[gId]),]
            clustInfo$clustId=clustInfo$featUniq
            clustInfo$good=0
            clustInfo$good[which(clustInfo$clustId%in%clustInfo$clustId[duplicated(clustInfo$clustId)])]=1
        }
        
        clustUniq=unique(clustInfo$clustId[which(clustInfo$good==1)])
        cellMetaMean=cellMeta=matrix(nrow=nrow(arrayData2),ncol=length(clustUniq),dimnames=list(rownames(arrayData2),clustUniq))
        for (cId in 1:length(clustUniq)) {
            cat("\n\n---------------",clustUniq[cId],"\n")
            j=which(colnames(arrayData)%in%clustInfo$feature[which(clustInfo$clustId==clustUniq[cId])])
            print(ann$feature[j])
            cellMetaMean[,clustUniq[cId]]=apply(arrayData2[,j],1,mean,na.rm=T)
    #		fit=prcomp(arrayData[,j], center=T, scale=T)
            fit=prcomp(arrayData2_1[,j], center=F, scale=F)
            cellMeta[,clustUniq[cId]]=fit$x[,1]
            png(paste("plot_meanVpca1_cellMeta","_",clustUniq[cId],cohort,"_",grpUniq[gId],".png",sep=""))
            plot(apply(arrayData2_1[,j],1,mean,na.rm=T),fit$x[,1],main=paste(grpUniq[gId],": ",clustUniq[cId],sep=""),xlab="mean",ylab="pca1")
            dev.off()
            png(paste("screeplot_cellMeta","_",clustUniq[cId],cohort,"_",grpUniq[gId],".png",sep=""))
            plot(fit,main=paste(grpUniq[gId],": ",clustUniq[cId],sep=""))
            dev.off()
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
            #	cell_rf=cbind(cell_rf,cell[,j],cellMeta)
                cell_rf=cbind(cell_rf,arrayData2_1[,j],cellMeta)
                ann_rf=rbind(ann_rf,ann[j,colId],data.frame(feature=colnames(cellMeta),type=typeThis,stringsAsFactors=F))
            #	cell_rfm=cbind(cell_rfm,arrayData2[,j],cellMetaMean)
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
        }
    )
}

save.image("tmp.RData")

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

####################################################################
####################################################################
## PCA

cohort="_mycRas"
cohort="_wt"

featureFlag=""
featureFlag="_reducedFeatMean"
featureFlag="_reducedFeatPCA"
featureFlag="_reducedBioFeatPCA"
featureFlag="_reducedBioFeatMean"

subsetList=c("","_zernike1Only")
subsetList=c("")

plotCutoffFlag=T
plotCutoffFlag=F

cohortName=ifelse(cohort=="_mycRas","Myc/Ras","Wildtype")
switch(cohort,
	   "_mycRas"={
	   cell=cellM
	   },
	   "_wt"={
	   cell=cellW
	   }
	   )

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
#		par(mfcol=c(2,2))
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
