library(preprocessCore)
library(WGCNA)


#Phenotypes and metabolome in Croatian datasets#

setwd("/common/WORK/SCHOOL2016/Data/Croatia/")
data<-read.table("rssso2015_combinedphenos_VIS_KORCULA_July2015.txt", header=T,stringsAsFactors=T)
glycansVis_LCMS<-read.table("nanoLCMS/Vis_IgG_nanoLCMS_raw_batchCorrected.txt", header=T,stringsAsFactors=T) #Glycans#
glycansVis_UPLC<-read.table("Vis_IgG_UPLC_raw_batchCorrected.txt", header=T,stringsAsFactors=T) #Glycans#

glycansVis_LCMS<-glycansVis_LCMS[,-2]
glycansVis_UPLC<-glycansVis_UPLC[,-c(2,3)]


#remove summary percentages of families of lipids#
index<-grep("p_",names(data))
names(data)[index]
data.cleaned<-data[,-index]

#remove rations between families#
index<-grep("BY",names(data.cleaned))
names(data.cleaned)[index]
data.cleaned<-data.cleaned[,-index]

#Add phenotypes and lipids to the database#
#DataGlycans<-merge(glycansVis_LCMS,data.cleaned[,c("id","sex","age","isVis","bmi","unrelated_gkin0625")],by.x="id",by.y="id") 
DataGlycans<-merge(glycansVis_UPLC,glycansVis_LCMS,by.x="id",by.y="id") #less common individuals, consider LCMS and UPLC separatly?
DataALL<-merge(DataGlycans,data.cleaned,by.x="id",by.y="id") #less common individuals, consider LCMS and UPLC separatly?


#Choose unrelated individuals#
#DataALLUnrelated<-DataALL[DataALL$unrelated_gkin0625==T,]


##############################################################################################
##############################################################################################

#First create matrix containing only Glycans/metabolites#
#Example with matrix of Glycans#

#OUTLIERS DETECTION#
sampleTree = hclust(dist(GLYCANS), method = "average");

sizeGrWindow(12,9)
pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)


#Normalization: use package preprocessing#
#Example: Quantile normalization#
GLYCANS.QN <- data.frame(normalize.quantiles(as.matrix(GLYCANS[,-1])))
colnames(GLYCANS.QN)<-names(GLYCANS)[-1]



