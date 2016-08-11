library(glmnet);
library(preprocessCore);
library(pROC);
library(SGL);
library(grplasso);
  

setwd("/common/WORK/SCHOOL/teachers/mar01")
source("/common/WORK/SCHOOL/teachers/mar01/codes/glmnet.2CV.R")
source("/common/WORK/SCHOOL/teachers/mar01/codes/MULTI.GAM.2CV.R")


index.metab<-101:279
index.glycans.UPLC<-2:24 
index.glycans.LCMS<-25:74
index.outcome<-78 #bmi:78
index.age.gender<-c(75,76)
data<-DataALL

index.fit<-c(index.outcome,index.glycans.LCMS,index.age.gender,index.metab)
X.fit<-data[complete.cases(data[,index.fit]),index.glycans.LCMS]
outcome.fit<-data[complete.cases(data[,index.fit]),index.outcome]
outcome.fit<-I(log(outcome.fit))
age.gender.fit<-data[complete.cases(data[,index.fit]),index.age.gender]
data.fit<-data.frame(outcome.fit,X.fit,age.gender.fit)

#Create outer CV partition#
nfolds=5
folds<-createFolds(1:nrow(X.fit), k = nfolds, list = T)


#Model based on age and gender#
formula<-outcome.fit~sex+age
BMI.AgeGender<-MULTI.GAM.2CV(formula,data=data.fit,nfolds=nfolds,folds=folds,family="gaussian")
p.BMI.AgeGender<-outcome.fit-BMI.AgeGender$res

#Model based on LCMS glycans#
BMI.LCMS.Rg<-glmnet.2CV(X=as.matrix(X.fit),Y=outcome.fit,alpha=0,folds=folds,nfolds=nfolds,family="gaussian");
BMI.LCMS.Lasso<-glmnet.2CV(X=as.matrix(X.fit),Y=outcome.fit,alpha=1,folds=folds,nfolds=nfolds,family="gaussian");


BMI.AgeGender.LCMS.Rg<-glmnet.2CV(X=as.matrix(X.fit),Y=outcome.fit,alpha=0,folds=folds,nfolds=nfolds,family="gaussian",offset=p.BMI.AgeGender);
BMI.AgeGender.LCMS.Lasso<-glmnet.2CV(X=as.matrix(X.fit),Y=outcome.fit,alpha=1,folds=folds,nfolds=nfolds,family="gaussian",offset=p.BMI.AgeGender);


#Introduce modification: joint model of sex and age and glycans introducing effect of age and sex as fixed effects#
#Hint: use penalty.factor in function glmnet#

X.fit<-data[complete.cases(data[,index.fit]),index.metab]
BMI.Lipids.Rg<-glmnet.2CV(X=as.matrix(X.fit),Y=outcome.fit,alpha=0,folds=folds,nfolds=nfolds,family="gaussian");
BMI.Lipids.Lasso<-glmnet.2CV(X=as.matrix(X.fit),Y=outcome.fit,alpha=1,folds=folds,nfolds=nfolds,family="gaussian");


BMI.AgeGender.Lipids.Rg<-glmnet.2CV(X=as.matrix(X.fit),Y=outcome.fit,alpha=0,folds=folds,nfolds=nfolds,family="gaussian",offset=p.BMI.AgeGender);
BMI.AgeGender.Lipids.Lasso<-glmnet.2CV(X=as.matrix(X.fit),Y=outcome.fit,alpha=1,folds=folds,nfolds=nfolds,family="gaussian",offset=p.BMI.AgeGender);



#Stacked datasets#

X.fit<-data[complete.cases(data[,index.fit]),c(index.glycans.LCMS,index.metab)]
BMI.Stacked.Rg<-glmnet.2CV(X=as.matrix(X.fit),Y=outcome.fit,alpha=0,folds=folds,nfolds=nfolds,family="gaussian");
BMI.Stacked.Lasso<-glmnet.2CV(X=as.matrix(X.fit),Y=outcome.fit,alpha=1,folds=folds,nfolds=nfolds,family="gaussian");





#Binary outcomes:#
roc(outcome,p.BMI.METAB.Rg,ci=T,ci.method="bootstrap")



#Biological clustering considered as gold-standard grouping (4 classes of LCMS glycans)
clustering<-c(3,3,3,4,4,4,3,3,4,4,2,1,1,2,2,2,1,1,2,2)



#Biological clustering considered as gold-standard grouping (6 classes of metabolites)
clustering<-c(rep(1,49),rep(2,15),rep(3,19),rep(4,15),rep(5,23),rep(6,9))


cluster<-list()
cluster<-lapply(1:nfolds,function(i)cluster[[i]]<-clustering[-1])

BMI.METAB.gglasso<-gglasso.2CV(X=X,Y=outcome,loss="ls",folds=folds,nfolds=nfolds,family="gaussian",cluster=cluster);
Q2.BMI.METAB.gglasso<-BMI.METAB.gglasso$Q2;
p.BMI.METAB.gglasso<-BMI.METAB.gglasso$p;




#Variable selection evaluation (LASSO)#
coef<-matrix(NA,ncol(CleanedMetabolites),nfolds)
for (i in 1:nfolds){
coef[,i]<-as.numeric(coef(BMI.METAB.LS$cv.fit.glmnet[[i]],s=BMI.METAB.LS$cv.fit.glmnet[[i]]$lambda.min))
}

coef<-coef[-1,]
selected.variables<-ifelse(coef>0,1,0)
selected.variables.frec<-apply(selected.variables,1,sum)
tab<-cbind(names(CleanedMetabolites)[-1],as.numeric(selected.variables.frec))

ttt = as.data.frame(tab[tab[,2]>=1,])
ttt[,2] = as.numeric(ttt[,2])
ttt = ttt[order(ttt[,2], decreasing = TRUE),]













