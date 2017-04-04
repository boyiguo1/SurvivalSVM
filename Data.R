#Download and Preprocessing Ovarian Data
##############################################################################
library(curatedOvarianData)
source(system.file("extdata","patientselection.config",package="curatedOvarianData"))
sapply(ls(), function(x) if(!x %in% c("remove.samples", "duplicates")) print(get(x)))
source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))

# Organize survival time and covariates
##############################################################################
surv.info = matrix(as.numeric(esets[[14]]$y),228/2,2)
data = as(phenoData(esets[[14]]),"data.frame")
probeset = exprs(esets[[14]])
dim(probeset)
t(probeset[1:5,1:5])
index = NULL
for(i in 1:length(unique(surv.info[,1]))){
    index=c(index,which(surv.info[,1]==unique(surv.info[,1])[i])[1])
}
data = data.frame(surv.info,t(probeset)[,1:100])
data=data[index,]

# Transform data
##############################################################################
tmp2 = transformation(data,n=dim(data)[1],p=100)
P=100
backup = tmp2
rownames(tmp2)=NULL
tmp2=data.frame(tmp2)
lab = data.matrix(unlist(tmp2[,1]))
Zi = data.matrix(tmp2[,4:(4+P-1)])
Zj = data.matrix(tmp2[,(4+P+2):(dim(tmp2)[2])])
Z = Zi-Zj
new.data.set = cbind(lab,Z)
test.indices = rbinom(dim(new.data.set)[1],1,0.5)
mean(test.indices)
train.set = new.data.set[which(test.indices==0),]
test.set = new.data.set[which(test.indices==1),]
registerDoParallel(cores=1)


# Train SVM
##############################################################################
ptm = proc.time()
scad.fix<- svm.fs(train.set[,2:dim(train.set)[2]], y=train.set[,1], fs.method="scad",
cross.outer= 0, grid.search = "discrete",
lambda1.set=c(0.01,0.1,0.3,0.5),
parms.coding = "none", show="none",
maxIter = 10, inner.val.method = "cv", cross.inner= 5,
seed=1, verbose=FALSE)
proc.time()-ptm

# Predict test data
##############################################################################
test.error.scad<-predict(scad.fix, newdata=test.set[,2:dim(train.set)[2]]
,newdata.labels=test.set[,1] )
print(test.error.scad$tab)