#require libraries
library(MASS)
library(foreach)
library(doParallel)
library(penalizedSVM)

##################################################################
#Function to simulate the dataset with exponential hazard function
#n is the survival sample size
#p is the number of variables(including noise variable)
#p.noise is the number of noise variables
#cen is if the data is censored data or not
#
#it reurn a list (coef,data) where coef is the coefficients
# of variables
##################################################################
simulation <- function(n=200, p=12, p.noise=5 ,cen=FALSE) {
    #build covariate matrix
    Sigma <- matrix(NA,nrow=p,ncol=p);
    rho <- 0.5;
    for(i in 1:p)
    {
        for(j in 1:p)
        {
            Sigma[i,j] <- rho^abs(i-j);
            Sigma[j,i] <- Sigma[i,j]
        }
    }
    Z <- mvrnorm(n,rep(0,p),Sigma)
    #having the result
    b <-c(runif(p-p.noise,-10,10),rep(0,p.noise))
    e <- log(rexp(n))
    T <- exp(-Z%*%b+e)
    if(cen==TRUE) C <- rexp(n)
    else C <- T
    U <- pmin(T,C)
    d <- as.numeric(T<=C)
    return(list("coef" = b, "data"=data.frame(U,d,Z)))
}

##################################################################
#Function to transform the survival data
#x is the survival data set
#x[,1] is U, x[,2] is delta, x[,3-] is covariates
#n is the survival sample size
#p is the number of variables(including noise variable)
#ncore is the number of computer core you want the program to use
#
#returns a 5+2p dimensonal data
##################################################################
transformation <- function (x,n,p,ncore=2) {
    registerDoParallel(cores=2)
    #ptm = proc.time()
    result <- foreach(i= 1:n,.combine='rbind') %:% 
        foreach(j= 1:n,.combine='rbind') %dopar% {
            if(i!=j) c(x[i,],x[j,])
        }
    #proc.time()-ptm

    lab = as.numeric(unlist(result[,1])<unlist(result[,1+p+2]))
    lab[which(lab==0)]=-1
    W=rep(0,dim(result)[1])
    indices = union(intersect(which(result[,2]==1),which(lab==1))
              ,intersect(which(result[,2+p+2]==1),which(lab==-1)))
    W[indices]=1
    result=cbind(lab,result)
    if(length(which(W==0))!=0) result=result[-which(W==0),]
    return(result)
}

#basic configuration for simulation study
##################################################################
P <- 12
P.noise <- 5
N <- 200

#Simulate the data and transform
##################################################################
tmp <- simulation(N,P,P.noise,cen=TRUE)
tmp2 <- transformation(tmp$data,N,P,ncore=2)

#further step to transform data
##################################################################
rownames(tmp2)=NULL
tmp2=data.frame(tmp2)
lab = data.matrix(unlist(tmp2[,1]))
Zi = data.matrix(tmp2[,4:(4+P-1)])
Zj = data.matrix(tmp2[,(4+P+2):(dim(tmp2)[2])])
Z = Zi-Zj
new.data.set = cbind(lab,Z)

#create train set and test set
##################################################################
test.indices = rbinom(dim(new.data.set)[1],1,0.5)
train.set = new.data.set[which(test.indices==0),]
test.set = new.data.set[which(test.indices==1),]
registerDoParallel(cores=1)

#feed train data to SCAD SVM
##################################################################
ptm = proc.time()
scad.fix<- svm.fs(train.set[,2:dim(train.set)[2]], y=train.set[,1], fs.method="scad",
cross.outer= 0, grid.search = "discrete",
lambda1.set=c(0.01,0.1,0.3,0.5),
parms.coding = "none", show="none",
maxIter = 10, inner.val.method = "cv", cross.inner= 5,
seed=1, verbose=FALSE)
proc.time()-ptm

#predict test data with previous fed SVM
##################################################################
test.error.scad<-predict(scad.fix, newdata=test.set[,2:dim(train.set)[2]]
,newdata.labels=test.set[,1] )
print(test.error.scad$tab)
