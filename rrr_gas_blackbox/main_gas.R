########
# set up
########

rm(list=ls())

library("foreign")
library("dplyr")
library("ggplot2")
library("quantreg")

library("MASS")
library("glmnet")	#lasso, group lasso, and ridge, for outcome models, LPM, mlogit. Also, unpenalized mlogit via glmnet is far faster than the mlogit package.
library("grplasso")	#glmnet group lasso requires the same n per group (ie multi-task learning), which is perfect for mlogit but wrong for the outcome model.
# library("mlogit")	#slower than unpenalized estimation in glmnet, but glmnet won't fit only an intercept
library("nnet")	#quicker multinomial logit
library("randomForest")
library("gglasso")
library("plotrix")
library("gridExtra")

setwd("~/Documents/research/rrr_gas_blackbox")

#######################
# clean and format data
#######################

for (quintile in 0:5){

#quintile=1

  print(paste0('quintile: '))
  print(paste0(quintile))

spec=1
#1 means Chernozhukov and Semenova
#2 means Chernozhukov and Semenova and additional interactions

data<-get_data(spec,quintile) #like Chernozhukov and Semenova

Y=data[[1]]
X=data[[2]]
X.up=data[[3]]
X.down=data[[4]]
delta=data[[5]]

##################
# helper functions
##################

source('primitives.R')

test<-get_MNG(Y,X,X.up,X.down,delta)
M_hat=test[[1]]
N_hat=test[[2]]
as.numeric(M_hat)
as.numeric(N_hat)

source('stage1.R')

# dictionary is identity
n=nrow(X)
p=ncol(X)

#p0=dim(X0) used in low-dim dictionary in the stage 1 tuning procedure
p0=ceiling(p/4) #p/2 works for low p
if (p>60){
  p0=ceiling(p/40)
}

D_LB=0 #each diagonal entry of \hat{D} lower bounded by D_LB
D_add=.2 #each diagonal entry of \hat{D} increased by D_add. 0.1 for 0, 0,.2 otw
max_iter=10 #max number iterations in Dantzig selector iteration over estimation and weights

###########
# algorithm
###########

alpha_estimator=0
gamma_estimator=3
bias=0
#alpha_estimator: 0 dantzig, 1 lasso
#gamma_estimator: 0 dantzig, 1 lasso, 2 rf, 3 nn

set.seed(1) # for sample splitting

source('stage2.R')
results<-rrr(Y,X,X.up,X.down,delta,p0,D_LB,D_add,max_iter,alpha_estimator,gamma_estimator,bias)
printer(results)
for_tex(results)

}