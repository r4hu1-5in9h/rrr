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

setwd("~/Documents/research/rrr_401k_blackbox")

#######################
# clean and format data
#######################

source('specifications.R')

for (quintile in 0:5){

#quintile=1

  print(paste0('quintile: '))
  print(paste0(quintile))
  
df  <- read.dta("sipp1991.dta")
spec=3 #spec in (1-3)
#quintile=0 #quintile in (1-5). 0 means all quintiles
data<-get_data(df,spec,quintile) #trimming like Farrell; different than Chernozhukov et al. 

Y=data[[1]]
T=data[[2]]
X=data[[3]] #no intercept

##################
# helper functions
##################

source('primitives.R')
source('stage1.R')

# dictionary
dict=b2 # b for partially linear model, b2 for interacted model. note that b2 appears in stage1.R for NN
p=length(b(T[1],X[1,]))

#p0=dim(X0) used in low-dim dictionary in the stage 1 tuning procedure
p0=ceiling(p/4) 
if (p>60){
  p0=ceiling(p/40)
  
}


D_LB=0 #each diagonal entry of \hat{D} lower bounded by D_LB
D_add=.2 #each diagonal entry of \hat{D} increased by D_add. 0.1 for 0, 0,.2 otw
max_iter=10 #max number iterations in Dantzig selector iteration over estimation and weights

###########
# algorithm
###########

set.seed(1) # for sample splitting

alpha_estimator=0
gamma_estimator=3
bias=0
#alpha_estimator: 0 dantzig, 1 lasso
#gamma_estimator: 0 dantzig, 1 lasso, 2 rf, 3 nn

source('stage2.R')
results<-rrr(Y,T,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator,bias)
printer(results)
for_tex(results)

}