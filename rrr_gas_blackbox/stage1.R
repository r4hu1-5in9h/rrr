#stage 1 function
l=0.1

RMD_dantzig <- function(M, G, D, lambda=0, sparse = TRUE) {
  
  p <- ncol(G)
  zp <- rep(0, p)
  L <-c(l,rep(1,p-1)) #dictionary is ordered (constant,T,X,TX)
  
  A <- solve(diag(D),G)
  R <- rbind(A, -A)
  
  a <- solve(diag(D),M)
  r <- c(a - lambda*L, -a - lambda*L)
  
  if(sparse) {
    Ip <- as(p, "matrix.diag.csr")
    R <- as.matrix.csr(R)
    f <- rq.fit.sfnc(Ip, zp, R = R, r = r)
  } else {
    Ip <- diag(p)
    f <- rq.fit.fnc(Ip, zp, R = R, r = r)
  }
  
  return(f)
}

RMD_lasso <- function(M, G, D, lambda=0, control = list(maxIter = 1000, optTol = 10^(-5), 
                                                        zeroThreshold = 10^(-6)), beta.start = NULL) {
  
  p <- ncol(G)
  
  # v1 - equivalent to v3
  #Gt<- solve(diag(D),G)
  #Mt <- solve(diag(D),M)
  
  # v2 - incorrect
  #D_inv<-diag(1/D)
  #Gt <- D_inv %*% G %*% D_inv
  #Mt <- solve(diag(D),M)
  
  # v3 - here and below. final
  Gt<-G
  Mt<-M
  
  L <-c(l,rep(1,p-1)) #dictionary is ordered (constant,...)
  lambda_vec=lambda*L*D #v3: insert D here
  
  if (is.null(beta.start)) {
    beta <- rep(0,p) #vs low-dimensional initialization
  }
  else {
    beta <- beta.start
  }
  wp <- beta
  mm <- 1
  while (mm < control$maxIter) {
    beta_old <- beta
    for (j in 1:p) {
      rho=Mt[j]-Gt[j,]%*%beta+Gt[j,j]*beta[j]
      z=Gt[j,j]
      
      if (sum(is.na(rho)) >= 1) {
        beta[j] <- 0
        next
      }
      if (rho < -1 * lambda_vec[j]) 
        beta[j] <- (rho+lambda_vec[j])/z
      if (abs(rho) <= lambda_vec[j]) 
        beta[j] <- 0
      if (rho > lambda_vec[j]) 
        beta[j] <- (rho-lambda_vec[j])/z
    }
    wp <- cbind(wp, beta)
    if (sum(abs(beta - beta_old), na.rm = TRUE) < control$optTol) {
      break
    }
    mm <- mm + 1
  }
  w <- beta
  w[abs(w) < control$zeroThreshold] <- 0
  return(list(coefficients = w, coef.list = wp, num.it = mm))
}


get_D <- function(Y,X,X.up,X.down,delta,m,rho_hat){
  n=dim(X)[1]
  p=dim(X)[2]
  
  df=matrix(0,p,n)
  for (i in 1:n){
    df[,i]=X[i,]*as.vector(rho_hat %*% X[i,])-m(Y[i],X[i,],X.up[i,],X.down[i,],delta,b)
  }
  df=df^2
  D2=rowMeans(df)
  
  D=sqrt(D2)
  return(D) #pass around D as vector
}

c=0.5
alpha=0.1
tol=1e-6

RMD_stable<-function(Y,X,X.up,X.down,delta,p0,D_LB,D_add,max_iter,is_alpha,is_lasso){
  
  n=dim(X)[1]
  p=dim(X)[2]
  
  # low-dimensional moments
  X0=X[,1:p0]
  X0.up=X.up[,1:p0]
  X0.down=X.down[,1:p0]
  MNG0<-get_MNG(Y,X0,X0.up,X0.down,delta)
  M_hat0=MNG0[[1]]
  N_hat0=MNG0[[2]]
  G_hat0=MNG0[[3]]
  
  # initial estimate
  rho_hat0=solve(G_hat0,M_hat0)
  rho_hat=c(rho_hat0,rep(0,p-ncol(G_hat0)))
  beta_hat0=solve(G_hat0,N_hat0)
  beta_hat=c(beta_hat0,rep(0,p-ncol(G_hat0)))
  
  # moments
  MNG<-get_MNG(Y,X,X.up,X.down,delta)
  M_hat=MNG[[1]]
  N_hat=MNG[[2]]
  G_hat=MNG[[3]]
  B=MNG[[4]]
  
  # penalty
  lambda=c*qnorm(1-alpha/(2*p))/sqrt(n) # snippet
  
  if(is_alpha){ 
    ###########
    # alpha_hat
    ###########
    diff_rho=1
    k=1
    while(diff_rho>tol & k<=max_iter){
      
      # previous values
      rho_hat_old=rho_hat+0
      
      # normalization
      D_hat_rho=get_D(Y,X,X.up,X.down,delta,m,rho_hat_old)
      D_hat_rho=pmax(D_LB,D_hat_rho)
      D_hat_rho=D_hat_rho+D_add
      
      # RMD estimate
      if(is_lasso){
        rho_hat=RMD_lasso(M_hat, G_hat, D_hat_rho, lambda)$coefficients
      }else{
        rho_hat=RMD_dantzig(M_hat, G_hat, D_hat_rho, lambda)$coefficients
      }
      
      # difference
      diff_rho=two.norm(rho_hat-rho_hat_old)
      k=k+1
      
    }
    
    print(paste0('k: '))
    print(paste0(k))
    return(rho_hat)
    
  } else { 
    ###########
    # gamma_hat
    ###########
    diff_beta=1
    k=1
    while(diff_beta>tol & k<=max_iter){
      
      # previous values
      beta_hat_old=beta_hat+0
      
      # normalization
      D_hat_beta=get_D(Y,X,X.up,X.down,delta,m2,beta_hat_old)
      D_hat_beta=pmax(D_LB,D_hat_beta)
      D_hat_beta=D_hat_beta+D_add
      
      # RMD estimate
      if(is_lasso){
        beta_hat=RMD_lasso(N_hat, G_hat, D_hat_beta, lambda)$coefficients
      }else{
        beta_hat=RMD_dantzig(N_hat, G_hat, D_hat_beta, lambda)$coefficients
      }
      
      # difference
      diff_beta=two.norm(beta_hat-beta_hat_old)
      k=k+1
      
    }
    
    print(paste0('k: '))
    print(paste0(k))
    return(beta_hat)
    
  }
}

arg_Forest<- list(clas_nodesize=1, reg_nodesize=5, ntree=1000, na.action=na.omit, replace=TRUE)
arg_Nnet<- list(size=8,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE)

get_stage1<-function(Y,X,X.up,X.down,delta,p0,D_LB,D_add,max_iter,alpha_estimator,gamma_estimator){
  
  n=dim(X)[1]
  p=dim(X)[2]
  MNG<-get_MNG(Y,X,X.up,X.down,delta)
  B=MNG[[4]]
  
  ###########
  # alpha hat
  ###########
  if(alpha_estimator==0){ # dantzig
    
    rho_hat=RMD_stable(Y,X,X.up,X.down,delta,p0,D_LB,D_add,max_iter,1,0)
    alpha_hat<-function(x){
      return(b(x)%*%rho_hat)
    }
    
  } else if(alpha_estimator==1){ # lasso
    
    rho_hat=RMD_stable(Y,X,X.up,X.down,delta,p0,D_LB,D_add,max_iter,1,1)
    alpha_hat<-function(x){
      return(b(x)%*%rho_hat)
    }
    
  }
  
  ###########
  # gamma hat
  ###########
  if(gamma_estimator==0){ # dantzig
    
    beta_hat=RMD_stable(Y,X,X.up,X.down,delta,p0,D_LB,D_add,max_iter,0,0)
    gamma_hat<-function(x){
      return(b(x)%*%beta_hat)
    }
    
  } else if(gamma_estimator==1){ # lasso
    
    beta_hat=RMD_stable(Y,X,X.up,X.down,delta,p0,D_LB,D_add,max_iter,0,1)
    gamma_hat<-function(x){ 
      return(b(x)%*%beta_hat)
    }
    
  } else if(gamma_estimator==2){ # random forest
    
    forest<- do.call(randomForest, append(list(x=B,y=Y), arg_Forest))
    gamma_hat<-function(x){
      return(predict(forest,newdata=b(x), type="response"))
    }
    
  } else if(gamma_estimator==3){ # neural net
    
    # scale down, de-mean, run NN, scale up, remean so that NN works well
    maxs_B <- apply(B, 2, max)
    mins_B <- apply(B, 2, min)
    
    maxs_Y<-max(Y)
    mins_Y<-min(Y)
    
    # hack to ensure that constant covariates do not become NA in the scaling
    const=maxs_B==mins_B
    keep=(1-const)*1:length(const)
    
    NN_B<-B
    NN_B[,keep]<-scale(NN_B[,keep], center = mins_B[keep], scale = maxs_B[keep] - mins_B[keep])
    
    NN_Y<-scale(Y, center = mins_Y, scale = maxs_Y - mins_Y)
    
    nn<- do.call(nnet, append(list(x=NN_B,y=NN_Y), arg_Nnet)) #why is it degenerate with fitted.values=1?
    gamma_hat<-function(x){
      
      test<-t(as.vector(x))
      NN_b<-test
      NN_b[,keep]<-scale(t(NN_b[,keep]), 
                         center = mins_B[keep], 
                         scale = maxs_B[keep] - mins_B[keep])
      
      NN_Y_hat<-predict(nn,newdata=NN_b)
      Y_hat=NN_Y_hat*(maxs_Y-mins_Y)+mins_Y
      
      return(Y_hat)
    }
    
  }
  
  return(list(alpha_hat,gamma_hat))
  
}
