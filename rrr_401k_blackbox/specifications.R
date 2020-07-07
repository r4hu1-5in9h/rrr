get_data<-function(df,spec,quintile){
  
  Y <- df[,"net_tfa"]
  T <- df[,"e401"]
  
  ## low dim specification
  X.L <- cbind(df[,"age"], df[,"inc"], df[,"educ"], df[,"fsize"], df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"])
  
  ## high dim specification. NOTE: original paper is this spec squared (pairwise interactions thereof?)
  X.H <- cbind(poly(df[,"age"], 6, raw=TRUE),
               poly(df[,"inc"], 8, raw=TRUE),
               poly(df[,"educ"], 4, raw=TRUE),
               poly(df[,"fsize"], 2, raw=TRUE),
               df[,"marr"], df[,"twoearn"], df[,"db"], df[,"pira"], df[,"hown"]) 
  
  # copied from EJ: "(poly(age, 6, raw=TRUE) + poly(inc, 8, raw=TRUE) + poly(educ, 4, raw=TRUE) + poly(fsize, 2, raw=TRUE) + marr + twoearn + db + pira + hown)^2"
  X.vH=model.matrix(~(poly(df[,"age"], 6, raw=TRUE) + 
                        poly(df[,"inc"], 8, raw=TRUE) + 
                        poly(df[,"educ"], 4, raw=TRUE) + 
                        poly(df[,"fsize"], 2, raw=TRUE) + 
                        df[,"marr"] + 
                        df[,"twoearn"] + 
                        df[,"db"] + 
                        df[,"pira"] + 
                        df[,"hown"])^2)
  X.vH=X.vH[,-1]
  
  if (spec==1){
    X=X.L
  } else if (spec==2) {
    X=X.H
  } else {
    X=X.vH
  }
  
  X <- scale(X,center=TRUE,scale=TRUE)
  
  #impose common support
  p.1 <- multinom(T~X-1, trace=FALSE)$fitted.values
  indexes.to.drop <- which(p.1 < min(p.1[T==1]) | max(p.1[T==1]) < p.1)
  if (length(indexes.to.drop)==0) {indexes.to.drop <- n+1}	#R throws a wobbly if [-indexes.to.drop] is negating an empty set. 
  n.per.treatment <- as.vector(table(T[-indexes.to.drop]))
  n.trim <- n.per.treatment[1]+n.per.treatment[2]
  
  Y.trimmed=Y[-indexes.to.drop]
  T.trimmed=T[-indexes.to.drop]
  X.trimmed=X[-indexes.to.drop,]
  
  if (spec==1){
    inc=X.trimmed[,2]
  } else if (spec==2) {
    inc=X.trimmed[,7]
  } else {
    inc=X.trimmed[,7]
  }
  
  if (quintile>0){
    q <- ntile(inc, 5)
    Y.q=Y.trimmed[q==quintile]
    T.q=T.trimmed[q==quintile]
    X.q=X.trimmed[q==quintile,]
  } else {
    Y.q=Y.trimmed
    T.q=T.trimmed
    X.q=X.trimmed
  }
  
  return(list(Y.q,T.q,X.q))
  
}