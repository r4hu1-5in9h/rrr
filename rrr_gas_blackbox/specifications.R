get_data<-function(spec,quintile){
  
  df  <- read.dta("gasoline_final_tf1.dta")
  names(df)
  
  N=nrow(df)
  
  # take logs of continuous vars
  cols.to.log<-c("gas","price","income","age","distance")
  df[,cols.to.log]<-lapply(df[,cols.to.log],log)
  
  # output
  Y=df$gas
  
  # construct vars
  df$age2<-((df$age)^2)
  df$income2<-((df$income)^2)
  
  # ensure diff memory location
  df.up<-data.frame(df)
  df.down<-data.frame(df)
  
  # construct df.up and df.down
  prices<-df$price
  delta=sd(prices)/4
  
  df.up$price<-prices+delta/2
  df.down$price<-prices-delta/2

  df$price2<-(prices)^2
  df.up$price2<-(df.up$price)^2
  df.down$price2<-(df.down$price)^2
  
  # specification
  if (spec==1){
    formula<- ~ (factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+price+price2+
      price:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year)))+
      price2:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year)))+urban+youngsingle+distance+distance^2+
      age+age2 + income +income2
  } else {
    formula<- ~ (factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+price+price2+
      price:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+age+age2+income+income2)+
      price2:((factor(driver)+factor(hhsize)+factor(month)+factor(prov)+factor(year))+age+age2+income+income2)+
      urban+youngsingle+distance+distance^2+
      age+age2 + income +income2
  }
  
  regressors<-model.matrix(formula,data=df)
  regressors.up<-model.matrix(formula, data=df.up)
  regressors.down<-model.matrix(formula, data=df.down)
  
  # check
  dim(regressors)
  dim(regressors.up)
  dim(regressors.down)
  
  if (quintile>0){
    q <- ntile(df$income, 5)
    Y.q=Y[q==quintile]
    regressors.q=regressors[q==quintile,]
    regressors.up.q=regressors.up[q==quintile,]
    regressors.down.q=regressors.down[q==quintile,]
  } else {
    Y.q=Y
    regressors.q=regressors
    regressors.up.q=regressors.up
    regressors.down.q=regressors.down
  }
  
  return(list(Y.q,regressors.q,regressors.up.q,regressors.down.q,delta))
  
}


