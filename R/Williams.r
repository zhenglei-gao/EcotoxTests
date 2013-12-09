### Note!
### Throughout this coding, n is a vector of sample size from 
### dose level  0 to k, and k/K is number of treatments not including
### control. 

##'@export
##'
powChow <- function(n=NULL,delta=NULL,sigma=1,tK=NULL,sig.level=0.05,power=NULL,
                    alternative = c("two.sided", "one-tailed")){
  ## powChow(G=5,power=0.8,tK=1.75,sigma=0.22,delta=0.11)
  ## powChow(n=c(4,4,4,4,4),sigma=0.22,delta=1:5)
  S2 <- sigma^2
  
  if(is.null(power)){
    n0 <- n[1]
    #G <- length(n)
    tK<- c(NA,interpWilliamTable(n=n,alpha=sig.level))
    names(tK)[1] <- 0
    power <- 1-pnorm(tK-delta/(sigma*sqrt(2/n)))
  }else{
    
    n0 <- 2*S2*(tK+qnorm(power))^2/(delta^2)
  }
  return(list(n0=n0,power=power))
}

##'@export
powWilliams <- function(data,trend="downward",Ha=2,method=c("approximate","tabulation","simulation","bootstrap"),alpha=0.01,simanz=9999,...){
  ## powWilliams(data=data1,alpha=0.05,method="approximate")
  ## powWilliams(data=data1,alpha=0.05,method="simulation")
  ## ## powWilliams(data=data1,alpha=0.05,method="tabulation")
  ## power calculation based on the alternative being true, which means the level 2 dosage is not the NOEC.
  method <- match.arg(method)
  ## calculate n and delta and S2
  dose <- data[,1]
  d <- by(data[,-1],data[,1],colMeans)
  d <- matrix(unlist(d),ncol=ncol(data)-1,byrow=TRUE)
  if(trend=="downward"){
    d <- data.frame(dose=levels(data[,1]),-d)} else d <- data.frame(dose=levels(data[,1]),d)
  names(d)<- colnames(data)
  n <- by(data[,-1],data[,1],nrow)
  n <- as.vector(unlist(n))
  k <- length(levels(dose))
  res <- matrix(NA,k,ncol(data)-1)
  rownames(res) <- 0:(k-1)
  colnames(res) <- colnames(data)[-1]
  for(coln in 2:ncol(data)){
    mu <- MLE_Williams(d[,coln],n)
    delta <- mu[Ha]-mu[1]
    if(trend=="downward") S2 <- sum((data[,coln]+rep(d[,coln],n))^2)/(sum(n)-k) else S2 <- sum((data[,coln]-rep(d[,coln],n))^2)/(sum(n)-k)
    if(method=="approximate"){
      res[,coln-1] <- powChow(n=n,delta=mu-mu[1],sigma=sqrt(S2))$power
    }else{
      if(method=="simulation"){
        ## simulation based method
        ## Generate data from the alternative
        browser()
        tobs <- (mu-d[1,coln])/sqrt(S2/n[1]+S2/n)      
        talt <- simTbar(n,mu=(mu-mu[1])/sqrt(S2*2/n[1]),simanz=simanz)
        beta <- sapply(1:k,function(k)sum(tobs[k] > talt[k,])/(1+simanz))
        res[,coln-1] <- 1-beta
      }else{
        if(method=="tabulation"){
          ### use both tabulation and simulations!
          ### count the number that we reject the NULL by tabulation!
          talt <- simTbar(n,mu=(mu-mu[1])/sqrt(S2*2/n[1]),simanz=simanz)
          ## calculate the observed MLEs and test statistics.
          taltobs <- apply(simTbar,1,function(x) x)
          ## check for each doselevel, if these <alpha
          
          ## for each dose level calculate the power.
        }
      }
    }
  }
  
  return(res)
}
##'Simualte the test statistic under hypothesis $\mu$
##'@param n sample size for all dose levels
##'@param mu mean at each level, default being 0
##'@param simanz number of simulations.
##'@return simulated t statistic for each dose level, a matrix of k*simanz.
##'@export
simTbar <- function(n,mu=NULL,simanz=9999){
  k <- length(n)
  nu <- sum(n)-k## here the k including the control level.
  maxn <- max(n)
  count <- 0
  if(is.null(mu)) mu <- rep(0,k)
  tbar <- sapply(1:simanz, function(x){
    data <- rnorm(sum(n),rep(mu,n)) ###mu 0, sd 1
    index <-c(0,cumsum(n))
    d <- rep(NA,k)
    for(i in 1:k) d[i] <- mean(data[(index[i]+1):index[i+1]])
    muMLE <- MLE_Williams(d,n)
    s2 <- sum((data-rep(d,n))^2)/nu
    #t <- (mu[k]-d[1])/(sqrt(s2/n[1]+s2/n[k]))
    tvec <- (muMLE-d[1])/(sqrt(s2/n[1]+s2/n))
    if(tvec[1]!=0) browser()
    return(tvec)
  })
  rownames(tbar) <- 0:(length(n)-1)
  return(tbar)
}

##'Function to calculate quantiles of the simulated tbar distribution.
##'
##' @param p the quantile to be calculated 
##'@param n sample size for all dose levels
##'@param mu mean at each level, default being 0
##'@param simanz number of simulations.
##'@return quantile of simulated distribution.
##'@export
qWilliamTbar <- function(p,n,mu=NULL,simanz=9999){
  ###@example \dontrun{
  ### qWilliamTbar(p=c(0.05,0.95),n=c(4,4,4))
  ### simTab <- sapply(1:1000,function(i){qWilliamTbar(p=c(0.05,0.95),n=c(4,4,4,4))})
  ### simTab_no4k4 <- sapply(1:1000,function(i){qWilliamTbar(p=c(0.05,0.95),n=c(4,4,4,4,4))})
  ###}
  tbar <- simTbar(n,mu,simanz)
  apply(tbar,1,function(x)quantile(x,p))
}
##'@export
percent <- function(x, digits = 2, format = "f", ...)
{
  paste(formatC(100 * x, format = format, digits = digits, ...), "%", sep = "")
}
####
##'@export
gCriticalval <- function(p=c(0.95,0.99),n=c(4,4,4,4,4),dof=NULL,ndose=NULL,method=c("tabulation","simulation"),...){
  ##example:gCriticalval(p,n=c(4,4,4),method="simulation",simanz=100000)
  ## example:gCriticalval(n=c(4,4,4))
  method <- match.arg(method)
  if(method=="tabulation"){
    if(all(p %in% c(0.99,0.95))){
        res <- sapply(p,function(x){
          if(x==0.99) {         
            return(interpWilliamTable(dof=dof,n,alpha=0.01,method="inverse"))
          }else{
            return(interpWilliamTable(dof=dof,n,alpha=0.05,method="inverse"))            
          }
        })
        colnames(res)<- percent(p,digits=0)
        return(res)
    }else{
      stop("There are only 5% and 1% tablulated values!")
    }
  }
  if(method=="simulation"){
    return(t(qWilliamTbar(p,n,...)))
  }
}
##'@export
calcWilliams <- function(data,trend=c("downward","upward","auto"),includecontrol= TRUE,method=c("tabulation","simulation","both"),critical=TRUE,...){
  ###@example \dontrun{
  ###res <- calcWiiliams(data=data1,trend="downward",includecontrol=TRUE)
  ### res <- calcWiiliams(data=data1,trend="downward",method="simulation")
  ###}
  trend <- match.arg(trend)
  method <- match.arg(method)
  ## Note that the data structure: 
  ## first column is dose level!
  ## the first level of dose has to be the control level, 
  ## in this case is 0!!
  dosecol <- 1
  dose <- data[,dosecol]
  if(!is.factor(dose)) dose <- factor(dose)
  d <- by(data[,-dosecol],data[,dosecol],colMeans)
  d <- matrix(unlist(d),ncol=ncol(data)-1,byrow=TRUE)
  if(trend=="downward"){
    d <- data.frame(dose=levels(dose),-d)} else d <- data.frame(dose=levels(dose),d)
  names(d)<- colnames(data)
  n <- by(data[,-dosecol],data[,dosecol],nrow)
  n <- as.vector(unlist(n))
  k <- length(levels(dose)) ## including the control level!!!!
  res <- list()
  
  for(coln in 2:ncol(data)){
    res[[coln-1]] <- data.frame(Dose=levels(dose),Mean=ifelse(trend=="downward",-1,1)*d[,coln],MLE=rep(NA,k),Tbar=rep(NA,k))
    mu <- rep(NA,k)
    for(i in 1:k){
      c<- rep(NA,i-1+1)
      for(u in 1:i){
        a <- rep(NA,k-i+1)
        for(v in i:k) a[v-i+1]<- sum(n[u:v]*d[u:v,coln])/sum(n[u:v])
        c[u]<- min(a)
      }
      mu[i] <- max(c)
    }
    ## # MLE is obtained by mu
    res[[coln-1]]$MLE <- mu
    ## For Every level in the dosage
    if(trend=="downward") S2 <- sum((data[,coln]+rep(d[,coln],n))^2)/(sum(n)-k) else S2 <- sum((data[,coln]-rep(d[,coln],n))^2)/(sum(n)-k)
    if(includecontrol==FALSE) {
      if(trend=="downward") S2 <- S2 - sum((data[1:n[1],coln]+d[1,coln])^2)/(sum(n)-k) else S2 <- S2 - sum((data[1:n[1],coln]-d[1,coln])^2)/(sum(n)-k)
    }
    index <- c(0,cumsum(n))
    d1 <- rep(d[,coln],n)
    ###
    tobs <- (mu[K]-d[1,coln])/sqrt(S2/n[1]+S2/n[K])
    res[[coln-1]]$Tbar[K] <- tobs
    if(method=="simulation"||method=="both"){
      tbar <- simTbar(n,mu,simanz)
      WillSim <- t(apply(tbar,1,function(x)quantile(x,p=c(0.95,0.99))))
      res[[coln-1]]$pval <- sapply(1:k,function(i) {1- sum(tbar[i,]<tobs[i])/(1+simanz)})
      res[[coln-1]]$WillSim <- WillSim
      if(critical==TRUE){
        res[[coln-1]]$critical_sim.95 <- WillSim[,1]*sqrt(S2/n[1]+S2/n)+d[1,coln]
        res[[coln-1]]$critical_sim.99 <- WillSim[,2]*sqrt(S2/n[1]+S2/n)+d[1,coln]
      }
    }
    if(method=="tabulation"||method=="both"){
      
      WillTab <- rbind(rep(NA,2),gCriticalval(p=c(0.95,0.99),n=n) )
      rownames(WillTab)[1] <- 0
      res[[coln-1]]$pvaltab <- sapply(1:k,function(i) getPval(tobs[i],q2=WillTab[i,]))
      res[[coln-1]]$WillTab <- WillTab
      if(critical==TRUE){
        res[[coln-1]]$critical_tab.95 <- WillTab[,1]*sqrt(S2/n[1]+S2/n)+d[1,coln]
        res[[coln-1]]$critical_tab.99 <- WillTab[,2]*sqrt(S2/n[1]+S2/n)+d[1,coln]
        
      }
          }
    for(K in k:2){
      tobs <- (mu[K]-d[1,coln])/sqrt(S2/n[1]+S2/n[K])
      res[[coln-1]]$Tbar[K] <- tobs
      ##########
      if(method=="simulation"||method=="both"){
        if(critical==TRUE){
          res[[coln-1]]$pval[K] <- pres$p
          
        }else{
          res[[coln-1]]$pval[K] <- pres
        }
      }else{
        if(method=="tabulation"||method=="both"){
          ptab <- try(gCriticalval(n=n[1:K],method="tabulation"),silent=TRUE)
          
          if(!class(ptab)=="try-error"){
            res[[coln-1]]$pvaltab[K] <- getPval(tobs,ptab)
            res[[coln-1]]$critical_tab.95[K] <- ptab[1]*sqrt(S2/n[1]+S2/n[K])+d[1,coln]
            res[[coln-1]]$critical_tab.99[K] <- ptab[2]*sqrt(S2/n[1]+S2/n[K])+d[1,coln]
          }
        }
      }
      
    }
  }
  names(res)<-colnames(data)[-1]
  return(res)
}
getPval <- function(tobs=2.00,q2=c(1.93,2.90)){
  if(is.na(tobs) || any(is.na(q2))){
    return("NA")
  }else{
  if(tobs>q2[2]){
    return("<1%")
  }else{
    if(tobs>q2[1]){
      return("<5%")
    }else{
      return(">5%")
    }
  }
}}
##'@export
MLE_Williams <- function(d,n){
  k <- length(n)
  mu <- rep(NA,k)
  for(i in 1:k){
    c<- rep(NA,i-1+1)
    for(u in 1:i){
      a <- rep(NA,k-i+1)
      for(v in i:k) a[v-i+1]<- sum(n[u:v]*d[u:v])/sum(n[u:v])
      c[u]<- min(a)
    }
    mu[i] <- max(c)}
  return(mu)
}
##'@export
simNull <- function(n,tobs,simanz=9999,critical=TRUE,qp=c(0.95,0.99),...){
  ### simulate the Null distribution to get the p-value.
  
  k <- length(n)
  if(length(tobs)==1) tobs <- rep(tobs,k)
  nu <- sum(n)-k
  maxn <- max(n)
  count <- 0
  mu <- rep(0,k)
  tnull <- sapply(1:simanz, function(x){
    data <- rnorm(sum(n),rep(mu,n))
    index <-c(0,cumsum(n))
    d <- rep(NA,k)
    for(i in 1:k) d[i] <- mean(data[(index[i]+1):index[i+1]])
    muMLE <- MLE_Williams(d,n)
    s2 <- sum((data-rep(d,n))^2)/nu
    ##t <- (muMLE[k]-d[1])/(sqrt(s2/n[1]+s2/n[k]))
    tvec <- (muMLE-d[1])/(sqrt(s2/n[1]+s2/n))
    return(tvec)
  })
  ###
  p <- sapply(1:k,function(i){1- sum(tnull[i,]<tobs[i])/(1+simanz)})
  if(critical==TRUE) {
    return(list(p=p,q=sapply(1:k,function(i)quantile(tnull[i,],qp))))
  }else{
    return(p)
  }
}

##'@export
##'@example interpWilliamTable(n=c(4,4,4,4,4),alpha=0.05)
interpWilliamTable <- function(dof=NULL,n=c(4,4,4,4),doselevel=NULL,alpha=0.05,method=c("inverse"),dofmethod=c("I","SAS")){
  ## http://stats.stackexchange.com/questions/64538/how-do-i-find-values-not-given-in-interpolate-in-statistical-tables
  if(alpha==0.05){
    data(tbar.tab1)
    tbarTab <- tbar.tab1
  }
  if(alpha==0.01){
    data(tbar.tab2)
    tbarTab <- tbar.tab2
  }
  dofmethod <- match.arg(dofmethod)
  ndose <- length(n)-1
  if(is.null(dof)){
    if(dofmethod=="I") v <- sum(n)-ndose-1 ### How to calculate degree of freedom???
    if(dofmethod=="SAS") {
      if(all(n==n[1])) v <- (n[1]-1)*ndose
    }
  }else v <- dof
  if(is.null(doselevel)) {
    calcdoselevel <- as.character(1:(length(n)-1))
  }else{
    calcdoselevel <- as.character(calcdoselevel)
  }
  df <- as.numeric(rownames(tbarTab))
  ## doselevel <- as.numeric(colnames(tbarTab))
  if(method=="inverse"){
    z <- rep(NA,length(calcdoselevel))
    names(z) <- calcdoselevel
    for(k in seq(calcdoselevel)){
      x <-   1/df
      y <- tbarTab[,calcdoselevel[k]]
      z[k] <- approx(x,y,xout=1/v)$y
    }
    
  }
  return(z)
}