##'@export
powWilliams <- function(data,trend="downward",Ha=2,method=c("approximate","simulation","bootstrap"),alpha=0.01,...){
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
  res <- rep(NA,ncol(data)-1)
  for(coln in 2:ncol(data)){
    mu <- MLE_Williams(d[,coln],n)
    delta <- mu[Ha]-mu[1]
    if(trend=="downward") S2 <- sum((data[,coln]+rep(d[,coln],n))^2)/(sum(n)-k) else S2 <- sum((data[,coln]-rep(d[,coln],n))^2)/(sum(n)-k)
    if(method=="approximate"){
      tnull <- simTbar(n,simanz=9999)
      ##plot(density(tnull))
      tquantile <- quantile(tnull,alpha)
      res[coln-1] <- (1-pnorm(tquantile-delta/sqrt(S2*2/n[1])))
    }else{
      if(method=="simulation"){
        ## simulation based method
        K <- Ha
        tobs <- (mu[K]-d[1,coln])/sqrt(S2/n[1]+S2/n[K])
        simanz=9999
        talt <- simTbar(n[1:K],mu=c(rep(0,Ha-1),(mu[Ha]-mu[1])/sqrt(S2*2/n[1])),simanz)
        beta <- sum(tobs > talt)/(1+simanz)
        res[coln-1] <- 1-beta
      }else{
        
      }
    }
  }
  names(res) <- colnames(data)[-1]
  return(res)
}
##'@export
simTbar <- function(n,mu=NULL,simanz=9999){
  k <- length(n)
  nu <- sum(n)-k## here the k including the control level.
  maxn <- max(n)
  count <- 0
  if(is.null(mu)) mu <- rep(0,k)
  tbar <- sapply(1:simanz, function(x){
    data <- rnorm(sum(n),rep(mu,n))
    index <-c(0,cumsum(n))
    d <- rep(NA,k)
    for(i in 1:k) d[i] <- mean(data[(index[i]+1):index[i+1]])
    mu <- MLE_Williams(d,n)
    s2 <- sum((data-rep(d,n))^2)/nu
    t <- (mu[k]-d[1])/(sqrt(s2/n[1]+s2/n[k]))
    return(t)
  })
  return(tbar)
}

##'@export
qWilliamTbar <- function(p,n,mu=NULL,simanz=9999){
  ###@example \dontrun{
  ###qWilliamTbar(p=c(0.05,0.95),n=c(4,4,4))
  ###}
  tbar <- simTbar(n,mu,simanz)
  quantile(tbar,p)
}
####
##'@export
gCriticalval <- function(p=c(0.95,0.99),n=c(4,4,4),dof=NULL,ndose=NULL,method=c("tabulation","simulation"),...){
  ##example:gCriticalval(p,n=c(4,4,4),method="simulation",simanz=100000)
  ## example:gCriticalval(n=c(4,4,4))
  method <- match.arg(method)
  if(method=="tabulation"){
    if(all(p %in% c(0.99,0.95))){
      vname <- c(5:20,22,24,26,28,30,35,40,60,120,Inf)
      kname <- 1:10
      
      if(!is.null(dof)&& !is.null(ndose)){
        v <- dof
        k <- ndose
      } else{
        k <- length(n)-1
        v <- sum(n)-k-1
      }
      
      if(all(v %in% vname) && all(k %in% kname)){
        rid <- vname==v
        cid <- kname==k
        
        res <- sapply(p,function(x){
          if(x==0.99) {
            data(tbar.tab2)
            return(tbar.tab2[rid,cid])
          }else{
            data(tbar.tab1)
            return(tbar.tab1[rid,cid])
          }
        })
        names(res)<- p
        return(res)
      }else{
        stop("The degree of freedom can only be in c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 35, 40, 60, 120, Inf) and the number of dose levels can only be in 1 to 10.")
      }
    }else{
      stop("There are only 5% and 1% tablulated values!")
    }
  }
  if(method=="simulation"){
    return(qWilliamTbar(p,n,...))
  }
}
##'@export
calcWiiliams <- function(data,trend="downward",includecontrol= TRUE,method=c("tabulation","simulation","both"),critical=TRUE,...){
  ###@example \dontrun{
  ###res <- calcWiiliams(data=data1,trend="downward",includecontrol=TRUE)
  ###}
  method <- match.arg(method)
  ## Note that the data structure: first column is dose level! the first level of dose has to be the control level, in      ## this case is 0!!
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
  k <- length(levels(dose))
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
    if(method=="simulation"||method=="both"){
      res[[coln-1]]$pval <- rep(NA,k)
      if(critical==TRUE){
        res[[coln-1]]$critical_sim.95 <- rep(NA,k)
        res[[coln-1]]$critical_sim.99 <- rep(NA,k)
      }
    }
    if(method=="tabulation"||method=="both"){
      res[[coln-1]]$pvaltab <- rep(NA,k)
      if(critical==TRUE){
        res[[coln-1]]$critical_tab.95 <- rep(NA,k)
        res[[coln-1]]$critical_tab.99 <- rep(NA,k)
      }
    }
    for(K in k:2){
      ## Estimate s2
      #if(trend=="downward") s2 <- sum((data[1:index[K],coln]+d1[1:index[K]])^2)/(sum(n[1:K])-K) else {
      #  s2 <- sum((data[1:index[K],coln]-d1[1:index[K]])^2)/(sum(n[1:K])-K)
      #}
      ##browser()
      ## tobs <- (mu[K]-d[1,coln])/sqrt(s2/n[1]+s2/n[K])
      tobs <- (mu[K]-d[1,coln])/sqrt(S2/n[1]+S2/n[K])
      res[[coln-1]]$Tbar[K] <- tobs
      if(method=="simulation"||method=="both"){
        
        pres<- simNull(n[1:K],tobs,simanz=9999,critical==critical)
        if(critical==TRUE){
          res[[coln-1]]$pval[K] <- pres$p
          res[[coln-1]]$critical_sim.95[K]<-pres$q[1]*sqrt(S2/n[1]+S2/n[K])+d[1,coln]
          res[[coln-1]]$critical_sim.99[K]<-pres$q[2]*sqrt(S2/n[1]+S2/n[K])+d[1,coln]
          
        }else{
          res[[coln-1]]$pval[K] <- pres
        }
      }else{
        if(method=="tabulation"||method=="both"){
          ptab <- try(gCriticalval(n=n[1:K],method=method),silent=TRUE)
          
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
  if(tobs>q2[2]){
    return("<1%")
  }else{
    if(tobs>q2[1]){
      return("<5%")
    }else{
      return(">5%")
    }
  }
}
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
  nu <- sum(n)-k
  maxn <- max(n)
  count <- 0
  mu <- rep(0,k)
  tnull <- sapply(1:simanz, function(x){
    data <- rnorm(sum(n),rep(mu,n))
    index <-c(0,cumsum(n))
    d <- rep(NA,k)
    for(i in 1:k) d[i] <- mean(data[(index[i]+1):index[i+1]])
    mu <- MLE_Williams(d,n)
    s2 <- sum((data-rep(d,n))^2)/nu
    t <- (mu[k]-d[1])/(sqrt(s2/n[1]+s2/n[k]))
    return(t)
  })
  ###
  p <- 1- sum(tnull<tobs)/(1+simanz)
  if(critical==TRUE) {
    return(list(p=p,q=quantile(tnull,qp)))
  }else{
    return(p)
  }
}
