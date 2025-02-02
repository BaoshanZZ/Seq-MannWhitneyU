source("Func_Tests.R")
##################################################################
############## code for categorical data #########################
##################################################################

### Calculate cumulative log-odds from class proportions
cum.logodds <- function(p) {
  nlev <- length(p)
  logodds <- head(sapply(1:nlev, function(i) log(sum(p[(i+1):length(p)])/sum(p[1:i]))),-1)
  return(logodds)
}

### Calculate class proportions from cumulative log-odds
cum.logodds.inv <- function(odds) {
  cum.p <- rev(c(1,exp(odds)/(1+exp(odds)),0))
  nlev <- length(cum.p)
  p <- sapply(2:nlev, function(i) cum.p[i] - cum.p[i-1])
  return(rev(p))
}

### quantities for asymptotic normality ################ 

calc.mean <- function(pX,pY,N,t,n.timepoints, rate.X) {
  
  ## Calculate mean of test statistic at time t (theta)
  k <- ifelse(length(pX) == length(pY), length(pX), NA)
  theta <- (0.5*pX[1])*pY[1]
  for (l in 2:k) {
    theta <- theta + (sum(pX[1:(l-1)]) + 0.5*pX[l])*pY[l]
  }
  
  
  varS_X <- ((0.5*pY[k])^2)*pX[k] - theta^2
  for (l in 1:(k-1)) {
    varS_X <- varS_X + ((sum(pY[(l+1):k]) + 0.5*pY[l])^2)*pX[l] 
  }
  varF_Y <- ((0.5*pX[1])^2)*pY[1] - theta^2
  for (l in 2:k) {
    varF_Y <- varF_Y + ((sum(pX[1:(l-1)]) + 0.5*pX[l])^2)*pY[l] 
  }
  varG_XY <- (0.25*pX[1])*pY[1] - theta^2
  for (l in 2:k) {
    varG_XY <- varG_XY + (sum(pX[1:(l-1)]) + 0.25*pX[l])*pY[l] 
  }
  N_k <- (t/n.timepoints)*N
  m <- rate.X*N_k
  n <- (1-rate.X)*N_k
  varUk <- (n-1)*(varS_X/(m*n)) + (m-1)*(varF_Y/(m*n)) + varG_XY/(m*n)
  
  ## Calculate test statistic at t
  Zk <- (theta-0.5)/sqrt(varUk)
  
  ## Create output data frame
  out <- data.frame(Zk = Zk, theta = theta, varUk = varUk)
  return(out)
}


### quantities for asymptotic normality with log odds transformation ################

calc.mean.log.ratio <- function(pX,pY,N,t,n.timepoints, rate.X) {
  
  ## Calculate mean of test statistic at time t (theta)
  k <- ifelse(length(pX) == length(pY), length(pX), NA)
  theta <- (0.5*pX[1])*pY[1]
  for (l in 2:k) {
    theta <- theta + (sum(pX[1:(l-1)]) + 0.5*pX[l])*pY[l]
  }
  
  
  varS_X <- ((0.5*pY[k])^2)*pX[k] - theta^2
  for (l in 1:(k-1)) {
    varS_X <- varS_X + ((sum(pY[(l+1):k]) + 0.5*pY[l])^2)*pX[l] 
  }
  varF_Y <- ((0.5*pX[1])^2)*pY[1] - theta^2
  for (l in 2:k) {
    varF_Y <- varF_Y + ((sum(pX[1:(l-1)]) + 0.5*pX[l])^2)*pY[l] 
  }
  varG_XY <- (0.25*pX[1])*pY[1] - theta^2
  for (l in 2:k) {
    varG_XY <- varG_XY + (sum(pX[1:(l-1)]) + 0.25*pX[l])*pY[l] 
  }
  N_k <- (t/n.timepoints)*N
  m <- rate.X*N_k
  n <- (1-rate.X)*N_k
  varUk <- (n-1)*(varS_X/(m*n)) + (m-1)*(varF_Y/(m*n)) + varG_XY/(m*n)
  
  varUk <- varUk / (theta^2 * (1-theta)^2)
  
  ## Calculate test statistic at t
  Zk <- (log(theta/(1-theta)))/sqrt(varUk)
  
  ## Create output data frame
  out <- data.frame(Zk = Zk, varUk = varUk)
  return(out)
}



####  critical boundaries ############################

ctz <- function(N, n.timepoints, ctl.p, or, rate.X, total.alpha, calc.mean) {
  B0 <- cum.logodds(ctl.p)
  B1 <- log(or)
  logodds <- B0 + B1
  trt.p <- cum.logodds.inv(logodds)
  
  
  
  varU <- calc.mean(ctl.p, trt.p, N, n.timepoints, n.timepoints, rate.X)$varUk
  t_k <- varU/(calc.mean(ctl.p, trt.p, N, c(1:n.timepoints), n.timepoints, rate.X)$varUk)
  
  # Create sequential covariance matrix of test statistic (for 1<i<j<n.timepoints, Cov_ij is ratio of variance at time j to variance at time i)
  sigma <- diag(n.timepoints)
  sigma[lower.tri(sigma)] <- unlist(sapply(1:(n.timepoints-1), 
                                           function(i) sqrt((calc.mean(ctl.p, trt.p, N, c((i+1):n.timepoints), n.timepoints, rate.X)$varUk)/
                                                              (calc.mean(ctl.p, trt.p, N, i, n.timepoints, rate.X)$varUk))))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  
  # Define critical values based on Kim and DeMets alpha spending function
  alpha <- c(total.alpha*(t_k[1])^2, diff(total.alpha*(t_k)^2))
  critical <- rep(0, n.timepoints)
  critical[1] <- qnorm(alpha[1], lower.tail = FALSE)
  
  if (n.timepoints > 1) {
    for (t in 2:n.timepoints) {
      fc <- function(ct) {
        pmvnorm(upper=critical[1:t-1], sigma = sigma[1:t-1,1:t-1]) - pmvnorm(upper = c(critical[1:t-1],ct), sigma = sigma[1:t,1:t]) - alpha[t]
      }
      critical[t] <- uniroot(fc,c(1,10))$root
    }
  }	 
  
  return(critical)
}



theor_TypeI <- function(N, n.timepoints, ctl.p, or, rate.X, total.alpha, calc.mean) {
  
  B0 <- cum.logodds(ctl.p)
  B1 <- log(or)
  logodds <- B0 + B1
  trt.p <- cum.logodds.inv(logodds)
  
  
  varU <- calc.mean(ctl.p, trt.p, N, n.timepoints, n.timepoints, rate.X)$varUk
  t_k <- varU/(calc.mean(ctl.p, trt.p, N, c(1:n.timepoints), n.timepoints, rate.X)$varUk)
  return(c(total.alpha * (t_k)^2))
  
  
}


theor_Power <- function(N, ctl.p, or, rate.X, n.timepoints, total.alpha, calc.mean, ctz)	{
  
  B0 <- cum.logodds(ctl.p)
  B1 <- log(or)
  logodds <- B0 + B1
  trt.p <- cum.logodds.inv(logodds)
  
  # Create sequential covariance matrix of test statistic (for 1<i<j<n.timepoints, Cov_ij is ratio of variance at time j to variance at time i)
  sigma <- diag(n.timepoints)
  sigma[lower.tri(sigma)] <- unlist(sapply(1:(n.timepoints-1), 
                                           function(i) sqrt((calc.mean(ctl.p, trt.p, N, c((i+1):n.timepoints), n.timepoints, rate.X)$varUk)/
                                                              (calc.mean(ctl.p, trt.p, N, i, n.timepoints, rate.X)$varUk))))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  
  alt.Z <- calc.mean(ctl.p, trt.p, N, c(1:n.timepoints), n.timepoints, rate.X)$Zk
  
  critical <- ctz(N, n.timepoints, ctl.p, or, rate.X, total.alpha, calc.mean)
  
  t.power <- cumsum(sapply(1:n.timepoints, function(t) {
    unname(pmvnorm(lower = c(rep(-Inf,t-1),critical[t]), upper = c(critical[1:t][-t],Inf), mean = alt.Z[1:t], 
                   sigma = sigma[1:t,1:t]))
  }))
  return(c(t.power))			 
  
}

theor_Power.log.ratio <- function(N, ctl.p, or, rate.X, n.timepoints, total.alpha, calc.mean.log.ratio, ctz)	{
  
  B0 <- cum.logodds(ctl.p)
  B1 <- log(or)
  logodds <- B0 + B1
  trt.p <- cum.logodds.inv(logodds)
  
  # Create sequential covariance matrix of test statistic (for 1<i<j<n.timepoints, Cov_ij is ratio of variance at time j to variance at time i)
  sigma <- diag(n.timepoints)
  sigma[lower.tri(sigma)] <- unlist(sapply(1:(n.timepoints-1), 
                                           function(i) sqrt((calc.mean.log.ratio(ctl.p, trt.p, N, c((i+1):n.timepoints), n.timepoints, rate.X)$varUk)/
                                                              (calc.mean.log.ratio(ctl.p, trt.p, N, i, n.timepoints, rate.X)$varUk))))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  
  alt.Z <- calc.mean.log.ratio(ctl.p, trt.p, N, c(1:n.timepoints), n.timepoints, rate.X)$Zk
  
  critical <- ctz(N, n.timepoints, ctl.p, or, rate.X, total.alpha, calc.mean.log.ratio)
  
  t.power <- cumsum(sapply(1:n.timepoints, function(t) {
    unname(pmvnorm(lower = c(rep(-Inf,t-1),critical[t]), upper = c(critical[1:t][-t],Inf), mean = alt.Z[1:t], 
                   sigma = sigma[1:t,1:t]))
  }))
  return(c(t.power))			 
  
}






### achieved rejection rates from different methods #######################

attained_cont <- function(N, ctl.p, or, rate.X, n.timepoints, B, total.alpha, ctz, calc.mean, calc.Zknew, wilcoxon.h0, brunner, brunner.log.ratio, ctz.adj){
  
  B0 <- cum.logodds(ctl.p)
  B1 <- log(or)
  logodds <- B0 + B1
  trt.p <- cum.logodds.inv(logodds)
  
  critical <- ctz(N, n.timepoints, ctl.p, or, rate.X, total.alpha, calc.mean)
  
  critical.log.ratio <- ctz(N, n.timepoints, ctl.p, or, rate.X, total.alpha, calc.mean.log.ratio)
  
  
  
  GS_results <- foreach(outind = 1 : B, .combine = cbind, .packages = c('splines', 'MASS'))%dorng% {
    
    X <- c(replicate((N*rate.X), sample(c(1 : length(ctl.p)),1,prob = ctl.p)))
    Y <- c(replicate((N*(1-rate.X)), sample(c(1 : length(trt.p)),1,prob = trt.p)))
    
    
    # Return Series of Zk
    
    seq.Z1 <- NULL
    seq.Z2 <- NULL
    seq.Z3 <- NULL
    seq.Z4 <- NULL
    
    for (ind in 1 : n.timepoints) {
      
      
      
      seq.Z1[ind] <- calc.Zknew(X, Y, N, n.timepoints, ind, rate.X)
      seq.Z2[ind] <- wilcoxon.h0(X, Y, N, n.timepoints, ind, rate.X)
      seq.Z3[ind] <- brunner(X, Y, N, n.timepoints, ind, rate.X) 
      seq.Z4[ind] <- brunner.log.ratio(X, Y, N, n.timepoints, ind, rate.X) 
      
    }
    if ( sum(seq.Z1 + seq.Z2 + seq.Z3 + seq.Z4) < 500000)  {
      
      critical.t <- ctz.adj(critical, X, Y, N, n.timepoints, rate.X)		
      
      reject1 <- rep(0, n.timepoints)
      reject1[Position(function(x) x == TRUE, seq.Z1 >= critical)] <- 1
      
      reject2 <- rep(0, n.timepoints)
      reject2[Position(function(x) x == TRUE, seq.Z1 >= critical.t)] <- 1
      
      reject3 <- rep(0, n.timepoints)
      reject3[Position(function(x) x == TRUE, seq.Z2 >= critical)] <- 1
      
      reject4 <- rep(0, n.timepoints)
      reject4[Position(function(x) x == TRUE, seq.Z3 >= critical)] <- 1
      
      reject5 <- rep(0, n.timepoints)
      reject5[Position(function(x) x == TRUE, seq.Z3 >= critical.t)] <- 1
      
      reject6 <- rep(0, n.timepoints)
      reject6[Position(function(x) x == TRUE, seq.Z4 >= critical.log.ratio)] <- 1
      
    } else {
      reject1 <- rep(0, n.timepoints)
      reject1[1] <- 1000
      
      reject2 <- rep(0, n.timepoints)
      reject2[1] <- 1000
      
      reject3 <- rep(0, n.timepoints)
      reject3[1] <- 1000
      
      reject4 <- rep(0, n.timepoints)
      reject4[1] <- 1000
      
      reject5 <- rep(0, n.timepoints)
      reject5[1] <- 1000
      
      reject6 <- rep(0, n.timepoints)
      reject6[1] <- 1000
      
    }
    
    c(reject1, reject2, reject3, reject4, reject5, reject6)
  }
  
  GS_colsums <- colSums(GS_results)
  GS_results <- GS_results[, (GS_colsums < 500)]
  
  
  return(rowMeans(GS_results))
  
  
}





B <- 1000

n.timepoints <- 3

total.alpha <- 0.01


ctl.p <- c(1/5, 1/5, 1/5, 1/5, 1/5)




ratexset <- c(1/2, 2/3)

orset <- c(1, 1.5, 2, 2.5, 3)

nset <- rev(c(108, 108*2, 108*4, 108*6))

res <- c("ctl.p <- c(1/5, 1/5, 1/5, 1/5, 1/5)","","","")

for (ratexind in 1 : 1) {
  for (orind in 1 : 5) {
    for (nind in 1 : 4) {
      
      or <- orset[orind]
      
      N <- nset[nind]
      
      rate.X <- ratexset[ratexind]
      
      cat("\n", "or = ", or)  
      
      cat("\n", "N = ", N)
      
      cat("\n", "rate.X = ", rate.X)
      
      type1 <- theor_TypeI(N, n.timepoints, ctl.p, or, rate.X, total.alpha, calc.mean) 
      
      type1 <- c("type1_error", type1)
      
      boundary <- ctz(N, n.timepoints, ctl.p, or, rate.X, total.alpha, calc.mean)
      
      boundary <- c("critical boundaries", boundary)
      
      th_p <- theor_Power(N, ctl.p, or, rate.X, n.timepoints, total.alpha, calc.mean, ctz)
      
      th_p <- c("theory power", th_p)
      
      th_p_logratio <- theor_Power.log.ratio(N, ctl.p, or, rate.X, n.timepoints, total.alpha, calc.mean.log.ratio, ctz)
      
      th_p_logratio <- c("theory power logratio", th_p_logratio)
      
      out <- round(attained_cont(N, ctl.p, or, rate.X, n.timepoints, B, total.alpha, ctz, calc.mean, calc.Zknew, wilcoxon.h0, brunner, brunner.log.ratio, ctz.adj), 5)
      
      #res <- c("","","","")
      
      res1 <- cumsum(out[1 : n.timepoints])
      
      res1 <- c("est", res1)
      
      res2 <- cumsum(out[(n.timepoints + 1) : (2 * n.timepoints)])
      
      res2 <- c("est.t", res2)
      
      res3 <- cumsum(out[(2 * n.timepoints + 1) : (3 * n.timepoints)])
      
      res3 <- c("wilcoxon.null", res3)
      
      res4 <- cumsum(out[(3 * n.timepoints + 1) : (4 * n.timepoints)])
      
      res4 <- c("brunner", res4)
      
      res5 <- cumsum(out[(4 * n.timepoints + 1) : (5 * n.timepoints)])
      
      res5 <- c("brunner.t", res5)
      
      res6 <- cumsum(out[(5 * n.timepoints + 1) : (6 * n.timepoints)])
      
      res6 <- c("log.ratio", res6)
      
      setting <- c(paste0("Total size=", N), paste0("or=", or), paste0("rate.X=", rate.X), "")
      
      res <- rbind(res, setting, type1, boundary, th_p, th_p_logratio, res1, res2, res3, res4, res5, res6)
      
      
    }
  }
}

res <- as.data.frame(res)


