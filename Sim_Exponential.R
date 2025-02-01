source("Func_Main.R")
######################################################################################
###################### code for exponential data########################################
######################################################################################

### quantities for asymptotic normality ################ 

calc.mean <- function(N, n.timepoints, t, rate.X, lambda.x, MR) {
  
  lambda.y <- lambda.x/MR
  
  mean.u.alt <- lambda.x / (lambda.x + lambda.y)
  
  mean.var.survx <- lambda.x / (lambda.x + 2*lambda.y) - (mean.u.alt)^2
  
  mean.var.cdfy <- lambda.y / (2*lambda.x + lambda.y) - (1 - mean.u.alt)^2
  
  
  N_k <- (t/n.timepoints)*N
  m <- rate.X*N_k
  n <- (1-rate.X)*N_k
  
  varUk <- (n-1) * (mean.var.survx/(m * n)) + (m-1)*(mean.var.cdfy/(m * n)) + (mean.u.alt - mean.u.alt^2)/(m * n)
  
  ## Calculate test statistic at t
  Zk <- (mean.u.alt - 0.5)/sqrt(varUk)
  
  ## Create output data frame
  out <- data.frame(Zk = Zk, varUk = varUk)
  return(out)
}

### quantities for asymptotic normality with log odds transformation ################

calc.mean.log.ratio <- function(N, n.timepoints, t, rate.X, lambda.x, MR) {
  
  lambda.y <- lambda.x/MR
  
  mean.u.alt <- lambda.x / (lambda.x + lambda.y)
  
  mean.var.survx <- lambda.x / (lambda.x + 2*lambda.y) - (mean.u.alt)^2
  
  mean.var.cdfy <- lambda.y / (2*lambda.x + lambda.y) - (1 - mean.u.alt)^2
  
  
  N_k <- (t/n.timepoints)*N
  m <- rate.X*N_k
  n <- (1-rate.X)*N_k
  
  varUk <- (n-1) * (mean.var.survx/(m * n)) + (m-1)*(mean.var.cdfy/(m * n)) + (mean.u.alt - mean.u.alt^2)/(m * n)
  
  varUk <- varUk / (mean.u.alt^2 * (1 - mean.u.alt)^2)
  
  ## Calculate test statistic at t
  Zk <- (log(mean.u.alt/(1-mean.u.alt)))/sqrt(varUk)
  
  ## Create output data frame
  out <- data.frame(Zk = Zk, varUk = varUk)
  return(out)
}

####  critical boundaries ############################

ctz <- function(N, n.timepoints, total.alpha, rate.X, lambda.x, MR, calc.mean) {
  #### variance for U under null for continuous data############
  varU <- calc.mean(N, n.timepoints, n.timepoints, rate.X, lambda.x, MR)$varUk
  t_k <- varU/(calc.mean(N, n.timepoints, c(1:n.timepoints), rate.X, lambda.x, MR)$varUk)
  
  sigma <- diag(n.timepoints)
  sigma[lower.tri(sigma)] <- unlist(sapply(1:(n.timepoints-1), 
                                           function(i) sqrt((calc.mean(N, n.timepoints, c((i+1) : n.timepoints), rate.X, lambda.x, MR)$varUk)/
                                                              (calc.mean(N, n.timepoints, i, rate.X, lambda.x, MR)$varUk))))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  
  ### Define critical values based on Kim and DeMets alpha spending function
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


theor_Power <- function(N, n.timepoints, total.alpha, rate.X, lambda.x, MR, calc.mean)	{
  
  
  sigma <- diag(n.timepoints)
  sigma[lower.tri(sigma)] <- unlist(sapply(1:(n.timepoints-1), 
                                           function(i) sqrt((calc.mean(N, n.timepoints, c((i+1) : n.timepoints), rate.X, lambda.x, MR)$varUk)/
                                                              (calc.mean(N, n.timepoints, i, rate.X, lambda.x, MR)$varUk))))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  
  alt.Z <- calc.mean(N, n.timepoints, c(1 : n.timepoints), rate.X, lambda.x, MR)$Zk
  
  critical <- ctz(N, n.timepoints, total.alpha, rate.X, lambda.x, MR, calc.mean)
  
  t.power <- cumsum(sapply(1:n.timepoints, function(t) {
    unname(pmvnorm(lower = c(rep(-Inf,t-1),critical[t]), upper = c(critical[1:t][-t],Inf), mean = alt.Z[1:t], 
                   sigma = sigma[1:t,1:t]))
  }))
  return(c(t.power))			 
  
}

theor_TypeI <- function(N, n.timepoints, rate.X, lambda.x, MR, total.alpha, calc.mean) {
  varU <- calc.mean(N, n.timepoints, n.timepoints, rate.X, lambda.x, MR)$varUk
  t_k <- varU/(calc.mean(N, n.timepoints, c(1:n.timepoints), rate.X, lambda.x, MR)$varUk)
  return(c(total.alpha * (t_k)^2))
}

### achieved rejection rates from different methods #######################

attained_cont <- function(N, MR, rate.X, n.timepoints, B, total.alpha,  ctz, calc.Zknew, wilcoxon.h0, brunner, brunner.log.ratio, lambda.x, calc.mean, calc.mean.log.ratio, ctz.adj){
  
  
  
  critical <- ctz(N, n.timepoints, total.alpha, rate.X, lambda.x, MR, calc.mean)
  critical.log.ratio <- ctz(N, n.timepoints, total.alpha, rate.X, lambda.x, MR, calc.mean.log.ratio)
  
  GS_results <- foreach(ind = 1 : B, .combine = cbind, .packages = c('splines', 'MASS'))%dorng% {
    
    X <- rexp((N*rate.X), lambda.x)
    
    Y<-rexp((N*(1-rate.X)), lambda.x/MR)
    
    
    
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





B <- 100000

n.timepoints <- 3

total.alpha <- 0.01

lambda.x <- 1

mrset <- c(1, 1.5, 2, 2.5, 3)

nset <- rev(c(108, 108*2, 108*4, 108*6))

ratexset <- c(1/2, 2/3)

res <- c("lambda.x <- 1","","","")

for (ratexind in 1 : 1) {
  for (mrind in 1 : 5) {
    for (nind in 1 : 4) {
      
      MR <- mrset[mrind]
      
      N <- nset[nind]
      
      rate.X <- ratexset[ratexind]
      
      cat("\n", "MR = ", MR)
      
      cat("\n", "N = ", N)
      
      cat("\n", "rate.X = ", rate.X)
      
      type1 <- theor_TypeI(N, n.timepoints, rate.X, lambda.x, MR, total.alpha, calc.mean)
      
      type1 <- c("type1_error", type1)
      
      boundary <- ctz(N, n.timepoints, total.alpha, rate.X, lambda.x, MR, calc.mean)
      
      boundary <- c("critical boundaries", boundary)
      
      th_p <- theor_Power(N, n.timepoints, total.alpha, rate.X, lambda.x, MR, calc.mean)
      
      th_p <- c("theory power", th_p)
      
      th_p_logratio <- theor_Power(N, n.timepoints, total.alpha, rate.X, lambda.x, MR, calc.mean.log.ratio)
      
      th_p_logratio <- c("theory power logratio", th_p_logratio)
      
      out <- round(attained_cont(N, MR, rate.X, n.timepoints, B, total.alpha,  ctz, calc.Zknew, wilcoxon.h0, brunner, brunner.log.ratio, lambda.x, calc.mean, calc.mean.log.ratio, ctz.adj), 5)
      
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
      
      setting <- c(paste0("Total size=", N), paste0("MR=", MR), paste0("rate.X=", rate.X), "")
      
      res <- rbind(res, setting, type1, boundary, th_p, th_p_logratio, res1, res2, res3, res4, res5, res6)
      
    }
  }
}	

res <- as.data.frame(res)
