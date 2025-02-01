
######################################################################################
###################### code for Poissson data#########################################
######################################################################################

# true thetas for MRs = 1, 1.5, 2, 2.5, 3
Utrue <- c(0.5, 0.6176168, 0.7114280, 0.7843031, 0.8405808)

# three parts of to get variances for different MRs = 1, 1.5, 2, 2.5, 3
vpart1 <- c(0.07454709, 0.06488253, 0.05133555, 0.03811935, 0.02685471)
vpart2 <- c(0.07450799, 0.07753859, 0.06985690, 0.05781059, 0.04509109)
vpart3 <- c(0.2439110, 0.2470775, 0.2247510, 0.1912077, 0.1551858)


### quantities for asymptotic normality ################ 

calc.mean <- function(Utrue, vpart1, vpart2, vpart3, mrind, N, t, n.timepoints, rate.X) {
  
  theta <- Utrue[mrind]
  
  N_k <- (t/n.timepoints)*N
  m <- rate.X*N_k
  n <- (1-rate.X)*N_k
  varUk <- (n-1)*(vpart1[mrind]/(m*n)) + (m-1)*(vpart2[mrind]/(m*n)) + vpart3[mrind]/(m*n)
  
  ## Calculate test statistic at t
  Zk <- (theta-0.5)/sqrt(varUk)
  
  ## Create output data frame
  out <- data.frame(Zk = Zk, varUk = varUk)
  return(out)
}

### quantities for asymptotic normality with log odds transformation ################

calc.mean.log.ratio <- function(Utrue, vpart1, vpart2, vpart3, mrind, N, t, n.timepoints, rate.X) {
  
  
  theta <- Utrue[mrind]
  
  N_k <- (t/n.timepoints)*N
  m <- rate.X*N_k
  n <- (1-rate.X)*N_k
  varUk <- (n-1)*(vpart1[mrind]/(m*n)) + (m-1)*(vpart2[mrind]/(m*n)) + vpart3[mrind]/(m*n)
  
  varUk <- varUk / (theta^2 * (1-theta)^2)
  
  ## Calculate test statistic at t
  Zk <- (log(theta/(1-theta)))/sqrt(varUk)
  
  ## Create output data frame
  out <- data.frame(Zk = Zk, varUk = varUk)
  return(out)
}

####  critical boundaries ############################

ctz <- function(Utrue, vpart1, vpart2, vpart3, N, n.timepoints, rate.X, total.alpha, calc.mean, mrind) {
  
  
  
  varU <- calc.mean(Utrue, vpart1, vpart2, vpart3, mrind, N,n.timepoints,n.timepoints, rate.X)$varUk
  t_k <- varU/(calc.mean(Utrue, vpart1, vpart2, vpart3, mrind, N, c(1:n.timepoints), n.timepoints, rate.X)$varUk)
  
  # Create sequential covariance matrix of test statistic (for 1<i<j<n.timepoints, Cov_ij is ratio of variance at time j to variance at time i)
  sigma <- diag(n.timepoints)
  sigma[lower.tri(sigma)] <- unlist(sapply(1:(n.timepoints-1), 
                                           function(i) sqrt((calc.mean(Utrue, vpart1, vpart2, vpart3, mrind, N, c((i+1):n.timepoints), n.timepoints, rate.X)$varUk)/
                                                              (calc.mean(Utrue, vpart1, vpart2, vpart3, mrind, N, i, n.timepoints, rate.X)$varUk))))
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


theor_TypeI <- function(N, n.timepoints, Utrue, vpart1, vpart2, vpart3, rate.X, total.alpha, calc.mean, mrind) {
  
  
  varU <- calc.mean(Utrue, vpart1, vpart2, vpart3, mrind, N, n.timepoints, n.timepoints, rate.X)$varUk
  t_k <- varU/(calc.mean(Utrue, vpart1, vpart2, vpart3, mrind, N, c(1:n.timepoints), n.timepoints, rate.X)$varUk)
  return(c(total.alpha * (t_k)^2))
  
  
}


theor_Power <- function(N, Utrue, vpart1, vpart2, vpart3, mrind, rate.X, n.timepoints, total.alpha, calc.mean, ctz)	{
  
  # Create sequential covariance matrix of test statistic (for 1<i<j<n.timepoints, Cov_ij is ratio of variance at time j to variance at time i)
  sigma <- diag(n.timepoints)
  sigma[lower.tri(sigma)] <- unlist(sapply(1:(n.timepoints-1), 
                                           function(i) sqrt((calc.mean(Utrue, vpart1, vpart2, vpart3, mrind, N, c((i+1):n.timepoints), n.timepoints, rate.X)$varUk)/
                                                              (calc.mean(Utrue, vpart1, vpart2, vpart3, mrind, N, i, n.timepoints, rate.X)$varUk))))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  
  alt.Z <- calc.mean(Utrue, vpart1, vpart2, vpart3, mrind, N, c(1:n.timepoints), n.timepoints, rate.X)$Zk
  
  critical <- ctz(Utrue, vpart1, vpart2, vpart3, N, n.timepoints, rate.X, total.alpha, calc.mean, mrind)
  
  t.power <- cumsum(sapply(1:n.timepoints, function(t) {
    unname(pmvnorm(lower = c(rep(-Inf,t-1),critical[t]), upper = c(critical[1:t][-t],Inf), mean = alt.Z[1:t], 
                   sigma = sigma[1:t,1:t]))
  }))
  return(c(t.power))			 
  
}

### achieved rejection rates from different methods #######################

attained_cont <- function(N, Utrue, vpart1, vpart2, vpart3, mrset, mrind, rate.X, n.timepoints, B, total.alpha, ctz, calc.mean, calc.Zknew, wilcoxon.h0, brunner, brunner.log.ratio, ctz.adj){
  
  
  critical <- ctz(Utrue, vpart1, vpart2, vpart3, N, n.timepoints, rate.X, total.alpha, calc.mean, mrind)
  
  critical.log.ratio <-ctz(Utrue, vpart1, vpart2, vpart3, N, n.timepoints, rate.X, total.alpha, calc.mean.log.ratio, mrind)
  
  GS_results <- foreach(outind = 1 : B, .combine = cbind, .packages = c('splines', 'MASS'))%dorng% {
    
    X <- rpois((N*rate.X), 1)
    
    Y <- rpois((N*(1-rate.X)), (1*mrset[mrind]))
    
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



ratexset <- c(1/2, 2/3)

mrset <- c(1, 1.5, 2, 2.5, 3)

nset <- rev(c(108, 108*2, 108*4, 108*6))

res <- c("mean.X <- 1","","","")

for (ratexind in 1 : 1) {
  for (mrind in 1 : 5) {
    for (nind in 1 : 4) {
      
      
      MR <- mrset[mrind]
      
      N <- nset[nind]
      
      rate.X <- ratexset[ratexind]
      
      cat("\n", "MR = ", MR)  
      
      cat("\n", "N = ", N)
      
      cat("\n", "rate.X = ", rate.X)
      
      type1 <- theor_TypeI(N, n.timepoints, Utrue, vpart1, vpart2, vpart3, rate.X, total.alpha, calc.mean, mrind)
      
      type1 <- c("type1_error", type1)
      
      boundary <- ctz(Utrue, vpart1, vpart2, vpart3, N, n.timepoints, rate.X, total.alpha, calc.mean, mrind)
      
      boundary <- c("critical boundaries", boundary)
      
      th_p <- theor_Power(N, Utrue, vpart1, vpart2, vpart3, mrind, rate.X, n.timepoints, total.alpha, calc.mean, ctz)
      
      th_p <- c("theory power", th_p)
      
      th_p_logratio <- theor_Power(N, Utrue, vpart1, vpart2, vpart3, mrind, rate.X, n.timepoints, total.alpha, calc.mean.log.ratio, ctz)
      
      th_p_logratio <- c("theory power logratio", th_p_logratio)
      
      
      
      
      out <- round(attained_cont(N, Utrue, vpart1, vpart2, vpart3, mrset, mrind, rate.X, n.timepoints, B, total.alpha, ctz, calc.mean, calc.Zknew, wilcoxon.h0, brunner, brunner.log.ratio, ctz.adj), 5)
      
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
      
      setting <- c(paste0("Total size=", N), paste0("MR=", MR), paste0("rate.X=", rate.X),"")
      
      res <- rbind(res, setting, type1, boundary, th_p, th_p_logratio, res1, res2, res3, res4, res5, res6)
      
      
    }
  }
}

res <- as.data.frame(res)