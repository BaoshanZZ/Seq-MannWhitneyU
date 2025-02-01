library(mvtnorm)
library(splines)
library(MASS)
library(foreach)
library(doRNG)
library(doParallel)
ncores<-4
registerDoParallel(cores=ncores)
library(readxl)
library(openxlsx)

### our proposed  test  #############

calc.Zknew <- function(X, Y, N, n.timepoints, t, rate.X) {
  
  m_k <- (t / n.timepoints) * N * rate.X
  n_k <- (t / n.timepoints) * N * (1-rate.X) 
  ### Calculate Uk
  X_k <- X[1 : m_k]
  Y_k <- Y[1 : n_k]
  
  Uk<-mean(sapply(Y_k, function(Y_k) length(X_k[X_k < Y_k]) / m_k + (length(X_k[X_k == Y_k]) / m_k)/2 ))
  
  survx.vec<-sapply(X_k, function(X_k) length(Y_k[Y_k > X_k]) / n_k + (length(Y_k[X_k == Y_k]) / n_k)/2)
  var.survx<-var(survx.vec)
  # CDF
  cdfy.vec<-sapply(Y_k, function(Y_k) length(X_k[X_k < Y_k]) / m_k + (length(X_k[X_k == Y_k]) / m_k)/2)
  var.cdfy<-var(cdfy.vec)
  
  meank <- mean(sapply(Y_k, function(Y_k) length(X_k[X_k < Y_k]) / m_k + (length(X_k[X_k == Y_k]) / m_k)/4 ))
  # Ind
  #ind.var <- var(X_k < Y_k)
  ind.var <- meank - meank^2
  
  # Var_Uk
  Var.Uk <- ((n_k-1) * var.survx + (m_k - 1) * var.cdfy + ind.var) / (n_k * m_k)
  
  if (abs(Var.Uk)<10^(-8)) {
    Zk <- 10^8
  } else {
    Zk <- (Uk - 0.5) / sqrt(Var.Uk)
  }
  
  return(Zk)
  
}

### Wilcoxon test with variance assuming X=Y ################

wilcoxon.h0 <- function(X, Y, N, n.timepoints, t, rate.X) {
  
  m_k <- (t / n.timepoints) * N * rate.X
  n_k <- (t / n.timepoints) * N * (1-rate.X) 
  ### Calculate Uk
  X_k <- X[1 : m_k]
  Y_k <- Y[1 : n_k]
  
  Uk<-mean(sapply(Y_k, function(Y_k) length(X_k[X_k < Y_k]) / m_k + (length(X_k[X_k == Y_k]) / m_k)/2 ))
  
  rankall.k <- rank(c(X_k, Y_k), ties.method = "average")
  
  sigmasq <- sum((rankall.k - (m_k + n_k +1)/2)^2)/(m_k + n_k -1)
  
  Var.Uk <- sigmasq/((m_k + n_k) * m_k * n_k)
  
  
  if (abs(Var.Uk)<10^(-8)) {
    Zk <- 10^8
  } else {
    Zk <- (Uk - 0.5) / sqrt(Var.Uk)
  }
  
  
  return(Zk)
  
}

### Brunner and Munzel test ##############

brunner <- function(X, Y, N, n.timepoints, t, rate.X) {
  
  m_k <- (t / n.timepoints) * N * rate.X
  n_k <- (t / n.timepoints) * N * (1-rate.X) 
  ### Calculate Uk
  X_k <- X[1 : m_k]
  Y_k <- Y[1 : n_k]
  
  Uk<-mean(sapply(Y_k, function(Y_k) length(X_k[X_k < Y_k]) / m_k + (length(X_k[X_k == Y_k]) / m_k)/2 ))
  
  rankall.k <- rank(c(X_k, Y_k), ties.method = "average")
  
  rank.Xk <- rank(X_k, ties.method = "average")
  
  rank.Yk <- rank(Y_k, ties.method = "average")
  
  sigmasq.X <- sum((rankall.k[1 : m_k] - rank.Xk - mean(rankall.k[1 : m_k]) + (m_k + 1)/2)^2)/((m_k-1)*n_k^2)
  
  sigmasq.Y <- sum((rankall.k[(m_k + 1) : (m_k + n_k)] - rank.Yk - mean(rankall.k[(m_k + 1) : (m_k + n_k)]) + (n_k + 1)/2)^2)/((n_k-1)*m_k^2)
  
  Var.Uk <- sigmasq.X/m_k + sigmasq.Y/n_k
  
  if (abs(Var.Uk)<10^(-8)) {
    Zk <- 10^8
  } else {
    Zk <- (Uk - 0.5) / sqrt(Var.Uk)
  }
  
  return(Zk)
  
}

### Brunner and Munzel with log odds transformation ###############

brunner.log.ratio <- function(X, Y, N, n.timepoints, t, rate.X) {
  
  
  
  m_k <- (t / n.timepoints) * N * rate.X
  n_k <- (t / n.timepoints) * N * (1-rate.X) 
  ### Calculate Uk
  X_k <- X[1 : m_k]
  Y_k <- Y[1 : n_k]
  
  Uk<-mean(sapply(Y_k, function(Y_k) length(X_k[X_k < Y_k]) / m_k + (length(X_k[X_k == Y_k]) / m_k)/2 ))
  
  rankall.k <- rank(c(X_k, Y_k), ties.method = "average")
  
  rank.Xk <- rank(X_k, ties.method = "average")
  
  rank.Yk <- rank(Y_k, ties.method = "average")
  
  sigmasq.X <- sum((rankall.k[1 : m_k] - rank.Xk - mean(rankall.k[1 : m_k]) + (m_k + 1)/2)^2)/((m_k-1)*n_k^2)
  
  sigmasq.Y <- sum((rankall.k[(m_k + 1) : (m_k + n_k)] - rank.Yk - mean(rankall.k[(m_k + 1) : (m_k + n_k)]) + (n_k + 1)/2)^2)/((n_k-1)*m_k^2)
  
  Var.Uk <- sigmasq.X/m_k + sigmasq.Y/n_k
  
  if (abs(Uk)<10^(-8) | abs(1-Uk)<10^(-8) | abs(Var.Uk)<10^(-8)) {
    Zk<-10^8
  } else {
    Var.Uk <- Var.Uk/(Uk^2 * (1 - Uk)^2)
    
    
    Zk <- (log(Uk/(1 - Uk)))/sqrt(Var.Uk)
  }
  return(Zk)
  
  
}

#### t-adjusted critical boundaries ##################

ctz.adj <- function(critical, X, Y, N, n.timepoints, rate.X) {
  
  res <- NULL
  
  for (ind in 1 : n.timepoints) {
    m_k <- (ind / n.timepoints) * N * rate.X
    n_k <- (ind / n.timepoints) * N * (1-rate.X) 
    ### Calculate Uk
    X_k <- X[1 : m_k]
    Y_k <- Y[1 : n_k]
    
    Uk<-mean(sapply(Y_k, function(Y_k) length(X_k[X_k < Y_k]) / m_k + (length(X_k[X_k == Y_k]) / m_k)/2 ))
    
    rankall.k <- rank(c(X_k, Y_k), ties.method = "average")
    
    rank.Xk <- rank(X_k, ties.method = "average")
    
    rank.Yk <- rank(Y_k, ties.method = "average")
    
    sigmasq.X <- sum((rankall.k[1 : m_k] - rank.Xk - mean(rankall.k[1 : m_k]) + (m_k + 1)/2)^2)/((m_k-1)*n_k^2)
    
    sigmasq.Y <- sum((rankall.k[(m_k + 1) : (m_k + n_k)] - rank.Yk - mean(rankall.k[(m_k + 1) : (m_k + n_k)]) + (n_k + 1)/2)^2)/((n_k-1)*m_k^2)
    
    Var.Uk <- sigmasq.X/m_k + sigmasq.Y/n_k
    
    df <- Var.Uk^2/(sigmasq.X^2/(m_k^2*(m_k-1)) + sigmasq.Y^2/(n_k^2*(n_k-1)))
    
    res[ind] <- qt(pnorm(critical[ind], 0, 1), df=df)
  }	
  
  return(res)
  
}	
