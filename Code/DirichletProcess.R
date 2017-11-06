rm(list = ls())
setwd("~/Documents/Research/Virus/Code")
source("allfcns.R")


logBeta <- function(alphas){
  lnum <- sapply(alphas, function(x){return(lgamma(x))})
  lden <- lgamma(sum(alphas))
  return(sum(lnum) - lden)
}

logstari <- function(Cs, Ys, idx, newci){
  J <- nrow(Ys)
  subYs <- Ys[ , -idx]
  subCs <- Cs[-idx]
  cidx <- which(subCs == newci)
  
  Yjs <- Ys[ , cidx]
  
  if(length(cidx) == 1){
    alphas <- Yjs + 1/J^2
  }else{
    alphas <- rowSums(Yjs) + 1/J^2
  }

  lstari <- -logBeta(alphas) + logBeta(alphas + Ys[ , idx])
  
  return(lstari)
}

logdiamondi <- function(Yi){
  J <- length(Yi)
  alphas <- rep(1/J^2, J)
  
  ldiamondi <- -logBeta(alphas) + logBeta(alphas + Yi)
  
  return(ldiamondi)
}

newciProb <- function(Cs, Ys, idx, alpha0){
  subYs <- Ys[ , -idx]
  subCs <- Cs[-idx]
  Yi <- Ys[ , idx]
  
  subCs_unique <- unique(subCs)
  K <- length(subCs_unique)
  
  lnum0 <- log(alpha0) + logdiamondi(Yi)
  nums <- sapply(subCs_unique, function(x){
    numc <- log(sum(subCs == x)) + logstari(Cs, Ys, idx, x) - lnum0
    return(exp(numc))
  })
  nums[K + 1] = 1
  
  prob <- nums/sum(nums)

  return(prob)
}

DPmixtureLP <- function(Cs, Ys, alpha0){
  Cs_unique <- unique(Cs)
  K <- length(Cs_unique)
  J <- nrow(Ys)
  
  Part1 <- -K * lgamma(alpha0/K)
  
  Part2 <- sapply(Cs_unique, function(x){
    xidx <- which(Cs == x)
    
    if(length(xidx) == 1){
      alphas <- Ys[ , xidx] + 1/J^2
    }else{
      alphas <- rowSums(Ys[ , xidx]) + 1/J^2
    }
   
    y <- logBeta(alphas)
    return(y)
    })
  Part2 <- sum(Part2)
  
  Part3 <- sapply(Cs_unique, function(x){
    y <- lgamma(sum(Cs == x) + alpha0/K)
    return(y)
  })
  Part3 <- sum(Part3)
  
  return(Part1 + Part2 + Part3)
}

DPmixtureGeweke <- function(Ys, initialC, s, a, bk, alpha0, savefile = NULL){
  n <- length(initialC)
  K <- length(unique(initialC))
  
  t <- 0
  C0 <- initialC
  LPs <- rep(NA, s)
  Cs <- matrix(NA, nrow = s, ncol = n)
  absZ <- 10; burn <- K * bk; z <- 0
  maxtotal <- burn + 1000000 * K
  
  while(absZ >= 2){
    for(idx in 1:n){
      ci <- C0[idx]
      probi <- newciProb(C0, Ys, idx, alpha0)
      cumprobi <- cumsum(probi)
      u <- runif(1)
      C0[idx] = sum(u > cumprobi) + 1
      K <- length(unique(C0))

      if(max(C0) > K){
        C0 <- sapply(C0, function(x){
          if(x > ci){return(x - 1)}
          else{return(x)}
        })
      }
      
      t <- t + 1
      
      if(t > burn){
        tb <- t - burn
        sidx <- tb %% s
        if(sidx == 0){sidx <- s}
        LPs[sidx] = DPmixtureLP(C0, Ys, alpha0) 
        Cs[sidx, ] <- C0
        
        if(sidx == s){
          print(t)
          MClp <- mcmc(LPs)
          geweke <- geweke.diag2(MClp)
          z <- geweke$z
          if(is.na(z)){
            absZ <- 0
          }else if(z == "Inf"){
            absZ <- 10
          }else{
            absZ <- abs(z)
          }
        }
        
        C <- convC(Cs, LPs)
        MC <- list(RunningC = Cs, LP = LPs, t = t, zscore = z, C = C)
        
        if(!is.null(savefile)){save(MC, file = savefile)}
      }
     
      if(t == maxtotal){
        absZ <- 0
        print(paste0("Stopped after ", maxtotal, " iterations"))
      }   
    }
    
  }
  
  C <- convC(Cs, LPs)
  MC <- list(RunningC = Cs, LP = LPs, t = t, zscore = z, C = C)
  if(!is.null(savefile)){save(MC, file = savefile)}
  return(MC)  
  
}

###############################################################################






