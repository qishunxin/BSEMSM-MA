rm(list=ls());
library(mvtnorm)#2t+1
#install.packages("HI");                                   
library(MASS)
library(HI)
library(coda)
library(MCMCpack)
##defined constants## 
rep_num <- 100
n <- 500
z <- 4
s <- 2
p <- 4
q <- 2
G <- 5
###############################################################################

n.mcmc = 20000
n.kept = 10000
n.burn = n.mcmc - n.kept
n.thin = 1
sig.mh = 3
sig.mh.g = 0.5

PI = 0   #for 1 latent mediator case
alpha = 0 

####### total number of acceptance for each variable in MH algorithm
##########################
# true parameter values #

Delta.true <- matrix(c(0, 0.2, 0, 0.5, 0, 0.5, 0.5,
                       0.2, 0, 0.5, 0, 0.5, 0, 0.5), 
              nrow = q, ncol = (1 + z + s), byrow = TRUE)

L.se.true <- cbind(Delta.true)
psd.true <- rep(0.3, q)  

B.true <- matrix(0, nrow = p, ncol = q) 
B.true[,1] <- c(1, 0.9, 0, 0)
B.true[,2] <- c(0, 0, 1, 0.9)

L.me.true <- B.true
psi.true <- rep(0.2, p)  
 
# for calculation

gam.z.true.a <- c(0.1, 0.3, 0.2, 0.3)
gam.s.true.a <- c(0.2, 0.3)
gam.m.true.a <- c(0.3, 0.2)
gam.true.a <- c(gam.z.true.a, gam.s.true.a, gam.m.true.a)
lam0.true.a <- rep(1, G)

gam.z.true.b <- c(0.3, 0.1, 0.3, 0.1)
gam.s.true.b <- c(0.3, 0.2)
gam.m.true.b <- c(0.2, 0.3)
gam.true.b <- c(gam.z.true.b, gam.s.true.b, gam.m.true.b)
lam0.true.b <- rep(1, G)

gam.z.true.c <- c(0.2, 0.2, 0.5, 0.3)
gam.s.true.c <- c(0.3, 0.4)
gam.m.true.c <- c(0.3, 0.4)
gam.true.c <- c(gam.z.true.c, gam.s.true.c, gam.m.true.c)
lam0.true.c <- rep(1, G)


osigmaT= 0.5;
omegaT = rnorm(n, 0, osigmaT); 
alpha1T = 0.4;
alpha2T = 0.5;

########## Model Identification #########
### for the loading matrix of CFA model

Id.B <- matrix(c(
  0,0,
  1,0,
  0,0,
  0,1
), nrow = p, ncol = q, byrow = TRUE )

Id.B <- (Id.B > 0)
## free loadings in each row of the loading matrix B
n.b.row <- rowSums(Id.B)  
n.B <- sum(Id.B)

###for the measurment equation
Id.me <- Id.B
n.me.row <- rowSums(Id.me)  
n.me <- sum(Id.me)
 
Id.psi <- rep(1, p)
Id.psi <- as.logical(Id.psi)
n.psi <- sum(Id.psi)
 
Id.Delta <- matrix(rep(1, q*(1 + s + z)), nrow  = q)
Id.Delta <- (Id.Delta > 0)

## for the structual equation
Id.se <- cbind(Id.Delta)
n.se.row <- rowSums(Id.se)  
n.se <- sum(Id.se)

Id.psd <- rep(1, q)
Id.psd <- as.logical(Id.psd)
n.psd <- sum(Id.psd)

### for the PH model
Id.gam.a <- as.logical(rep(1, s + z + q))
n.gam.a <- sum(Id.gam.a)
Id.gam.b <- as.logical(rep(1, s + z + q))
n.gam.b <- sum(Id.gam.b)
Id.gam.c <- as.logical(rep(1, s + z + q))
n.gam.c <- sum(Id.gam.c)

##########prior setting

rho.scale <- 3
rho0 <- rho.scale + 1
Delta0 <- matrix(rep(0, 1 + z + s), nrow = q, ncol = (1 + z + s), byrow = TRUE)

sig.delta <- rep(0.001,n.se)         # prior precision for β

L.se0 <- cbind(Delta0)
 
alpha.psd <- 9      # α0_epsilon
beta.psd <- 4       # β0_epsilon
 
B0 <- matrix(c(
  0,0,
  0,0,
  0,0,
  0,0
), nrow = p, ncol = q, byrow = TRUE)
sig.b <- 0.001                     

L.me0 <- B0

alpha.psi <- 9      # α0_ζ
beta.psi <- 4       # β0_ζ

gam0.a <- rep(0, s + z + q)
sig.gamma.a <- rep(0.001,n.gam.a)   

gam0.b <- rep(0, s + z + q)
sig.gamma.b <- rep(0.001,n.gam.b)         
 
gam0.c <- rep(0, s + z + q)
sig.gamma.c <- rep(0.001,n.gam.c)     

#### For the piecewise constant λ0
#### we choose G = 5, thus λ = (λ1, λ2，λ3，λ4，λ5), corresponding to time interval
####                          [d1 = 0,d2],(d2,d3],(d3,d4],(d4,d5],(d5,d6 = max(T)]; 
####                          with dj = (j - 1)/5th quantile of T, vector d = (d1, d2, ..., d6), 
####                          interval length corresponding to λj is d[j+1] - d[j] 
 
alpha.lam0.a <- 1        # prior shape parameter for λj in piecewise constant baseline
beta.lam0.a <- 0.01      # prior rate parameter for λj in piecewise constant baseline

alpha.lam0.b <- 1        # prior shape parameter for λj in piecewise constant baseline
beta.lam0.b <- 0.01      # prior rate parameter for λj in piecewise constant baseline

alpha.lam0.c <- 1        # prior shape parameter for λj in piecewise constant baseline
beta.lam0.c <- 0.01      # prior rate parameter for λj in piecewise constant baseline


##### some useful functions
##### calculate log-likelihood of the full conditional dist. of M|X, V, theta, & observed data
Loglike <- function(B, psi, lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, dd, uu, omage, alpha1t, alpha2t, osigma, Delta, PI = matrix(0, ncol = q, nrow = q), psd, V, M, Z, S, X){  
  
    short <- array()
    invP  <- solve(diag(psi))
    gam.z.a <- gam.a[1:z]
    gam.s.a <- gam.a[(z + 1):(z + s)]
    gam.m.a <- gam.a[(z + s + 1):(z + s + q)]

    gam.z.b <- gam.b[1:z]
    gam.s.b <- gam.b[(z + 1):(z + s)]
    gam.m.b <- gam.b[(z + s + 1):(z + s + q)]

    gam.z.c <- gam.c[1:z]
    gam.s.c <- gam.c[(z + 1):(z + s)]
    gam.m.c <- gam.c[(z + s + 1):(z + s + q)]

    for (i in 1:n) {
      short[i] <- (V[i, ] - M[i, ] %*% t(B)) %*% invP %*% t(V[i,] - M[i, ] %*% t(B))
    }
    log.like1 <- -0.5*(p*log(2*pi) + log(abs(det(diag(psi)))) + short)
 
 ######  
   
    shorta1 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      parta1 <- u[,j]*(X[,3]+X[,4])*(log(lam0.a[j]) + Z %*% gam.z.a + S %*% gam.s.a + M %*% gam.m.a + omage)
      
      tempa1 <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempa1 <- tempa1 + as.numeric(lam0.a[g]*(d[g + 1] - d[g]))
        }
      }
      parta2 <- -(u[,j]*(lam0.a[j]*(X[,1] - d[j]) + tempa1))*exp(Z %*% gam.z.a + S %*% gam.s.a + M %*% gam.m.a + omage)
        
      shorta1[, j] <-  parta1 + parta2

    }
    log.likea1 <- apply(shorta1, 1, sum)
#sum(log.likea1)
##### 

    shortb1 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      partb1 <- u[,j]*X[,5]*(log(lam0.b[j]) + Z %*% gam.z.b + S %*% gam.s.b + M %*% gam.m.b + alpha1t*omage)
      
      tempb1 <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempb1 <- tempb1 + as.numeric((lam0.b[g])*(d[g + 1] - d[g]))
        }
      }
      partb2 <- -(u[,j]*((lam0.b[j])*(X[,1] - d[j]) + tempb1))*exp(Z %*% gam.z.b + S %*% gam.s.b + M %*% gam.m.b + alpha1t*omage)
        
      shortb1[, j] <-  partb1 + partb2

    }
    log.likeb1 <- apply(shortb1, 1, sum)
#sum(log.likeb1)
##### 

    shortc1 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      partc1 <- uu[,j]*X[,4]*(log(lam0.c[j]) + Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alpha2t*omage)
      
      tempc1 <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempc1 <- tempc1 + as.numeric((lam0.c[g])*(dd[g + 1] - dd[g]))
        }
      }
      partc2 <- -(X[,3]+X[,4])*(uu[,j]*((lam0.c[j])*(X[,2] - dd[j]) + tempc1))*exp(Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alpha2t*omage)
        
      shortc1[, j] <-  partc1 + partc2

    }
    log.likec1 <- apply(shortc1, 1, sum)
    
    log.like2 <- log.likea1 + log.likeb1 + log.likec1 
 
#sum(log.likec1)
# 4.91123e+175
 ###### 
 
    const <- rep(1,n)
    Dat1 <- cbind(const, Z, S)
    Coe1 <- cbind(Delta)
    Mcen <- Dat1 %*% t(Coe1) 
    invPI0 <- solve(diag(1, q) - PI)
    short2 <- array()
    SIG <- invPI0 %*% diag(psd,q) %*% t(invPI0)
    ISIG <- solve(SIG)
    for (i in 1:n) {
      short2[i] <- t(M[i, ] - invPI0 %*% Mcen[i, ]) %*% ISIG %*% (M[i, ] - invPI0 %*% Mcen[i, ])
    }
    log.like3 <- -0.5*(q*log(2*pi) + log(abs(det(SIG))) + short2)
    
    log.like <- log.like1 + log.like2 + log.like3 
    return(log.like)
}
#Loglike(B, psi, lam0.t, gam.t, lam0.c, gam.c, omage, alphat, d, u, Delta,PI, psd, V, M, Z, S, X, W)

Loglike_GA.a <- function(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, omage, alpha1t, alpha2t, osigma, M, Z, S, X){

  gam.z.a <- gam.a[1:z]
  gam.s.a <- gam.a[(z + 1):(z + s)]
  gam.m.a <- gam.a[(z + s + 1):(z + s + q)]

 ######  
   sum1 <- 0
   # shorta2 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      parta3 <- u[,j]*(X[,3]+X[,4])*(log(lam0.a[j]) + Z %*% gam.z.a + S %*% gam.s.a + M %*% gam.m.a + omage)
      
      tempa2 <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempa2 <- tempa2 + as.numeric(lam0.a[g]*(d[g + 1] - d[g]))
        }
      }
      parta4 <- -(u[,j]*(lam0.a[j]*(X[,1] - d[j]) + tempa2))*exp(Z %*% gam.z.a + S %*% gam.s.a + M %*% gam.m.a + omage)
        
     # shorta2[, j] <- parta3 + parta4
    sum1 <- sum1 + sum(parta3) + sum(parta4)  

    }
    log.like1 <- sum1
#apply(shorta2, 1, sum)
#sum(log.like1)
#######
 
  invP.a  <- diag(sig.gamma.a, n.gam.a)
  gam.a <- matrix(gam.a, ncol = 1)
  short2.a <- t(gam.a - gam0.a) %*% invP.a %*% (gam.a - gam0.a)
  log.like2.a <- -0.5*((n.gam.a)*log(2*pi) + log(abs(det(chol2inv(chol(invP.a))))) + short2.a)

  log.like <- log.like1 + log.like2.a
  return(log.like)
}

# Loglike_GA.a(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, M, Z, S, X, W)

Loglike_GA.b <- function(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, omage, alpha1t, alpha2t, osigma, M, Z, S, X){

  gam.z.b <- gam.b[1:z]
  gam.s.b <- gam.b[(z + 1):(z + s)]
  gam.m.b <- gam.b[(z + s + 1):(z + s + q)]

 ###### 
   sum1 <- 0
    #shortb2 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      partb3 <- u[,j]*X[,5]*(log(lam0.b[j]) + Z %*% gam.z.b + S %*% gam.s.b + M %*% gam.m.b + alpha1t*omage)
      
      tempb2 <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempb2 <- tempb2 + as.numeric((lam0.b[g])*(d[g + 1] - d[g]))
        }
      }
      partb4 <- -(u[,j]*((lam0.b[j])*(X[,1] - d[j]) + tempb2))*exp(Z %*% gam.z.b + S %*% gam.s.b + M %*% gam.m.b + alpha1t*omage)
        
      #shortb2[, j] <-  partb3 + partb4
    sum1 <- sum1 + sum(partb3) + sum(partb4)  
    }
    log.like1 <- sum1
   #apply(shortb2, 1, sum)
#sum(log.like1)
 ###### 

  invP.b  <- diag(sig.gamma.b, n.gam.b)
  gam.b <- matrix(gam.b, ncol = 1)
  short2.b <- t(gam.b - gam0.b) %*% invP.b %*% (gam.b - gam0.b)
  log.like2.b <- -0.5*((n.gam.b)*log(2*pi) + log(abs(det(chol2inv(chol(invP.b))))) + short2.b)

  log.like <- log.like1 + log.like2.b
  return(log.like)
}


Loglike_GA.c <- function(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, dd, uu, omage, alpha1t, alpha2t, osigma, M, Z, S, X){

  gam.z.c <- gam.c[1:z]
  gam.s.c <- gam.c[(z + 1):(z + s)]
  gam.m.c <- gam.c[(z + s + 1):(z + s + q)]


 ###### 
   sum1 <- 0
   # shortc2 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      partc3 <- uu[,j]*X[,4]*(log(lam0.c[j]) + Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alpha2t*omage)
      
      tempc3 <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempc3 <- tempc3 + as.numeric((lam0.c[g])*(dd[g + 1] - dd[g]))
        }
      }
      partc4 <- -(X[,3]+X[,4])*(uu[,j]*((lam0.c[j])*(X[,2] - dd[j]) + tempc3))*exp(Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alpha2t*omage)
       # 
      #shortc2[, j] <-  partc3 + partc4
    sum1 <- sum1 + sum(partc3) + sum(partc4)  
    }
    log.like1 <- sum1 
   #apply(shortc2, 1, sum)
 #sum(log.like1)   
 
  invP.c  <- diag(sig.gamma.c, n.gam.c)
  gam.c <- matrix(gam.c, ncol = 1)
  short2.c <- t(gam.c - gam0.c) %*% invP.c %*% (gam.c - gam0.c)
  log.like2.c <- -0.5*((n.gam.c)*log(2*pi) + log(abs(det(chol2inv(chol(invP.c))))) + short2.c)

  log.like <- log.like1 + log.like2.c
  return(log.like)
}
 
#Loglike_GA.c(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, omage, alpha1t, alpha2t, osigma, M, Z, S, X)

 

Loglike_alpha1t <- function(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, omage, alpha1t, alpha2t, osigma, M, Z, S, W){

 ###### 
   sum1 <- 0
    #shortb2 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      partb3 <- u[,j]*X[,5]*(log(lam0.b[j]) + Z %*% gam.z.b + S %*% gam.s.b + M %*% gam.m.b + alpha1t*omage)
      
      tempb2 <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempb2 <- tempb2 + as.numeric((lam0.b[g])*(d[g + 1] - d[g]))
        }
      }
      partb4 <- -(u[,j]*((lam0.b[j])*(X[,1] - d[j]) + tempb2))*exp(Z %*% gam.z.b + S %*% gam.s.b + M %*% gam.m.b + alpha1t*omage)
        
      #shortb2[, j] <-  partb3 + partb4
    sum1 <- sum1 + sum(partb3) + sum(partb4)  
    }
    log.like1 <- sum1

  log.like2.c <- dnorm(alpha1t, 0, 0.5, log=TRUE)

  log.like <- log.like1 + log.like2.c
  return(log.like)
}
 
#Loglike_alpha1t(lam0.t, gam.t, lam0.c, gam.c, omage, alpha1t, alpha2t, osigma, M, Z, S, X, W)

Loglike_alpha2t <- function(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, dd, uu, omage, alpha1t, alpha2t, osigma, M, Z, S){

  ###### 
   sum1 <- 0
   # shortc2 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      partc3 <- uu[,j]*X[,4]*(log(lam0.c[j]) + Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alpha2t*omage)
      
      tempc3 <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempc3 <- tempc3 + as.numeric((lam0.c[g])*(dd[g + 1] - dd[g]))
        }
      }
      partc4 <- -(X[,3]+X[,4])*(uu[,j]*((lam0.c[j])*(X[,2] - dd[j]) + tempc3))*exp(Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alpha2t*omage)
       # 
      #shortc2[, j] <-  partc3 + partc4
    sum1 <- sum1 + sum(partc3) + sum(partc4)  
    }
    log.like1 <- sum1 

  log.like2.c <- dnorm(alpha2t, 0, 0.5, log=TRUE)

  log.like <- log.like1 + log.like2.c
  return(log.like)
}
 

Loglike.omage <- function(B, psi, lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, dd, uu, omage, alpha1t, alpha2t, osigma, Delta, PI = matrix(0, ncol = q, nrow = q), psd, V, M, Z, S){
  
    gam.z.a <- gam.a[1:z]
    gam.s.a <- gam.a[(z + 1):(z + s)]
    gam.m.a <- gam.a[(z + s + 1):(z + s + q)]

    gam.z.b <- gam.b[1:z]
    gam.s.b <- gam.b[(z + 1):(z + s)]
    gam.m.b <- gam.b[(z + s + 1):(z + s + q)]

    gam.z.c <- gam.c[1:z]
    gam.s.c <- gam.c[(z + 1):(z + s)]
    gam.m.c <- gam.c[(z + s + 1):(z + s + q)]

   
    shorta1 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      parta1 <- u[,j]*(X[,3]+X[,4])*(log(lam0.a[j]) + Z %*% gam.z.a + S %*% gam.s.a + M %*% gam.m.a + omage)
      
      tempa1 <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempa1 <- tempa1 + as.numeric(lam0.a[g]*(d[g + 1] - d[g]))
        }
      }
      parta2 <- -(u[,j]*(lam0.a[j]*(X[,1] - d[j]) + tempa1))*exp(Z %*% gam.z.a + S %*% gam.s.a + M %*% gam.m.a + omage)
        
      shorta1[, j] <-  parta1 + parta2

    }
    log.likea1 <- apply(shorta1, 1, sum)
#sum(log.likea1)
##### 

    shortb1 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      partb1 <- u[,j]*X[,5]*(log(lam0.b[j]) + Z %*% gam.z.b + S %*% gam.s.b + M %*% gam.m.b + alpha1t*omage)
      
      tempb1 <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempb1 <- tempb1 + as.numeric((lam0.b[g])*(d[g + 1] - d[g]))
        }
      }
      partb2 <- -(u[,j]*((lam0.b[j])*(X[,1] - d[j]) + tempb1))*exp(Z %*% gam.z.b + S %*% gam.s.b + M %*% gam.m.b + alpha1t*omage)
        
      shortb1[, j] <-  partb1 + partb2

    }
    log.likeb1 <- apply(shortb1, 1, sum)
#sum(log.likeb1)
##### 

    shortc1 <- array(0, dim = c(n, G))
    for (j in 1:G) {
      partc1 <- uu[,j]*X[,4]*(log(lam0.c[j]) + Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alpha2t*omage)
      
      tempc1 <- 0
      if (j > 1) {
        for (g in 1:(j - 1)) {
          tempc1 <- tempc1 + as.numeric((lam0.c[g])*(dd[g + 1] - dd[g]))
        }
      }
      partc2 <- -(X[,3]+X[,4])*(uu[,j]*((lam0.c[j])*(X[,2] - dd[j]) + tempc1))*exp(Z %*% gam.z.c + S %*% gam.s.c + M %*% gam.m.c + alpha2t*omage)
        
      shortc1[, j] <-  partc1 + partc2

    }
    log.likec1 <- apply(shortc1, 1, sum)
    
    log.like2 <- log.likea1 + log.likeb1 + log.likec1 
 
###########    
    log.like4 <- sum(-(t(omage)%*%omage)/(2*osigma))        
    log.like <- log.like2 + log.like4 
    return(log.like)
} 
  

####### Store the mean and sd. of each parameter in all replications
mean.B <- sd.B <- array(0, dim = c(rep_num, n.B))
mean.psi <- sd.psi <- array(0, dim = c(rep_num, n.psi))
mean.L.se <- sd.L.se <- array(0, dim = c(rep_num, n.se))
mean.psd <- sd.psd <- array(0, dim = c(rep_num, n.psd))
mean.gam.a <- sd.gam.a <- array(0, dim = c(rep_num, n.gam.a))
mean.lam0.a <- sd.lam0.a <- array(0, dim = c(rep_num, G))
mean.gam.b <- sd.gam.b <- array(0, dim = c(rep_num, n.gam.b))
mean.lam0.b <- sd.lam0.b <- array(0, dim = c(rep_num, G))
mean.gam.c <- sd.gam.c <- array(0, dim = c(rep_num, n.gam.c))
mean.lam0.c <- sd.lam0.c <- array(0, dim = c(rep_num, G))
mean.alpha1t <- sd.alpha1t <- array(0, dim = c(rep_num, 1))
mean.alpha2t <- sd.alpha2t <- array(0, dim = c(rep_num, 1))
mean.omage <- sd.omage <- array(0, dim = c(rep_num, n))

####### Store the quantile of each parameter in all replications
q.B <- array(0, dim = c(rep_num, n.B*2))
q.psi <- array(0, dim = c(rep_num, n.psi*2))
q.L.se <- array(0, dim = c(rep_num, n.se*2))
q.psd <- array(0, dim = c(rep_num, n.psd*2))
q.gam.a <- array(0, dim = c(rep_num, n.gam.a*2))
q.lam0.a <- array(0, dim = c(rep_num, G*2))
q.gam.b <- array(0, dim = c(rep_num, n.gam.b*2))
q.lam0.b <- array(0, dim = c(rep_num, G*2))
q.gam.c <- array(0, dim = c(rep_num, n.gam.c*2))
q.lam0.c <- array(0, dim = c(rep_num, G*2))
q.alpha1t <- array(0, dim = c(rep_num, 1*2))
q.alpha2t <- array(0, dim = c(rep_num, 1*2))


##############################################################################
 
##### for storing results
result <- list('B' = array(NA, dim = c(n.mcmc, n.B)),
               'alpha1t' = array(NA, dim = c(n.mcmc, 1)),
               'alpha2t' = array(NA, dim = c(n.mcmc, 1)),
               'omage' = array(NA, dim = c(n.mcmc, n)),
               'L.se' = array(NA, dim = c(n.mcmc, n.se)),
               'gam.a' = array(NA, dim = c(n.mcmc, n.gam.a)),
               'gam.b' = array(NA, dim = c(n.mcmc, n.gam.b)),
               'gam.c' = array(NA, dim = c(n.mcmc, n.gam.c)),
               'psd' = array(NA, dim = c(n.mcmc, q)),
               'psi' = array(NA, dim = c(n.mcmc, n.psi)), 
               'lam0.a' = array(NA, dim = c(n.mcmc, G)),
               'lam0.b' = array(NA, dim = c(n.mcmc, G)),
               'lam0.c' = array(NA, dim = c(n.mcmc, G)))

t0 <- Sys.time()
cat('Generating', rep_num, 'sets of data...  Please wait \n')
#############
#simulate data with censoring

 
for (crep in 1:rep_num) {
  set.seed(200 + crep *10)
  ################################
  ######  data generation   ######
  S <- t(rmultinom(n, 1, prob = c(0.5, 0.2, 0.3)))
  S <- S[, 2:3]
  Z <- rmvnorm(n, mean = rep(0,z), sigma = diag(1,z))
 
  X <- array(NA, dim = c(n, 6))
 
  Dat <- cbind(rep(1,n), Z, S)
  Coe <- cbind(Delta.true)
  M <- Dat %*% t(Coe) + rmvnorm(n, mean = rep(0, q), sigma = diag(psd.true, q))
  M.true <- M <- M %*% t(solve((diag(1, q) - PI)))
  const <- rep(1,n)
  Dat1 <- cbind(M)
  V <- Dat1 %*% t(L.me.true) + rmvnorm(n, mean = rep(0, p), sigma = diag(psi.true))
  
  # censoring time and survival time
  Dat2 <- cbind(Z, S, M)
  temp.a <- exp(Dat2 %*% gam.true.a + omegaT)
  temp.b <- exp(Dat2 %*% gam.true.b + alpha1T*omegaT)
  temp.c <- exp(Dat2 %*% gam.true.c + alpha2T*omegaT)

  # with lambda_0(t) = 1
  TTP <- (-1/(1*temp.a))*log(runif(n))
  OS_orig <- (-1/(1*temp.b))*log(runif(n))
  OSP <- (-1/(1*temp.c))*log(runif(n))
  # with lamda_0(t) = 2t + 1
  #TTP <- sqrt((-1/temp.a)*log(runif(n)) + 1/4) - 1/2
  #OS_orig <- sqrt((-1/temp.b)*log(runif(n)) + 1/4) - 1/2
  #OSP <- sqrt((-1/temp.c)*log(runif(n)) + 1/4) - 1/2
  
  
  C <- rep(4, n); 
#runif(n, min = 0, max = c[2])
  PFS = pmin(TTP, OS_orig);

  del=rep(0,n);
  OS=rep(0,n);
  delta1=rep(0,n);
  delta2=rep(0,n);
  delta3=rep(0,n);
  delta4=rep(0,n);
  ti1=rep(0,n);
  ti2=rep(0,n);

  for(i in 1:n){
    if (PFS[i] == TTP[i]) {
      ti1[i]=PFS[i]
      ti2[i]=pmin((C-TTP)[i],OSP[i])
      if(ti2[i]==(C-TTP)[i]){
        delta1[i]=1
        del[i]=0
		OS[i]=C[i]   
      }else {
        delta2[i]=1
	  del[i]=1
		OS[i]=(TTP+OSP)[i]
      }
    }
   else {
      ti1[i]=pmin(C[i],OS_orig[i])
      if(ti1[i]==OS_orig[i]){
        delta3[i]=1
        del[i]=1
		OS[i]=OS_orig[i]
      }else{
        delta4[i]=1
	  del[i]=0
		OS[i]=C[i]
      }
    }
  }
  data=data.frame(ti1,ti2,delta1,delta2,delta3,delta4);  

  X[, 1] <- ti1
  X[, 2] <- ti2
  X[, 3] <- delta1
  X[, 4] <- delta2
  X[, 5] <- delta3
  X[, 6] <- delta4
 
length(which(X[,3]==1))/n
length(which(X[,4]==1))/n
length(which(X[,5]==1))/n
length(which(X[,6]==1))/n

#delta1+delta2+delta3+delta4
print(c(length(which(delta1==1)), length(which(delta2==1)), length(which(delta3==1)), length(which(delta4==1))))
print(c(length(which(ti1<0)), length(which(ti2<0))))

### generate initial values for latent mediator M1
###init value 1 for parameters

B <- matrix(c(
  1,0,
  0,0,
  0,1,
  0,0
), nrow = p, ncol = q, byrow = TRUE)

L.me <- B

psi <- rep(1, p)  

Delta <- matrix(rep(0, q*(1 + s + z)), nrow = q, ncol = (1 + z + s), byrow = TRUE)

L.se <- cbind(Delta)

psd <- rep(1, q)  

gam.z.a <- c(0, 0, 0, 0)
gam.s.a <- c(0, 0)
gam.m.a <- c(0, 0)
gam.a <- c(gam.z.a, gam.s.a, gam.m.a)
lam0.a <- rep(1, G)

gam.z.b <- c(0, 0, 0, 0)
gam.s.b <- c(0, 0)
gam.m.b <- c(0, 0)
gam.b <- c(gam.z.b, gam.s.b, gam.m.b)
lam0.b <- rep(1, G)

gam.z.c <- c(0, 0, 0, 0)
gam.s.c <- c(0, 0)
gam.m.c <- c(0, 0)
gam.c <- c(gam.z.c, gam.s.c, gam.m.c)
lam0.c <- rep(1, G)

alpha1t = 0.1
alpha2t = 0.1
osigma = 0.1
omage = omegaT

# to facilitate computation later 
iv.psi <- 1/psi
iv.sqrt.psi <- sqrt(iv.psi)

iv.psd <- 1/psd
iv.sqrt.psd <- sqrt(iv.psd)

if (q > 0) {
  const <- rep(1,n)
  Dat <- cbind(const, Z, S)
  Coe <- cbind(Delta)
  M <- Dat %*% t(Coe) + rmvnorm(n, mean = rep(0, q), sigma = diag(psd, q))
  M <- M %*% t(solve((diag(1, q) - PI)))   
}
 
   GC1  = (X[ , 1])[!duplicated((X[ , 1]))]

##将所有样本量N中不重复的按升序排序##
   GG1  = GC1[order(GC1,decreasing=F)]  
  
  d <- c(quantile(GG1, probs = seq(0, 1, 1/G)))
  d[1] <- 0
  u <- array(0,dim = c(n, G))
  for (j in 1:G) {
    al <- as.logical(d[j] < X[,1])*(X[,1] <= d[j + 1])
    u[, j] <- al
  }

 
  GC2  = (X[ , 2])[!duplicated((X[ , 2]))]
##将所有样本量N中不重复的按升序排序##
  GG2  = GC2[order(GC2,decreasing=F)] 

  dd <- c(quantile(GG2, probs = seq(0, 1, 1/G)))
    dd[1] <- 0

  uu <- array(0,dim = c(n, G))
  for (j in 1:G) {
    a2 <- as.logical(dd[j] < X[ , 2])*(X[ , 2] <= dd[j + 1])
    uu[, j] <- a2
  }


n.accept.gam.a = 1
n.accept.gam.b = 1
n.accept.gam.c = 1
n.accept.M = rep(1,n) 
n.accept.lam0.c= rep(1,G)
n.accept.alpha1t = 1
n.accept.alpha2t = 1
n.accept.omage = rep(1,n) 

   ## iteration
   it=1
   while (it<n.mcmc+1)
{
cat("loopp", crep, "Iteration", it, fill=TRUE);  
 

  # step1: update latent mediators M
  if (q > 0) {
    PI0 <- (diag(1, q) - PI)
    ISG <- crossprod(iv.sqrt.psi * B) + crossprod(iv.sqrt.psd * PI0)     # ∑_ω^-1
    SIG <- chol2inv(chol(ISG))                     
    SIG <- sig.mh*SIG
    cSIG <- chol(SIG)
  }
  #### Random Walk Metropolis, M_t + N[0, σ_mh^2*∑_ω_t]S %*% gam.s.a
  M.new <- M + t(crossprod(cSIG, matrix(rnorm(q*n), nrow = q)))
  ll1 <- Loglike(B, psi, lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, dd, uu, omage, alpha1t, alpha2t, osigma, Delta, PI, psd, V, M, Z, S, X)
  ll2 <- Loglike(B, psi, lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, dd, uu, omage, alpha1t, alpha2t, osigma, Delta, PI, psd, V, M.new, Z, S, X)
  ## acceptance ratio
  p.accept <- exp(ll2 - ll1)
  accept <- (runif(n) < p.accept)  
  M[accept, ] <- M.new[accept, ]
  n.accept.M <- n.accept.M + accept
  #n.accept.M/n.mcmc
  # step2: update loading matrix B and residual variance psi
  count.B <- 1
  L.me <- B
  omg.me <- M
  for (k in 1:p) {
    free <- Id.me[k, ]                  # which column on row K is free/fixed
    len <- n.me.row[k]
    Vcen <- t(L.me[k, !free, drop = F] %*% t(omg.me)[!free, , drop = F])
    #Mk <- t(omg.me)               # data of the latent variable that the kth indicator loads on
    Psiginv <- rep(sig.b, len)          # H_0yk
 
    Vk.star <- V[ , k] - Vcen
    alpha.psi.star <- alpha.psi + 0.5*n
    beta.psi.star <- beta.psi + 0.5*sum(Vk.star^2)
    if (len > 0) {
	if(len==1){Mk=matrix(omg.me[,free],nrow=1)}
      if(len>1){Mk<-omg.me[,free]}
      A_vk <- chol2inv(chol(diag(Psiginv, len) + tcrossprod(Mk)))
      temp <- Psiginv*L.me0[k, free] + Mk %*% Vk.star
      a_vk <- A_vk %*% temp
      beta.psi.star <- beta.psi.star + 0.5*(sum(L.me0[k, free]*Psiginv*L.me0[k, free]) - sum(temp*a_vk))
    }

    iv.psi[k] <- rgamma(1, shape = alpha.psi.star, rate = beta.psi.star)
    psi[k] <- 1/iv.psi[k]                  
    iv.sqrt.psi[k] <- sqrt(iv.psi[k])

    if (len > 0) {
      L.me[k, free] <- rmvnorm(1, a_vk, sigma = psi[k]*A_vk)
      if (n.b.row[k] > 0) {
        B[k, ] <- L.me[k, ]
      }
 
        if (n.b.row[k] > 0) {
          result$B[it, count.B:(count.B + n.b.row[k] - 1)] <- B[k, Id.B[k, ]]
        }
 
      count.B <- count.B + n.b.row[k]
    }   
  }

  # step3: update regression coefficients delta and mean_m μ and error variance psd
  ## Note!! Need revise if the number of latent mediator changed!
  ## Note!! Need revise if the number of latent mediator changed!
  count.se <- 1
  L.se <- cbind(Delta)
  omg.se <- cbind(rep(1,n), Z, S)
  for (k in 1:q) {
    free <- Id.se[k, ]  # which column on row K is free/fixed
    len <- n.se.row[k]
    Mcen <- t(L.se[k, !free, drop = F] %*% t(omg.se)[!free, , drop = F])
    Mk.star <- M[ , k] - as.vector(Mcen)
    alpha.psd.star <- alpha.psd + 0.5*n
    beta.psd.star <- beta.psd + 0.5*sum(Mk.star^2)

    if (len > 0) {
      Yk <- omg.se[, free, drop = F]  
      iH0dk <- diag(sig.delta, len)
      L.se0k <- L.se0[k, free]
      A_dk <- chol2inv(chol(iH0dk + crossprod(Yk)))
      temp <- iH0dk %*% L.se0k + t(Yk) %*% Mk.star
      a_dk <- A_dk %*% temp
      beta.psd.star <- beta.psd.star + 0.5*(crossprod(crossprod(iH0dk, L.se0k), L.se0k)
                                    - crossprod(a_dk, temp))
    }

    iv.psd[k] <- rgamma(1, shape = alpha.psd.star, rate = beta.psd.star)
    psd[k] <- 1/iv.psd[k]
    iv.sqrt.psd[k] <- sqrt(iv.psd[k])

    if (len > 0) {
      L.se[k, free] <- rmvnorm(1, a_dk, sigma = psd[k]*A_dk)

        result$L.se[it, count.se:(count.se + len - 1)] <- L.se[k, free]

      count.se <- count.se + len
    }
  }

  Delta <- L.se[, 1:(1 + z + s), drop = F]
 
    result$psd[it, ] <- psd
    result$psi[it, ] <- psi
  

  # step4: update regression coefficients gam, M-H I guess
  #### Random Walk Metropolis, gam_t + N[0, σ_mh^2*I}

 gam.new.a <- gam.a + rmvnorm(1, rep(0, n.gam.a), sigma = diag(rep(0.0015, n.gam.a)))
 llga1 <- Loglike_GA.a(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, omage, alpha1t, alpha2t, osigma, M, Z, S, X)
 llga2 <- Loglike_GA.a(lam0.a, gam.new.a, lam0.b, gam.b, lam0.c, gam.c, d, u, omage, alpha1t, alpha2t, osigma, M, Z, S, X)
 ## acceptance ratio
 p.accept <- exp(llga2 - llga1)
 accept <- (runif(1) < p.accept)
 if (accept) gam.a <-  gam.new.a
 gam.z.a <- gam.a[1:z]
 gam.s.a <- gam.a[(z + 1):(z + s)]
 gam.m.a <- gam.a[(z + s + 1):(z + s + q)]
 n.accept.gam.a <- n.accept.gam.a + accept
 #n.accept.gam.a/n.mcmc
    result$gam.a[it, ] <- gam.a

  # step4: update regression coefficients gam, M-H I guess
  #### Random Walk Metropolis, gam_t + N[0, σ_mh^2*I}

 gam.new.b <- gam.b + rmvnorm(1, rep(0, n.gam.b), sigma = diag(rep(0.0015, n.gam.b)))
 llgb1 <- Loglike_GA.b(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, omage, alpha1t, alpha2t, osigma, M, Z, S, X)
 llgb2 <- Loglike_GA.b(lam0.a, gam.a, lam0.b, gam.new.b, lam0.c, gam.c, d, u, omage, alpha1t, alpha2t, osigma, M, Z, S, X)
 ## acceptance ratio
 p.accept <- exp(llgb2 - llgb1)
 accept <- (runif(1) < p.accept)
 if (accept) gam.b <-  gam.new.b
 gam.z.b <- gam.b[1:z]
 gam.s.b <- gam.b[(z + 1):(z + s)]
 gam.m.b <- gam.b[(z + s + 1):(z + s + q)]
 n.accept.gam.b <- n.accept.gam.b + accept
 #n.accept.gam.b/n.mcmc
    result$gam.b[it, ] <- gam.b
 
  # step4: update regression coefficients gam, M-H I guess
  #### Random Walk Metropolis, gam_t + N[0, σ_mh^2*I}

 gam.new.c <- gam.c + rmvnorm(1, rep(0, n.gam.c), sigma = diag(rep(0.0015, n.gam.c)))
 llgc1 <- Loglike_GA.c(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, dd, uu, omage, alpha1t, alpha2t, osigma, M, Z, S, X)
 llgc2 <- Loglike_GA.c(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.new.c, dd, uu, omage, alpha1t, alpha2t, osigma, M, Z, S, X)
 ## acceptance ratio
 p.accept <- exp(llgc2 - llgc1)
 accept <- (runif(1) < p.accept)
 if (accept) gam.c <-  gam.new.c
 gam.z.c <- gam.c[1:z]
 gam.s.c <- gam.c[(z + 1):(z + s)]
 gam.m.c <- gam.c[(z + s + 1):(z + s + q)]
 n.accept.gam.c <- n.accept.gam.c + accept
 #n.accept.gam.c/n.mcmc
    result$gam.c[it, ] <- gam.c 
 
  # step5: update piecewise constant hazard lam0, conjugate gamma prior?

  alpha.lam0.star.a <- rep(alpha.lam0.a, G) + as.vector(t(X[,3]+X[,4]) %*% u)

  temp2.a <- array(0, dim = c(n, G))

  for (j in 1:G) {
    temp3.a <- 0
        if (j < G) {
          for (g in (j + 1):G) {
            temp3.a <- temp3.a + u[, g]*(d[j + 1] - d[j])
          }
        }
    temp2.a[,j] <- exp(cbind(Z, S, M) %*% as.vector(gam.a) + omage)*(u[,j]*(X[,1] - d[j]) + temp3.a)
  }

  beta.lam0.star.a <- rep(beta.lam0.a, G) + colSums(temp2.a)

  for (j in 1:G) {
    lam0.a[j] <- rgamma(1, shape = alpha.lam0.star.a[j], rate = beta.lam0.star.a[j])
  }
 
     result$lam0.a[it, ] <- lam0.a
  
 
  # step5: update piecewise constant hazard lam0, conjugate gamma prior?
   
  alpha.lam0.star.b <- rep(alpha.lam0.b, G) + as.vector(t(X[,5]) %*% u)

  temp2.b <- array(0, dim = c(n, G))

  for (j in 1:G) {
    temp3.b <- 0
        if (j < G) {
          for (g in (j + 1):G) {
            temp3.b <- temp3.b + u[, g]*(d[j + 1] - d[j])
          }
        }
    temp2.b[,j] <- exp(cbind(Z, S, M) %*% as.vector(gam.b) + alpha1t*omage)*(u[,j]*(X[,1] - d[j]) + temp3.b)
  }

  beta.lam0.star.b <- rep(beta.lam0.b, G) + colSums(temp2.b)

  for (j in 1:G) {
    lam0.b[j] <- rgamma(1, shape = alpha.lam0.star.b[j], rate = beta.lam0.star.b[j])
  }
 
     result$lam0.b[it, ] <- lam0.b



  # step5: update piecewise constant hazard lam0, conjugate gamma prior?
 
################################################################## 
 
  alpha.lam0.star.c <- rep(alpha.lam0.c, G) + as.vector(t(X[,4]) %*% uu)

  temp2.c <- array(0, dim = c(n, G))

  for (j in 1:G) {
    temp3.c <- 0
        if (j < G) {
          for (g in (j + 1):G) {
            temp3.c <- temp3.c + uu[, g]*(dd[j + 1] - dd[j])
          }
        }
    temp2.c[,j] <- exp(cbind(Z, S, M) %*% as.vector(gam.c) + alpha2t*omage)*(uu[,j]*(X[,2] - dd[j]) + temp3.c)*(X[,3] + X[,4])
  }

  beta.lam0.star.c <- rep(beta.lam0.c, G) + colSums(temp2.c)

  for (j in 1:G) {
    lam0.c[j] <- rgamma(1, shape = alpha.lam0.star.c[j], rate = beta.lam0.star.c[j])
  }
 
     result$lam0.c[it, ] <- lam0.c 
 

 
 #lam0.new.c <- as.vector(lam0.c + rmvnorm(1, rep(0, G), sigma = diag(rep(0.0015, G))))
 #llla1 <- Loglike_LA.c(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, dd, uu, omage, alpha1t, alpha2t, osigma, M, Z, S, X, W)
 #llla2 <- Loglike_LA.c(lam0.a, gam.a, lam0.b, gam.b, lam0.new.c, gam.c, dd, uu, omage, alpha1t, alpha2t, osigma, M, Z, S, X, W)
 #acceptance ratio
 #la.accept <- exp(llla2 - llla1)
 #accept <- (runif(1) < la.accept)
 #if (accept) lam0.c <- lam0.new.c
 #n.accept.lam0.c <- n.accept.lam0.c + accept
 #n.accept.lam0.c/n.mcmc
 #result$lam0.c[it, ] <- lam0.c

################################################################## 
 

 alpha1t.new <- alpha1t + rnorm(1, 0, 0.25)
 lla1 <- Loglike_alpha1t(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, omage, alpha1t, alpha2t, osigma, M, Z, S)
 lla2 <- Loglike_alpha1t(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, omage, alpha1t.new, alpha2t, osigma, M, Z, S)
 ## acceptance ratio
 a.accept <- exp(lla2 - lla1)
 accept <- (runif(1) < a.accept)
 if (accept) alpha1t <-  alpha1t.new
  alpha1t[accept] <- alpha1t.new[accept]  
 n.accept.alpha1t <- n.accept.alpha1t + accept
  result$alpha1t[it] <- alpha1t
 #n.accept.alpha1t/n.mcmc 
 #alpha1T 


 alpha2t.new <- alpha2t + rnorm(1, 0, 0.25)
 llaa1 <- Loglike_alpha2t(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, dd, uu, omage, alpha1t, alpha2t, osigma, M, Z, S)
 llaa2 <- Loglike_alpha2t(lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, dd, uu, omage, alpha1t, alpha2t.new, osigma, M, Z, S)
 ## acceptance ratio
 a.accept <- exp(llaa2 - llaa1)
 accept <- (runif(1) < a.accept)
 if (accept) alpha2t <-  alpha2t.new
  alpha2t[accept] <- alpha2t.new[accept]  
 n.accept.alpha2t <- n.accept.alpha2t + accept
  result$alpha2t[it] <- alpha2t
 #n.accept.alpha2t/n.mcmc 
 #alpha2T 

# step8
 

  #### Random Walk Metropolis, M_t + N[0, σ_mh^2*∑_ω_t]
  omage.new <- omage + rnorm(n, 0, 0.035)
  lo1 <- Loglike.omage(B, psi, lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, dd, uu, omage, alpha1t, alpha2t, osigma, Delta, PI, psd, V, M, Z, S)
  lo2 <- Loglike.omage(B, psi, lam0.a, gam.a, lam0.b, gam.b, lam0.c, gam.c, d, u, dd, uu, omage.new, alpha1t, alpha2t, osigma, Delta, PI, psd, V, M, Z, S)
  ## acceptance ratio
  o.accept <- exp(lo2 - lo1)
  accept <- (runif(n) < o.accept)
  omage[accept] <- omage.new[accept]  
  n.accept.omage <- n.accept.omage + accept
result$omage[it, ] <- omage
  #n.accept.omage/n.mcmc 
 #omegaT -  omage
#
library(invgamma)
#actovr
aa=0.05
bb=0.05
osigma <- rinvgamma(1, aa+n/2, bb+0.5*t(omage)%*%(omage))
#omage-omegaT

  it=it+1
} # end of mcmc
 
  mean.B[crep, ] <- colMeans(result$B[n.burn:n.mcmc,])
  mean.psi[crep, ] <- colMeans(result$psi[n.burn:n.mcmc,])
  mean.L.se[crep, ] <- colMeans(result$L.se[n.burn:n.mcmc,])
  mean.psd[crep, ] <- colMeans(result$psd[n.burn:n.mcmc,])
  mean.gam.a[crep, ] <- colMeans(result$gam.a[n.burn:n.mcmc,])
  mean.lam0.a[crep, ] <- colMeans(result$lam0.a[n.burn:n.mcmc,])
  mean.gam.b[crep, ] <- colMeans(result$gam.b[n.burn:n.mcmc,])
  mean.lam0.b[crep, ] <- colMeans(result$lam0.b[n.burn:n.mcmc,])
  mean.gam.c[crep, ] <- colMeans(result$gam.c[n.burn:n.mcmc,])
  mean.lam0.c[crep, ] <- colMeans(result$lam0.c[n.burn:n.mcmc,])
  mean.omage[crep,] <- colMeans(result$omage[n.burn:n.mcmc,])
#var(apply(result$omage[,,n.burn:n.mcmc],c(1,2),mean))
  mean.alpha1t[crep] <- mean(result$alpha1t[n.burn:n.mcmc])
  mean.alpha2t[crep] <- mean(result$alpha2t[n.burn:n.mcmc])
 

####################################################################################

  
  sd.B[crep, ] <- apply(result$B[n.burn:n.mcmc,], 2, sd)
  sd.psi[crep, ] <- apply(result$psi[n.burn:n.mcmc,],2, sd)
  sd.L.se[crep, ] <- apply(result$L.se[n.burn:n.mcmc,], 2, sd)
  sd.psd[crep, ] <- apply(result$psd[n.burn:n.mcmc,], 2, sd)
  sd.gam.a[crep, ] <- apply(result$gam.a[n.burn:n.mcmc,], 2, sd)
  sd.lam0.a[crep, ] <- apply(result$lam0.a[n.burn:n.mcmc,], 2, sd)
  sd.gam.b[crep, ] <- apply(result$gam.b[n.burn:n.mcmc,], 2, sd)
  sd.lam0.b[crep, ] <- apply(result$lam0.b[n.burn:n.mcmc,], 2, sd)
  sd.gam.c[crep, ] <- apply(result$gam.c[n.burn:n.mcmc,], 2, sd)
  sd.lam0.c[crep, ] <- apply(result$lam0.c[n.burn:n.mcmc,], 2, sd)
  sd.alpha1t[crep] <- sd(result$alpha1t[n.burn:n.mcmc])
  sd.alpha2t[crep] <- sd(result$alpha2t[n.burn:n.mcmc])
  #sd.omage[crep] <- sd(result$omage[n.burn:n.mcmc])

  
  q.B[crep, ] <- as.vector(apply(result$B[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.psi[crep, ] <- as.vector(apply(result$psi[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.L.se[crep, ] <- as.vector(apply(result$L.se[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.psd[crep,] <- as.vector(apply(result$psd[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.gam.a[crep, ] <- as.vector(apply(result$gam.a[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.lam0.a[crep, ] <- as.vector(apply(result$lam0.a[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.gam.b[crep, ] <- as.vector(apply(result$gam.b[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.lam0.b[crep, ] <- as.vector(apply(result$lam0.b[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.gam.c[crep, ] <- as.vector(apply(result$gam.c[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  q.lam0.c[crep, ] <- as.vector(apply(result$lam0.c[n.burn:n.mcmc,], 2, quantile, c(0.025, 0.975), na.rm = TRUE))
  #q.omage[crep, ] <- as.vector(quantile(result$omage[n.burn:n.mcmc], c(0.025, 0.975), na.rm = TRUE))
  q.alpha1t[crep, ] <- as.vector(quantile(result$alpha1t[n.burn:n.mcmc], c(0.025, 0.975), na.rm = TRUE))
  q.alpha2t[crep, ] <- as.vector(quantile(result$alpha2t[n.burn:n.mcmc], c(0.025, 0.975), na.rm = TRUE))

  
}
 

B.list <- list(Mean = colMeans(mean.B), Sd = colMeans(sd.B),Se = apply(mean.B,2,sd),
               Bias = colMeans(mean.B) - B.true[Id.B],
               Rms = sqrt(rowSums((t(mean.B) - B.true[Id.B])^2) / rep_num),
               p.CI95 = colMeans((q.B[ , (1:n.B)*2 - 1] - t(replicate(
                 rep_num, B.true[Id.B])) > 0) * (q.B[, (1:n.B)*2] - t(
                   replicate(rep_num, B.true[Id.B])) > 0) == 0))

psi.list <- list(Mean = colMeans(mean.psi), Sd = colMeans(sd.psi),Se = apply(mean.psi,2,sd),
                 Bias = colMeans(mean.psi) - psi.true[Id.psi],
                 Rms = sqrt(rowSums((t(mean.psi) - psi.true[Id.psi])^2) / rep_num),
                 p.CI95 = colMeans((q.psi[ , (1:n.psi)*2 - 1] - t(replicate(
                   rep_num, psi.true[Id.psi])) > 0) * (q.psi[, (1:n.psi)*2] - t(
                     replicate(rep_num, psi.true[Id.psi])) > 0) == 0))

L.se.list <- list(Mean = colMeans(mean.L.se), Sd = colMeans(sd.L.se),Se = apply(mean.L.se,2,sd),
               Bias = colMeans(mean.L.se) - t(L.se.true)[t(Id.se)],
               Rms = sqrt(rowSums((t(mean.L.se) - t(L.se.true)[t(Id.se)])^2) / rep_num),
               p.CI95 = colMeans((q.L.se[ , (1:n.se)*2 - 1] - t(replicate(
                 rep_num, t(L.se.true)[t(Id.se)])) > 0) * (q.L.se[, (1:n.se)*2] - t(
                   replicate(rep_num, t(L.se.true)[t(Id.se)])) > 0) == 0))

psd.list <- list(Mean = colMeans(mean.psd), Sd = colMeans(sd.psd),Se = apply(mean.psd,2,sd),
                 Bias = colMeans(mean.psd) - t(psd.true),
                 Rms = sqrt(rowSums((t(mean.psd) - psd.true[Id.psd])^2) / rep_num),
                 p.CI95 = colMeans((q.psd[ , (1:n.psd)*2 - 1] - t(replicate(
                   rep_num, psd.true[Id.psd])) > 0) * (q.psd[, (1:n.psd)*2] - t(
                     replicate(rep_num, psd.true[Id.psd])) > 0) == 0))

gam.a.list <- list(Mean = colMeans(mean.gam.a), Sd = colMeans(sd.gam.a),Se = apply(mean.gam.a,2,sd),
                 Bias = colMeans(mean.gam.a) - gam.true.a,
                 Rms = sqrt(rowSums((t(mean.gam.a) - gam.true.a)^2) / rep_num),
                 p.CI95 = colMeans((q.gam.a[ , (1:n.gam.a)*2 - 1] - t(replicate(
                   rep_num, gam.true.a)) > 0) * (q.gam.a[, (1:n.gam.a)*2] - t(
                     replicate(rep_num, gam.true.a)) > 0) == 0))

gam.b.list <- list(Mean = colMeans(mean.gam.b), Sd = colMeans(sd.gam.b),Se = apply(mean.gam.b,2,sd),
                 Bias = colMeans(mean.gam.b) - gam.true.b,
                 Rms = sqrt(rowSums((t(mean.gam.b) - gam.true.b)^2) / rep_num),
                 p.CI95 = colMeans((q.gam.b[ , (1:n.gam.b)*2 - 1] - t(replicate(
                   rep_num, gam.true.b)) > 0) * (q.gam.b[, (1:n.gam.b)*2] - t(
                     replicate(rep_num, gam.true.b)) > 0) == 0))

gam.c.list <- list(Mean = colMeans(mean.gam.c), Sd = colMeans(sd.gam.c),Se = apply(mean.gam.c,2,sd),
                 Bias = colMeans(mean.gam.c) - gam.true.c,
                 Rms = sqrt(rowSums((t(mean.gam.c) - gam.true.c)^2) / rep_num),
                 p.CI95 = colMeans((q.gam.c[ , (1:n.gam.c)*2 - 1] - t(replicate(
                   rep_num, gam.true.c)) > 0) * (q.gam.c[, (1:n.gam.c)*2] - t(
                     replicate(rep_num, gam.true.c)) > 0) == 0))


est.list <- list(B = B.list, psi = psi.list, L.se = L.se.list,
                 psd = psd.list, gam.a = gam.a.list, gam.b = gam.b.list,  
                 gam.c = gam.c.list)
 

for (i in 1:length(est.list)) {
  write.csv(est.list[[i]], paste(names(est.list)[i],".csv",sep = ""))
}


save.image(paste(crep, '.Rdata', sep = ''))

 Sys.time()-t0; 
