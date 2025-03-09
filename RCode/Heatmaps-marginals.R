rm(list = ls(all.names = TRUE))

###
### Maximum likelihood estimation functions
###

# Binomial probability density
pn_binom <- function(x,xR,alpha,nu,nA0,t){
  theta <- alpha*xR
  p0 <- (1+exp(-(sum(theta)+nu)*t)*sum(theta)/nu)/(1+sum(theta)/nu)
  prob <- c(p0)
  return(dbinom(x,size=nA0,prob))
}

# Log-likelihood for a single time
lk_binom_one_time <- function(x,xR,alpha,nu,nA0,t){
  return(-log(pn_binom(x,xR,alpha,nu,nA0,t)))
}

# Additive log-likelihood for a dataset evaluated at different times
lk_binom <- function(dataset,xR,alpha,nu,nA0){
  times <- dataset[1,]
  ntimes <- length(times)
  logL <- 0
  S <- length(dataset[,1])
  for (i in 1:ntimes){
    x <- dataset[2,i]
    logL <- logL + lk_one_time(x,xR,alpha,nu,nA0,times[i])
  }
  return(logL)
}

# Likelihood optimization function
lkopt_binom <- function(data, par){
  alphap <- par[1]
  nup <- par[2]
  S <- data[1][[1]]
  nA0 <- data[2][[1]]
  xR <- data[3][[1]]
  dataset <- data[4][[1]]
  times <- dataset[1,]
  ntimes <- length(times)
  logL <- 0
  for (i in 1:ntimes){
    fin <- S+2
    x <- dataset[2,i]
    logL <- logL + lk_binom_one_time(x,xR,alphap,nup,nA0,times[i])
  }
  return(logL)
}

###
### Single realization of Gillespie Stochastic Simulation Algorithm
###

library(GillespieSSA)

# Parameters
S <- 5          # Number of resources
alpha <- 5      # Attack rate
nup <- 10       # Relaxation rate
nA0 <- 100      # Total number of consumers
xR <- seq(0.9,0.1,length.out=S)   # Vector of resource densities

# Model definitions
simName <- "Multi-resource HII"   # Name
x0 <- c(nA0, rep(0,S))            # Initially, all predators are free
names(x0) <- c("A","AR1","AR2","AR3","AR4","AR5")

theta <- alpha*xR
parms <- c(t1 = theta[1], t2 = theta[2], t3 = theta[3], 
           t4 = theta[4], t5 = theta[5], nup = nup, nA0 = nA0)

# Propensity matrix
nu <- matrix(c(-1, 1, -1, 1, -1, 1, -1, 1, -1, 1,
               1, -1,  0, 0,  0, 0,  0, 0,  0, 0,
               0,  0,  1, -1, 0, 0,  0, 0,  0, 0,
               0,  0,  0, 0,  1, -1, 0, 0,  0, 0,
               0,  0,  0, 0,  0, 0,  1, -1, 0, 0,
               0,  0,  0, 0,  0, 0,  0, 0,  1, -1), 
             nrow = (S+1), byrow = TRUE)
a <- c("t1*(nA0-(AR1+AR2+AR3+AR4+AR5))", "nup*AR1",
       "t2*(nA0-(AR1+AR2+AR3+AR4+AR5))", "nup*AR2",
       "t3*(nA0-(AR1+AR2+AR3+AR4+AR5))", "nup*AR3",
       "t4*(nA0-(AR1+AR2+AR3+AR4+AR5))", "nup*AR4",
       "t5*(nA0-(AR1+AR2+AR3+AR4+AR5))", "nup*AR5")

# Characteristic time
chartime <- 1/(sum(theta)+nup)    # tau

# simulations of ML parameter estimation
simdata <- c()
numtimes <- 100          # Number of temporal repetitions in each data set

# Interval of sampling times (out-of-equilibrium sampling: [0, 10*tau])
# Comment the following three lines and uncomment the next three lines for steady-state sampling
tini <- 0.01
tend <- 10*chartime
flag <- 0                # Out-of-equilibrium sampling

# Interval of sampling times (steady_state sampling: [10*tau,20*tau])
#tini <- 10*chartime
#tend <- 20*chartime
#flag <- 1                # Steady-state sampling

tvec <- runif(numtimes,tini,tend)

for (i in 1:numtimes){
  
  tf <- tvec[i]
  out <- ssa(
    x0 = x0,
    a = a,
    nu = nu,
    parms = parms,
    tf = tf,
    method = ssa.d(),
    simName = simName,
    verbose = FALSE,
    consoleInterval = 1
  )
  
  # Keep only the number of free predators to estimate parameters
  simdata <- cbind(simdata,c(tail(out$data,n=1)[1,1:2]))
}

# Log-likelihood heatmap (loop in alpha and nu)

logL <- c()
alpha_vec <- seq(2.5,7.5,length.out=200)
nu_vec <- seq(7.5,12.5,length.out=200)

for (nu0 in nu_vec){
  row <- c()
  
  for (alpha0 in alpha_vec){

    row <- c(row, lkopt_binom(data = list(S,nA0,xR,simdata),par = c(alpha0,nu0)))
  }
  
  logL <- rbind(logL, row)
}

heat_df = data.frame(nu=rep(nu_vec, ncol(logL)),
                     alpha=rep(alpha_vec, each=nrow(logL)),
                     logL=c(logL))

# Figure

library(ggpubr)
library(ggsci)
library(ggplot2)
library(grid)
library(gridExtra)
library(latex2exp)
library(patchwork)
library(ggExtra)
library(cowplot)
library(viridis)
library(lattice)
library(viridisLite)
library("sp")
library("latticeExtra")
library("unmarked")

source("mytheme_ggplot.R")

coul <- viridis(100)

# Out-of-equilibrium sampling: single optima

if (flag == 0) {
  
  # Exact parameter values
  x <- c(5)
  y <- c(10)
  name <- c("exact")
  dummy <- data.frame(x, y, name)
  coordinates(dummy) <- ~ x + y
  
  # Estimate optimal values from log-lokelihood
  opt <- optim(par = c(1,1), fn = lkopt_binom, data = list(S,nA0,xR,simdata), 
               method = "L-BFGS-B", lower=c(0.001,0.001))
  
  x <- opt$par[1]
  y <- opt$par[2]
  name <- c("estimated")
  dummy1 <- data.frame(x, y, name)
  coordinates(dummy1) <- ~ x + y
  
  # Plot heatmap as well as optimal and estimated values
  g <- levelplot(logL~alpha*nu, data=heat_df, col.regions = coul, contour = TRUE, 
                  cuts = 15, font = 2,
                  xlab = list(expression(paste('Attack rate, ',symbol('a'))), fontsize=16),
                  ylab = list(expression(paste('Relaxation rate, ',symbol('n'))), fontsize=16),
                  scales=list(tck=c(1,0), x=list(cex=1), y=list(cex=1)),
                  colorkey = list(labels=list(cex=1))) +
    layer(sp.points(dummy, pch = 1, cex = 0.8, col = "white")) +
    layer(sp.points(dummy1, pch = 19, cex = 0.8, col = "white")) +
    layer(panel.arrows(5.5, 10, 5, 10, length = 0.1, col = 'white')) 

  print(g)
}

# Steady-state sampling

if (flag == 1) {
  
  # Estimate line of degenerate optima
  data_hist <- c()
  
  for (j in 1:500){
    
    # Optimization algorithm for varying initial conditions
    alpha0 <- runif(1,0,20)
    nu0 <- runif(1,0,20)
    
    opt <- optim(par = c(alpha0, nu0), fn = lkopt_binom, data = list(S,nA0,xR,simdata), 
                 method = "L-BFGS-B", lower=c(0.001,0.001))

    data_hist <- rbind(data_hist,c(opt$par,opt$value))
  }

  hist_df <- as.data.frame(data_hist)
  names(hist_df)<-c("alpha","nu","-logL")
  
  # Estimate line from linear model
  linmod <- lm(nu ~ alpha, data = hist_df)
  
  # Plot heatmap as well as the line of optimal values
  g <- levelplot(logL~alpha*nu, data=heat_df, col.regions = coul, contour = TRUE, 
                 cuts = 15, font = 2,
                 xlab = list(expression(paste('Attack rate, ',symbol('a'))), fontsize=16),
                 ylab = list(expression(paste('Relaxation rate, ',symbol('n'))), fontsize=16),
                 scales=list(tck=c(1,0), x=list(cex=1), y=list(cex=1)),
                 colorkey = list(labels=list(cex=1))) +
    layer(panel.abline(0,linmod$coefficients[2], col="white"))
  
  print(g)
}

