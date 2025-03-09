rm(list = ls(all.names = TRUE))

###
### Maximum likelihood estimation functions
###

# Multinomial probability density
pn <- function(x,xR,alpha,nu,nA0,t){
  theta <- alpha*xR
  p0 <- (1+exp(-(sum(theta)+nu)*t)*sum(theta)/nu)/(1+sum(theta)/nu)
  prob <- c(p0)
  aux <- (1-exp(-(sum(theta)+nu)*t))/(1+sum(theta)/nu)
  prob <- c(prob,theta*aux/nu)
  return(dmultinom(x,size=nA0,prob))
}

# Log-likelihood for a single time
lk_one_time <- function(x,xR,alpha,nu,nA0,t){
  return(-log(pn(x,xR,alpha,nu,nA0,t)))
}

# Additive log-likelihood for a dataset evaluated at different times
lk <- function(dataset,xR,alpha,nu,nA0){
  times <- dataset[1,]
  ntimes <- length(times)
  logL <- 0
  S <- length(dataset[,1])
  for (i in 1:ntimes){
    x <- dataset[2:S,i]
    logL <- logL + lk_one_time(x,xR,alpha,nu,nA0,times[i])
  }
  return(logL)
}

# Likelihood optimization function
lkopt <- function(data, par){
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
    x <- dataset[2:fin,i]
    logL <- logL + lk_one_time(x,xR,alphap,nup,nA0,times[i])
  }
  return(logL)
}

###
### Gillespie's Stochastic Simulation Algorithm
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
data_hist <- c()
numsim <- 100
numtimes <- 20          # Number of temporal repetitions in each data set

# Interval of sampling times (out-of-equilibrium sampling: [0, 10*tau])
tini <- 0.01
tend <- 10*chartime

for (j in 1:numsim){
  simdata <- c()
  
  # Sample end times in interval [tini, tend]
  tvec <- runif(numtimes,tini,tend)
  seqtimes <- sort(tvec, decreasing=FALSE)
  
  for (i in 1:numtimes){
    
    # Single SSA realization
    tf <- seqtimes[i]
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
    
    simdata <- cbind(simdata,c(tail(out$data,n=1)))
  }
  
  # Optimization algorithm to estimate alpha and nu
  alpha0 <- 5
  nu0 <- 10

  opt <- optim(par = c(alpha0, nu0), fn = lkopt, data = list(S,nA0,xR,simdata))
  
  # Append estimates to compute distributions
  data_hist <- rbind(data_hist,opt$par)
  
  print(paste("Estimation ", j, " completed"))
}

data_hist_df <- as.data.frame(data_hist)
names(data_hist_df)<-c("alpha","nu")

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
source("mytheme_ggplot.R")

# Dataframe for true values 
vline1 <- as.data.frame(t(c(5,10)))
names(vline1) <- c("alpha","nu")
vline.data <- stack(vline1) %>%
  group_by(ind) %>%
  summarize(z = mean(values))
par_names = c("Attack rate", "Relaxation rate")
ind_labeller <- function(variable,value){
  return(par_names[value])
}

# Histogram pre-processing
hist_df <- stack(data_hist_df)
hist_df <- cbind(hist_df,rep(100,nrow(hist_df)))
names(hist_df) <- c("data","ind","nA0")

# Plot
g1 <- ggplot(hist_df, aes(data, fill=factor(ind))) +
  geom_density(color="black", alpha=0.4) + mytheme +
  xlab(TeX(r'(Value)')) +
  ylab(TeX(r'(Density)')) +
  scale_fill_discrete(name = "", labels = c(TeX(r'($\alpha$)'),TeX(r'($\nu$)'))) +
  geom_vline(aes(xintercept = z), vline.data, colour = "paleturquoise4", linewidth = 0.7, linetype = "dashed") +
  facet_grid(. ~ nA0, scales = "free_x") +
  theme(axis.title=element_text(size=16))

hist_df_d <- cbind(data_hist_df,as.factor(rep("100",nrow(data_hist_df))))
names(hist_df_d)[3] <- "nA0"

g2 <- ggplot(hist_df_d) +
  stat_density_2d(aes(x = alpha, y = nu, fill = ..level..),
                  geom = "polygon", contour_var = "ndensity") + mytheme +
  scale_fill_viridis() +
  xlim(c(1.9999,8.0001)) +
  ylim(c(3.9999,17.0001)) +
  xlab(TeX(r'(Attack rate, $\alpha$)')) +
  ylab(TeX(r'(Relaxation rate, $\nu$)')) +
  facet_grid(. ~ nA0) +
  geom_vline(aes(xintercept = z), vline.data[1,], colour = "paleturquoise4", linewidth = 0.7, linetype = "dashed") +
  geom_hline(aes(yintercept = z), vline.data[2,], colour = "paleturquoise4", linewidth = 0.7, linetype = "dashed") +
  theme(axis.title=element_text(size=16), legend.title=element_blank())

g1/g2
