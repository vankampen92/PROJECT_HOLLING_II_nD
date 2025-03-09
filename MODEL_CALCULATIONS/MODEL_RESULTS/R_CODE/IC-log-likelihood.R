rm(list = ls(all.names = TRUE))

###
### Log-likelihood profile computation
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

# Optimization function for profiles: alpha
lkopt_prof_alpha <- function(data, par){
  nup <- par[1]
  S <- data[1][[1]]
  nA0 <- data[2][[1]]
  xR <- data[3][[1]]
  dataset <- data[4][[1]]
  alphap <- data[5][[1]]
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

# Optimization function for profiles: alpha
lkopt_prof_nu <- function(data, par){
  alphap <- par[1]
  S <- data[1][[1]]
  nA0 <- data[2][[1]]
  xR <- data[3][[1]]
  dataset <- data[4][[1]]
  nup <- data[5][[1]]
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

tini <- 0.01
tend <- 10*chartime

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
  
  simdata <- cbind(simdata,c(tail(out$data,n=1)))
}

# Optimal estimates using the complete log-likelihood
opt <- optim(par = c(1,1), fn = lkopt, data = list(S,nA0,xR,simdata), 
             method = "L-BFGS-B", lower=c(0.001,0.001))

# Log-likelihood profile for parameter alpha
logL <- c()
alphaseq <- seq(4,7,length.out=200)

for (alphap in alphaseq){
  
  opt1 <- optim(par = c(1), fn = lkopt_prof_alpha, data = list(S,nA0,xR,simdata,alphap), 
                method = "L-BFGS-B", lower=c(0.001,0.001))
  aux <- lkopt(c(alphap,opt1$par),data = list(S,nA0,xR,simdata))
  logL <- c(logL,aux)
}

plot_df <- data.frame(cbind(alphaseq,logL,rep("alpha",length(alphaseq))))
names(plot_df) <- c("x","y","ind")
plot_df$x <- as.numeric(plot_df$x)
plot_df$y <- as.numeric(plot_df$y)

# Evaluate loh-likelihood to define CI based on the chi-square distribution
logL_min <- lkopt(opt$par,data = list(S,nA0,xR,simdata))
idx <- which(logL<logL_min+qchisq(0.95,1)/2 & logL>logL_min)

# Log-likelihood profile for parameter nu
logL <- c()
nuseq <- seq(8,14,length.out=200)

for (nup in nuseq){
  
  opt1 <- optim(par = c(1), fn = lkopt_prof_nu, data = list(S,nA0,xR,simdata,nup), 
                method = "L-BFGS-B", lower=c(0.001,0.001))
  aux <- lkopt(c(opt1$par,nup),data = list(S,nA0,xR,simdata))
  logL <- c(logL,aux)
}

plot_df1 <- data.frame(cbind(nuseq,logL,rep("nu",length(alphaseq))))
names(plot_df1) <- c("x","y","ind")
plot_df1$x <- as.numeric(plot_df1$x)
plot_df1$y <- as.numeric(plot_df1$y)

plot_df <- rbind(plot_df,plot_df1)

# Evaluate loh-likelihood to define CI based on the chi-square distribution
logL_min <- lkopt(opt$par,data = list(S,nA0,xR,simdata))
idx <- which(logL<logL_min+qchisq(0.95,1)/2 & logL>logL_min)

###
### Figure
###

# Set vertical lines at the minima of negative log likelihood
vline1 <- as.data.frame(t(opt$par))
names(vline1) <- c("alpha", "nu")
vline.data <- stack(vline1) %>%
  group_by(ind) %>%
  summarize(z = mean(values))

# Set horizontal lines at the extremes of confidence intervals
hline1 <- as.data.frame(t(c(logL_min+qchisq(0.95,1)/2,logL_min+qchisq(0.95,1)/2)))
names(hline1) <- c("alpha", "nu")
hline.data <- stack(hline1) %>%
  group_by(ind) %>%
  summarize(z = mean(values))

# Plot
source("mytheme_ggplot.R")

labels <- c(`alpha` = "Attack rate",`nu` = "Relaxation rate") 
g <- ggplot(plot_df, aes(x=x, y=y, color=factor(ind))) +
  geom_line(linewidth = 1) + mytheme +
  ylab(TeX(r'(Negative log-likelihood)')) +
  geom_vline(aes(xintercept = z), vline.data, colour = "#2980B9", linewidth = 1, linetype = "dashed") +
  geom_hline(aes(yintercept = z), hline.data, colour = "#2980B9", linewidth = 1, linetype = "dashed") +
  facet_wrap(. ~ ind, scales = "free_x", labeller = labeller(ind = labels), strip.position = "bottom") +
  xlab(NULL) +
  theme(axis.title=element_text(size=16),legend.position = "none",
        strip.background = element_blank(),strip.placement = "outside",
        strip.text.x=element_text(size=16))
g

