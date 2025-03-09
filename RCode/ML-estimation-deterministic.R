rm(list = ls(all.names = TRUE))

###
### Parameter estimation based on the deterministic model
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
numsim <- 500

# End time for stochastic realizations
tend <- 50*chartime

data_hist <- c()

for (j in 1:numsim){

  # Single SSA realization
  tf <- tend
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

  # Parameter alpha estimation as the slope at t = 0
  dt <- out$data[1:20,1:2]
  df <- as.data.frame(dt)
  names(df) <- c("t","n")
  df$n <- df$n/nA0
  fit<-lm(n~t,df)
  sfit <- summary(fit)
  
  # Estimate of parameter alpha
  alpha_est <- -fit$coefficients[2]/sum(xR)
  
  # Estimating the ration nu/alpha using steady-state samples
  steady <- out$data[500:nrow(out$data),]
  averages <- apply(steady[,2:(S+2)],2,mean)
  ratio <- averages[1]/(sum(averages[2:(S+1)]))*sum(xR)
  
  # Append estimates to compute distributions
  data_hist <- rbind(data_hist, c(alpha_est,ratio*alpha_est))
  
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
  xlim(c(0.0001,10.0001)) +
  ylim(c(0.0001,20.0001)) +
  xlab(TeX(r'(Attack rate, $\alpha$)')) +
  ylab(TeX(r'(Relaxation rate, $\nu$)')) +
  facet_grid(. ~ nA0) +
  geom_vline(aes(xintercept = z), vline.data[1,], colour = "paleturquoise4", linewidth = 0.7, linetype = "dashed") +
  geom_hline(aes(yintercept = z), vline.data[2,], colour = "paleturquoise4", linewidth = 0.7, linetype = "dashed") +
  theme(axis.title=element_text(size=16), legend.title=element_blank())

g1/g2
