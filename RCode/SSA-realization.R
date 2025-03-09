rm(list = ls(all.names = TRUE))

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
chartime <- 1/(sum(theta)+nup)

# Final time
tf <- 10*chartime

# Gillespie SSA single realization
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

# Plot trajectories
ssa.plot(out, show.title = TRUE, show.legend = FALSE)

