N0 <- 23560000 #Population size of the Republic of Congo in 1976

beta.calc <- function(Ro,mu=1/0.99,sigma=1/5.99,gamma=0.1){
  Ro/((sigma/(mu+sigma))*(1/(mu+gamma)))
}

library(deSolve)

# Create parameters setup of ODEs for simulation
seir <- function(t,y,params){
  S <- y[1]
  E <- y[2]
  I <- y[3]
  R <- y[4]
  
  beta <- params["beta"]
  N <- params["N"]
  mu <- params["mu"]
  sigma <- params["sigma"]
  gamma <- params["gamma"]
  #nu <- mu*N
  
  
  dSdt <- -beta*I/N
  dEdt <- beta*I/N-sigma*E 
  dIdt <- sigma*E-mu*I-gamma*I
  dRdt <- gamma*I
  
  return(list(c(dSdt,dEdt,dIdt,dRdt)))
}

# Value of beta that yields a reproduction number of 1.3 used in the project.
# Initialisation of remaining parameters.
param.vals <- c(
  beta=beta.calc(1.3),

  N=N0,
  mu=1/0.99,
  sigma=1/5.99,
  gamma=0.1
)

# Time interval for plots
times <- seq(0,60,1)

#Initial conditions.
S0 <- 200000
E0 <- 5000
I0 <- 255
#R0 <- 0

# Vector of initial conditions
init <- c(sus=S0,exp=E0,inf=I0,rec=N0-S0-E0-I0)

# Data frame to show the changes in the individuals in the different compartments.
tc <- data.frame(lsoda(init,times,seir,param.vals))
par('ps'=16,lwd=2)

# Plot of the number of cases against time.
plot(tc$time,tc$inf,type="l",xlab="Time (days)",ylab="Number of cases",main = "Without isolation",
     bty="n",ylim = c(0,600))

