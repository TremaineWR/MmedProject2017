N0 <- 23560000   #Population size of the Republic of Congo in 1976

# beta.calc <- function(Ro,mu=1/0.99,sigma=1/5.99,gamma=0.1){
#   Ro/((sigma/(mu+sigma))*(1/(mu+gamma)))
# }

library(deSolve)

# Create parameters setup of ODEs for simulation
seir <- function(t,y,params){
  S <- y[1]
  E <- y[2]
  I <- y[3]
  H <- y[4]
  R <- y[5]
  
  beta <- params["beta"]
  N <- params["N"]
  mu <- params["mu"]
  sigma <- params["sigma"]
  gamma <- params["gamma"]
  kappa <- params["kappa"]
  lambda <- params["lambda"]
  #nu <- mu*N
  
  
  dSdt <- -beta*I/N
  dEdt <- beta*I/N-sigma*E
  dIdt <- sigma*E-mu*I-gamma*I-lambda*I
  dHdt <- lambda*I - kappa*H
  dRdt <- gamma*I + kappa*H
  
  
  return(list(c(dSdt,dEdt,dIdt,dHdt,dRdt)))
}

# Value of beta that yields a reproduction number of 1.3 used in the project.
# Initialisation of remaining parameters.
param.vals <- c(
  beta=10.1748,
  
  N=N0,
  mu=1/0.99,
  sigma=1/5.99,
  gamma=0.1,
  kappa = 0.21,
  lambda = 1/3
)

# Time interval for plots
times <- seq(0,60,1)

#Initial conditions.
S0 <- 20000
E0 <- 5000
I0 <- 255
H0 <- 150   #Isolate 150 infected individuals

# Vector of initial conditions
init <- c(sus=S0,exp=E0,inf=I0,H0,rec=N0-S0-E0-I0-H0)

# Data frame to show the changes in the individuals in the different compartments.
tc <- data.frame(lsoda(init,times,seir,param.vals))
par('ps'=16,lwd=2)

# Plot of the number of cases against time.
plot(tc$time,tc$inf,type="l",xlab="Time (days)",ylab="Number of cases",main = "With isolation",
     bty="n",ylim = c(0,600))
