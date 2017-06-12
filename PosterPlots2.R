N0 <- 23560000

# beta.calc <- function(Ro,mu=1/0.99,sigma=1/5.99,gamma=0.1){
#   Ro/((sigma/(mu+sigma))*(1/(mu+gamma)))
# }

library(deSolve)

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

## Note that N0 is defined earlier in the script as 262.

param.vals <- c(
  beta=10.1748,
  
  N=N0,
  mu=1/0.99,
  sigma=1/5.99,
  gamma=0.1,
  kappa = 0.21,
  lambda = 1/3
)

times <- seq(0,60,1)

S0 <- 20000
E0 <- 5000
I0 <- 255
H0 <- 150
#R0 <- 0

## Stop and think about this. What is the level of immunity in
## the population for these initial conditions and a population
## size N0?

init <- c(sus=S0,exp=E0,inf=I0,H0,rec=N0-S0-E0-I0)

tc <- data.frame(lsoda(init,times,seir,param.vals))
par('ps'=16,lwd=2)
plot(tc$time,tc$inf,type="l",xlab="Time (days)",ylab="Infected population size",main = "With control measure",
     bty="n",ylim = c(0,600))
