library(deSolve)

# Figure 5. Example of isotopically non-stationary MFA. 
# We are calculating the labeling kinetics of FBP.
# The pool size of Gluc6P will affect how fast FBP get labeled.
# Eventually, when the labeling reaches steady-state, the Gluc6P pool size does not affect the labeling.


# function calculating the labeling pattern when combining x and y
CP <- function(x,y) {
  nx <- length(x)
  ny <- length(y)
  n <- nx+ny-1
  Tx <- matrix(0,ncol=ny,nrow=n)
  Product <- rep(0,n)
  for (i in 1:ny) {
    Tx[i:(i+nx-1),i] <- x 
  }
  Product <- Tx %*% y
  return(t(Product))
}

# The kinetics model for ODE
MODEL <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {

    #The state variables are the labeling patterns    
    G6P <- state[c(1:3)]
    FBP <- state[c(4:6)]
    DHAP <- state[c(7,8)]
    GAP <- state[c(9,10)]
    G6P.123 <- state[c(11,12)]
    G6P.456 <- state[c(13,14)]
    FBP.123 <- state[c(15,16)]
    FBP.456 <- state[c(17,18)]

    #The rate of change follows Equation 28
    dG6P <- 100*(Tracer-G6P)/parameters[1]
    dG6P.123 <- 100*(Tracer.123-G6P.123)/parameters[1]
    dG6P.456 <- 100*(Tracer.456-G6P.456)/parameters[1]
    dFBP <- (100*(G6P-FBP)+50*(CP(DHAP,GAP)-FBP))/parameters[2]
    dFBP.123 <- (100*(G6P.123-FBP.123)+50*(DHAP-FBP.123))/parameters[2]
    dFBP.456 <- (100*(G6P.456-FBP.456)+50*(GAP-FBP.456))/parameters[2]
    dDHAP <- (150*(FBP.123-DHAP)+150*(GAP-DHAP))/parameters[3]
    dGAP <- (150*(FBP.456-GAP)+250*(DHAP-GAP))/parameters[4]
    
    #Return the labeling pattern as the state variables
    return(list(c(dG6P,dFBP,dDHAP,dGAP,dG6P.123,dG6P.456,dFBP.123,dFBP.456)))
  })
}

# Describe the glucose tracer
# The tracer is 1,2-13C-glucose, which is M+2 labeled in C-123 and M+0 labeled in C-456 
Tracer <- c(0,1,0)
Tracer.123 <- c(0,1)
Tracer.456 <- c(1,0)

# Setting the initial labeling pattern
# The metabolites are Gluc6P, FBP, DHAP, GAP, Gluc6P.123, Gluc6P.456, FBP.123, FBP.456
# Each metabolite is described as c(M+0, M+2, M+4) or c(M+0, M+2)
init <- c(rep(c(1,0,0),2),rep(c(1,0),6))

# Setting the concentration of the metabolites
# The metabolites are Gluc6P, FBP, DHAP, GAP
# The unit is nmol
Concentration <- c(3000,1000,1000,1000)

# Setting the time range of the experiment. The total time is 200 seconds. 
times      <- seq(0, 200, by = 1)

# Calculating and plotting the kinetics of FBP labeling
out.slow <- ode(y=init, times=times, func=MODEL, parms=Concentration)
plot(times, out.slow[,5],col="blue",type="l",lwd=2,xlab="Time (s)", ylab="Labeled Fraction",ylim = c(0,1))
lines(x=times, y=out.slow[,6],col="green",lwd=2)
lines(x=times, y=out.slow[,7],col="orange",lwd=2)
legend(150,0.5,legend=c("M+0","M+2","M+4"), col=c("blue","green","orange"),ncol=1,lty=1,lwd=2)
title("FBP Labeling (Glc6P=3000 nmol)")


# Change the pool size of Gluc6P to be 400 nmols
Concentration <- c(400,1000,1000,1000)
# Calculating and plotting the kinetics of FBP labeling
out.fast <- ode(y=init, times=times, func=MODEL, parms=Concentration)
plot(times, out.fast[,5],col="blue",type="l",lwd=2,xlab="Time (s)", ylab="Labeled Fraction",ylim = c(0,1))
lines(x=times, y=out.fast[,6],col="green",lwd=2)
lines(x=times, y=out.fast[,7],col="orange",lwd=2)
legend(150,0.5,legend=c("M+0","M+2","M+4"), col=c("blue","green","orange"),ncol=1,lty=1,lwd=2)
title("FBP Labeling (Glc6P=400 nmol)")






