# Name=Didar Singh

# s2903783  

#here is the function to implement SEIR model simulation with diffrent beta

##here we have a function which have inputs like a number of people and 
##initial defective and the number of days, and then we have  beta parameter parameters

seir<-function(n=10000,ni=10,nt=100,gamma=1/3,delta=2/3,bmu=0.5,bsc=0.05){
  ## SEIR stochastic simulation model.
  ## n = population size; ni = initially infective; nt = number of days
  ## gamma = daily prob E -> I; delta = daily prob I -> R;
  ## bmu = mean beta; bsc = var(beta) = bmu * bsc
  
  x<-rep(0,n) ## initialize to susceptible state
  x[1:ni]<-2 ## create some infectives
  beta <- rgamma(n,shape=bmu/bsc,scale=bsc) ## individual infection rates
  S<-E<-I<-R<-rep(0,nt) ## set up storage for pop in each state
  S[1] <- n-ni;I[1] <- ni ## initialize
  ##loop over days 
  for(i in 2:nt) {
    u <- runif(n) ## uniform random deviates
    x[x==2&u<delta] <- 3 ## I -> R with prob delta
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
    x[x==0&u<beta*I[i-1]] <- 1 ## S -> E with prob beta*I[i-1]
    S[i] <- sum(x==0); E[i] <- sum(x==1)
    I[i] <- sum(x==2); R[i] <- sum(x==3)
  }
  list(S=S,E=E,I=I,R=R,beta=beta)
 ##seir   
}
seir(n=10000,ni=10,nt=100,gamma=1/3,delta=1/5,bmu=5e-5,bsc=1e-5)$S

## step-1

#-----------------------


## n=total number of people 
n<-10000

##here we create a h vector of length and which contains the household ID for each person

## in single line by help of rep and sample

## maximum number of people in  household is five

## here we set.seed to generates a same number,Every time we run the code.

set.seed(9783)
hmax<-5
h<-rep(c(1:n),sample(1:hmax,n,replace=TRUE))[1:n]
length(h)
head(h,50)


##---------------------------------------------------------------
## Function: get.net
##---------------------------------------------------------------
## Purpose:
##   Create a network of non-household contacts based on 
##   individuals’ activity levels (Beta) and household structure (h).
##
## Inputs:
##   Beta : vector of activity levels for each person
##   h    : vector indicating household membership
##   nc   : expected number of non-household contacts per person
##
## Output:
##   A list of length n, where each element contains the indices 
##   of that person’s non-household contacts.
##---------------------------------------------------------------

get.net<-function(beta,h,nc){
  ##Number of individuals in population
  n<-length(beta)
  
  # Create an empty list where each element will store contacts for one person
  mylist<-vector("list",n)
  
  # Initialize each list element as an empty integer vector
  for(i in 1:n) mylist[[i]]<-integer(0)
  
  # Mean activity level across the population (used in probability formula)
  B<-sum(beta)/n
  
  
  # Loop through all unique pairs (i, j) with i < j to avoid double checking
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      # Only consider non-household contacts
      if (h[i] != h[j]) {
        
        # Compute probability of forming a link between i and j
        p <- (nc * beta[i] * beta[j]) / ((B^2) * (n - 1))
        
        # Ensure probability does not exceed 1
        p <- min(1, p)
        
        # Draw one random number to decide if link is formed
        if (runif(1, 0, 1) < p) {
          
          # Record the link in both individuals' contact lists (symmetric)
          mylist[[i]] <- c(mylist[[i]], j)
          mylist[[j]] <- c(mylist[[j]], i)
        }
      }
    }
  }
  
  # Return the full network as a list
  return(mylist)
}

##----------------------

## step-3 

##--------------------

## function for finding number of prople in S,E,I,R for each day
x<-c(0,0,1,1,2,3,4,2)
b<-(which(x==2))
b

# nseir: Individual-based SEIR epidemic model with household, network, and random mixing
#
# Inputs:
#   beta   : numeric vector of length n, sociability/activity of each individual
#   h      : integer/vector of length n, household IDs
#   alink  : list of length n, each element contains regular contacts (vector of indices)
#   alpha  : vector of 3 infection probabilities c(alpha_h, alpha_c, alpha_r)
#            alpha_h: household, alpha_c: network, alpha_r: random mixing
#   delta  : I -> R daily probability
#   gamma  : E -> I daily probability
#   nc     : average community contacts
#   nt     : number of simulation days
#   pinf   : initial fraction infected (0-1)
#
# Output:
#   list with elements:
#     S, E, I, R : daily counts of susceptible, exposed, infectious, recovered
#     beta       : returned sociability vector for reference

nseir <- function(beta, h, alink, alpha=c(.1,.01,.01), delta=.2, gamma=.4, nc=15, nt=100, pinf=.005){
  
  n <- length(beta)  # population size
  stopifnot(length(h) == n, length(alink) == n)  # basic checks
  alpha_h <- alpha[1]; alpha_c <- alpha[2]; alpha_r <- alpha[3]
  
  x <- rep(0, n)   # initialize all individuals as Susceptible (0)
  
  ni <- max(1, round(pinf * n))  # initial number of infected
  B <- mean(beta)                 # mean sociability for random mixing
  x[sample(1:n, ni)] <- 2         # randomly infect initial individuals
  
  # Storage for daily counts
  S <- E <- I <- R <- rep(0, nt)
  S[1] <- sum(x == 0); E[1] <- sum(x == 1)
  I[1] <- sum(x == 2); R[1] <- sum(x == 3)
  
  # Main simulation loop over days
  for(day in 2:nt){
    
    inf_idx <- which(x == 2)  # indices of infectious individuals
    n_inf <- length(inf_idx)
    
    # --- Progression transitions ---
    # I -> R (recovery)
    u <- runif(n)
    x[x==2 & u < delta] <- 3
    
    # E -> I (incubation to infectious)
    u <- runif(n)
    x[x==1 & u < gamma] <- 2
    
    # Susceptibles
    S_idx <- which(x == 0)
    
    if(length(S_idx) > 0 && n_inf > 0){
      
      # Precompute coefficient for random mixing term
      coef_pair <- (alpha_r * nc) / (B^2 * (n - 1))
      
      # Loop over susceptibles
      for(i in S_idx){
        
        # 1) Household infection
        IH <- sum(h[inf_idx] == h[i])  # number of infectious household members
        p_h <- if(IH > 0) 1 - (1 - alpha_h)^IH else 0
        
        # 2) Network infection
        nb <- alink[[i]]
        IN <- if(length(nb) == 0) 0 else sum(nb %in% inf_idx)
        p_c <- if(IN > 0) 1 - (1 - alpha_c)^IN else 0
        
        # 3) Random mixing infection
        p_ij <- coef_pair * beta[i] * beta[inf_idx]  # infection prob from each infectious individual
        p_ij[p_ij < 0] <- 0; p_ij[p_ij > 1] <- 1    # clip probabilities
        survive <- prod(1 - p_ij)  # probability of avoiding infection from all infectives
        p_r <- 1 - survive
        
        # --- Combine all sources of infection ---
        p_total <- 1 - ((1 - p_h) * (1 - p_c) * (1 - p_r))
        p_total <- min(max(p_total, 0), 1)  # safety for numerical errors
        
        # Infection trial
        if(runif(1,0,1) < p_total) x[i] <- 1
      }
      
    } # end susceptibles loop
    
    # Store daily counts
    S[day] <- sum(x == 0); E[day] <- sum(x == 1)
    I[day] <- sum(x == 2); R[day] <- sum(x == 3)
  }
  
  # Return results
  list(S = S, E = E, I = I, R = R, beta = beta)
}

# plot_SEIR: Plot SEIR dynamics from nseir simulation
#
# Input:
#   sim  : list returned by nseir() containing S, E, I, R vectors
#   main : optional title for the plot
#
# Output:
#   Generates a line plot of S, E, I, R over time

plot_SEIR<- function(sim, main="SEIR Dynamics (Scatter Plot)") {
  
  # Check input
  stopifnot(all(c("S","E","I","R") %in% names(sim)))
  
  nt <- length(sim$S)          # number of days
  time <- 1:nt
  
  # Create an empty plot with proper y-limits
  plot(time, sim$S, type="n", ylim=c(0, max(sim$S, sim$E, sim$I, sim$R)),
       xlab="Day", ylab="Number of individuals", main=main)
  
  # Add scatter points for each compartment (all as round circles)
  points(time, sim$S, col="black", pch=16)   # S: blue
  points(time, sim$E, col="orange", pch=16) # E: orange
  points(time, sim$I, col="red", pch=16)    # I: red
  points(time, sim$R, col="green", pch=16)  # R: green
  
  
  # Add legend
  legend("topright", legend=c("S","E","I","R"), col=c("blue","orange","red","green"), pch=16, bty="n")
  
  # Optional: add grid for clarity
  grid()
}



# Example: simulate one scenario with nseir
set.seed(9983)
n <- 1000
beta <- runif(n, 0, 1)
hmax <- 5
h<-rep(c(1:n),sample(1:hmax,n,replace=TRUE))[1:n]
alink <- get.net(beta, h, nc=15)

# Set plot window: 2 rows, 4 columns, with margins
par(mfcol=c(2,4), mar=c(4,4,2,1))

# -------------------------------
# Scenario 0: Full model (default)
# -------------------------------
epi_0 <- nseir(beta, h, alink)

# Plot beta distribution
hist(epi_0$beta, xlab="beta", main="Scenario 0: Full model")

plot_SEIR(epi_0)
# ------------------------------------------------
# Scenario 1: Remove household/network, random mixing
# alpha_h = alpha_c = 0, alpha_r = 0.04
# ------------------------------------------------
epi_1 <- nseir(beta, h, alink, alpha=c(0,0,0.04))

hist(epi_1$beta, xlab="beta", main="Scenario 1: Random mixing")

plot_SEIR(epi_1)


# ------------------------------------------------
# Scenario 2: Constant beta, full model
# ------------------------------------------------
constant_beta <- rep(mean(beta), length(beta))
epi_2 <- nseir(beta=constant_beta, h, alink)

hist(epi_2$beta, xlab="beta", main="Scenario 2: Constant beta")

plot_SEIR(epi_2)


# ------------------------------------------------
# Scenario 3: Constant beta, random mixing
# alpha_h = alpha_c = 0, alpha_r = 0.04
# ------------------------------------------------
epi_3 <- nseir(beta=constant_beta, h, alink, alpha=c(0,0,0.04))

hist(epi_3$beta, xlab="beta", main="Scenario 3: Constant beta + random mixing")

plot_SEIR(epi_3)



