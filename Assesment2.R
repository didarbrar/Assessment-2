###############################################################
# GROUP PROJECT - GROUP 23
# -------------------------------------------------------------
# Members:
#   Member 1: Didar Singh (s2903783)
#   Member 2: Liyan Yang()
#
# Contribution:
#   • Didar Singh – 50% (Steps 1, 2, 3: Generating nSEIR function and contact network)
#   • Liyan Yang  – 50% (Steps 4, 5: Plot function and scenario visualizations)
#
# Description:
#   This project simulates the SEIR (Susceptible–Exposed–Infectious–Recovered)
#   epidemic model for a population using both structured and random contacts.
#   The `nseir()` function runs the simulation, while `plot_SEIR()` visualizes
#   the results along with a histogram of individual contact rates (β).
###############################################################



## step-1

#-----------------------


## n=total number of people 
n<-10000

##here we create a h vector of length and which contains the household ID for each person

## in single line by help of rep and sample

## maximum number of people in  household is five

## here we set.seed to generates a same number,Every time we run the code.

set.seed(9983)
hmax<-5
h<-rep(c(1:n),sample(1:hmax,n,replace=TRUE))[1:n]
length(h)

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
  stopifnot(length(h) == n, length(alink) == n)  # basic input checks
  alpha_h <- alpha[1]; alpha_c <- alpha[2]; alpha_r <- alpha[3]
  
  # State codes: 0 = Susceptible, 1 = Exposed, 2 = Infectious, 3 = Recovered
  
  x <- rep(0, n)   # initialize all individuals as Susceptible (0)
  
  # Seed initial infections randomly: at least 1 infected
  
  ni <- max(1, round(pinf * n))  # initial number of infected
  B <- mean(beta)                 # mean sociability for random mixing
  x[sample(1:n, ni)] <- 2         # randomly infect initial individuals
  
  # Storage vector for daily counts of S,E,I,R
  S <- E <- I <- R <- rep(0, nt)
  S[1] <- sum(x == 0); E[1] <- sum(x == 1)
  I[1] <- sum(x == 2); R[1] <- sum(x == 3)
  
  # Main simulation loop over days
  for(day in 2:nt){
    
    inf_idx <- which(x == 2)  # indices of infectious individuals
    n_inf <- length(inf_idx)
    
    # --- Progression transitions (state changes that do not depend on susceptibles) ---
    # I -> R : infectious individuals recover with probability 'delta' each day
    u <- runif(n)
    x[x==2 & u < delta] <- 3
    
    # E -> I : exposed individuals become infectious with probability 'gamma' each day
    u <- runif(n)
    x[x==1 & u < gamma] <- 2
    
    # Identify susceptibles for possible new infections this day
    S_idx <- which(x == 0)
    
    if(length(S_idx) > 0 && n_inf > 0){
      
      # Precompute coefficient used in the pairwise random-mixing infection probability.
      # This term scales probabilities so that average contact rate / force of infection is controlled
      # by alpha_r and nc, and normalizes by mean(beta)^2 and (n-1).
      coef_pair <- (alpha_r * nc) / (B^2 * (n - 1))
      
      # Loop over susceptibles
      for(i in S_idx){
        
        # --------------------
        # 1) Household infection
        # --------------------
        # Count number of infectious household members (same household id)
        IH <- sum(h[inf_idx] == h[i])
        # If there are IH infectious members, the probability of infection from household is:
        # probability of not infecting by any of one is (1 - alpha_h)^IH
        # so probability of infecting by atleat one of them is 
        # p_h = 1 - (1 - alpha_h)^IH   (assuming independent per-member transmission attempts)
        p_h <- if(IH > 0) 1 - (1 - alpha_h)^IH else 0
        
        # --------------------
        # 2) Regular network infection (alink)
        # --------------------
        # alink[[i]] lists the neighbors of i; IN = number of those neighbors who are infectious
        nb <- alink[[i]]
        IN <- if(length(nb) == 0) 0 else sum(nb %in% inf_idx)
        # Per-contact independent risk assumption as above
        p_c <- if(IN > 0) 1 - (1 - alpha_c)^IN else 0
        
        # --------------------
        # 3) Random mixing infection (mass-action / pairwise)
        # --------------------
        # For random mixing we compute an infection probability from each infectious j:
        #   p_ij = coef_pair * beta[i] * beta[j]
        # This models heterogeneity via individual betas (sociability).
        # Then the probability that i avoids infection from *all* infectives is prod(1 - p_ij),
        # so infection probability is 1 - product(...).
        p_ij <- coef_pair * beta[i] * beta[inf_idx]  # vector over infectious individuals
        # Clip to [0,1] for numerical safety before taking product
        p_ij[p_ij < 0] <- 0; p_ij[p_ij > 1] <- 1
        survive <- prod(1 - p_ij)  # prob of avoiding infection from every infectious individual
        p_r <- 1 - survive
        
        # --------------------
        # Combine independent sources of infection
        # --------------------
        # Assuming independence between household, network, and random-mixing infection attempts:
        #   probability of staying uninfected = (1-p_h)*(1-p_c)*(1-p_r)
        # so total infection probability from atleast one contact is  = 1 - product(...)
        p_total <- 1 - ((1 - p_h) * (1 - p_c) * (1 - p_r))
        p_total <- min(max(p_total, 0), 1)  # clip to [0,1] to avoid tiny numerical errors
        
        # Final Bernoulli trial: if infected, move susceptible -> exposed (0 -> 1)
        if(runif(1,0,1) < p_total) x[i] <- 1
      }
      
    } # end susceptibles loop
    
    # Store daily compartment counts for this day
    S[day] <- sum(x == 0); E[day] <- sum(x == 1)
    I[day] <- sum(x == 2); R[day] <- sum(x == 3)
  }
  
  # Return results and the input beta for reference
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

# ==========================================================
# Function: plot_SEIR
# Purpose : Visualize SEIR results along with beta histogram
# ----------------------------------------------------------
# Inputs:
#   sim   : output list from nseir() containing S,E,I,R,beta
#   main1 : title for histogram plot
#   main2 : title for SEIR scatter plot
#
# Output:
#   Produces two plots:
#     (1) Histogram of beta (contact rate distribution)
#     (2) Scatter plot showing SEIR dynamics over time
# ==========================================================
plot_SEIR <- function(sim,
                      main1="Distribution of Beta (Contact Rates)",
                      main2="SEIR Dynamics (Scatter Plot)") {
  
  # ensure required components are present
  stopifnot(all(c("S","E","I","R") %in% names(sim)))
  
  nt <- length(sim$S)    # number of days
  time <- 1:nt           # x-axis time vector
  
  # ---- Step 1: Plot histogram of beta ----
  hist(sim$beta, col="skyblue", border="white",
       main=main1,
       xlab=expression(beta), ylab="Count")
  grid()
  
  # ---- Step 2: Create SEIR scatter plot ----
  plot(time, sim$S, type="n",
       ylim=c(0, max(sim$S, sim$E, sim$I, sim$R)),
       xlab="Day", ylab="Number of individuals",
       main=main2)
  
  # add scatter points for each compartment
  points(time, sim$S, col="black")   # Susceptible
  points(time, sim$E, col="orange")  # Exposed
  points(time, sim$I, col="red")     # Infectious
  points(time, sim$R, col="green")   # Recovered
  
  # add legend for color coding
  legend("topright", legend=c("S","E","I","R"),
         col=c("black","orange","red","green"),
         pch=16, bty="n")
  
  # add grid lines for better readability
  grid()
}


# ==========================================================
# Example: Simulate and plot one SEIR scenario
# ----------------------------------------------------------
# This example demonstrates how to:
#   1) Define population and contact structure
#   2) Run nseir() to simulate epidemic spread
#   3) Plot results using plot_SEIR()
# ==========================================================

par(mfcol=c(1,1), mar=c(4,4,2,1))  # set single plot layout

# ---- Step 1: Define model inputs ----
set.seed(9983)     # reproducibility
n <- 1000          # total population
beta <- runif(n, 0, 1)  # individual contact rates
hmax <- 5          # maximum household size

# assign household IDs
h <- rep(c(1:n), sample(1:hmax, n, replace=TRUE))[1:n]

# construct social contact network
alink <- get.net(beta, h, nc=30)

# ---- Step 2: Create constant-beta comparison scenario ----
constant_beta <- rep(mean(beta), length(beta))

# ---- Step 3: Run SEIR simulation ----
epi_4 <- nseir(constant_beta, h, alink, alpha=c(0.1,0.01,0.01))

# ---- Step 4: Plot SEIR results ----
plot_SEIR(epi_4)

# Output interpretation:
#   - Top: Histogram showing contact rate distribution (beta)
#   - Bottom: Scatter plot of susceptible, exposed, infectious, recovered counts over time
#   The plot allows visual comparison of how the population tra



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

plot_SEIR(epi_0,main2="Default full model")
# ------------------------------------------------
# Scenario 1: Remove household/network, random mixing
# alpha_h = alpha_c = 0, alpha_r = 0.04
# ------------------------------------------------

epi_1 <- nseir(beta, h, alink, alpha=c(0,0,0.04))

plot_SEIR(epi_1,main2="compared_model1(only random mixing)")

# ------------------------------------------------
# Scenario 2: Constant beta, full model
# ------------------------------------------------
constant_beta <- rep(mean(beta), length(beta))
epi_2 <- nseir(beta=constant_beta, h, alink)

plot_SEIR(epi_2,main2="compared_model2(constant beta)")

# ------------------------------------------------
# Scenario 3: Constant beta, random mixing
# alpha_h = alpha_c = 0, alpha_r = 0.04
# ------------------------------------------------
epi_3 <- nseir(beta=constant_beta, h, alink, alpha=c(0,0,0.04))

plot_SEIR(epi_3,main2="compared_model2(constant beta,only random mixing)")






