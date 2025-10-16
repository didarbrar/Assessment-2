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
n<-1000
##here we create a h vector of length and which contains the household ID for each person

## in single line by help of rep and sample

## maximum number of people in  household is five

## here we set.seed to generates a same number,Every time we run the code.

set.seed(9783)
hmax<-5
h<-rep(c(1:n),sample(1:hmax,n,replace=TRUE))[1:n]
length(h)
head(h,50)
h[1000]

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
for (i in b){
  print(i)
}
nseir<-function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005){
  
  n<-length(beta)
  
  x <- rep(0,n) ## initialize to susceptible state
  
  ni<-pinf*n
  
  B<-sum(beta)/n
  
  for (i in 1:ni) {
    x[sample(1:n,1)]<-2
  }
  
  S<-E<-I<-R<-rep(0,nt) ## set up storage for pop in each state
  S[1] <- n-ni;I[1] <- ni ## initialize
  ##loop over days 
  for(i in 2:nt) {
    u <- runif(n) ## uniform random deviates
    b<-which(x==2)
    x[x==2&u<delta] <- 3 ## I -> R with prob delta
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
    a<-which(x==0)
    for (i in a){
      IH<-0
      for (j in b){
        if(h(i)==h(j)){
          IH<-IH+1
        }
      }
      a1<-IH*alpha[1]
      IN<-0
      for (j in b){
        if(alink[[i]][j]==j){
          In<-In+1
        }
      } 
      a2<-IN*alpha[2]
      a3<-0
      for (j in b){
        a3<-a3+((alpha[3]*nc * beta[i] * beta[j]) / ((B^2) * (n - 1)))
      }
      if (u<a1|u<a2|u<a3){
        x[i]<-1
      }
    }
    S[i] <- sum(x==0); E[i] <- sum(x==1)
    I[i] <- sum(x==2); R[i] <- sum(x==3)
  }
  list(S=S,E=E,I=I,R=R,beta=beta)
  ##seir
}


