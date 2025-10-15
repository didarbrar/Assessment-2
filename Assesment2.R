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
  
  
  
