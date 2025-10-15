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

## step-2 ___________________

##-------------------

get.net<-function(Beta,h,nc){
  mylist<-vector("list",n)
  B<-sum(Beta)/n
  for (i in  1:n){
    link<-c()
    for (j in 1:n){
      u<-runif(1,0,1)
      p=(nc*Beta[i]*Beta[j])/((B**2)*(n-1))
      if(h[i]!=h[j]&u<p){
        link<-c(link,j)
      }
    }
    mylist[[i]]<-link
  }
  return (mylist)
}

set.seed(42)
n <- 10
Beta <- runif(n, 0.5, 1.5)
h <- c(1,1,2,2,3,3,4,4,5,5)

contacts <- get.net(Beta, h, nc = 15)
contacts


