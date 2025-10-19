# SOURCE CODE FOR IOPTOOFA

# libraries
library(combinat)
library(matrixcalc)
library(ggplot2)
library(reshape2)
############ Functions


### generates the full design with all m! permutations
generate_full_design <- function(m){
  
  Plist <- permn(1:m) ## get all permutations of 1:m
  return(matrix(unlist(Plist),byrow=TRUE,nrow=factorial(m)))
  
}


# design_to_PWO
# takes a design matrix D whose rows are permutations of (1:m)
# transforms to matrix of PWO variables
design_to_PWO <- function(D){
  m <- ncol(D)
  n <- nrow(D)
  q <- choose(m, 2)    #C^m_2
  X <- matrix(NA, nrow = n, ncol = q)
  cnames <- numeric(q+1)
  cnames[1] <- "Intercept"
  ind0 <- 2
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      cnames[ind0] <- paste("z", i, j, sep = "")
      ind0 <- ind0 + 1
    }
  }
  for(row in 1:n){
    index <- 1
    for(i in 1:(m-1)){
      for(j in (i+1):m){
        k <- which(D[row,] == i) #pairwise sign of each elecment of each row
        l <- which(D[row,] == j)
        X[row, index] <- sign(l-k)
        index <- index + 1
      }
    }
  }
  # cbind(1, X) adds a column of 1s for the intercept
  X <- cbind(1, X)
  colnames(X) <- cnames
  return(X)
  
}

trace <- function(X){
  return(sum(diag(X)))
}

I_efficiency <- function(X, B){ #X is the permutation matrix;
                                #B is the matrix to optimize
  M <- t(X)%*%X/nrow(X)
  if(is.singular.matrix(M)){ #singular matrix does not an inverse
    return(1e10)
  }
  
  return(trace(solve(M)%*%B))
  
}

# random_exchange
# randomly selects adjacent pairs of components to swap
# performs the swap if it improves I optimality
# TODO: Try to improve this. Move of the iterations are being wasted..
# inputs: 
# D: a n by m matrix of permutations of (1:m)
# B: the B matrix for the I opt criterion
# max_iter: max number of iterations, default is 1e4
# only neighbor swap is considered
random_exchange <- function(D, B, max_iter = 1e4){
  Dold <- D
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- design_to_PWO(Dold)
  Iold <- I_efficiency(Xold, B)
  delta <- Inf              # What is delta?
  Ieffs <- numeric(max_iter)
  for(iter in 1:max_iter){
    
    print(paste("Iteration", iter, sep = " "))
    
    
    i <- sample(1:n, 1)  # random row     sample(n,1)?
    j<- sample(1:(m-1),1) # random column     sample(m-1,1)?

    Dnew <- Dold
    # exchange component j with component (j+1)
    Dnew[i,j] <- Dold[i,j+1]
    Dnew[i,j+1] <- Dold[i, j]
    Xnew <- design_to_PWO(Dnew)
    Inew <- I_efficiency(Xnew, B)
    if(Inew < Iold){
      Iold <- Inew
      Dold <- Dnew
      Xold <- Xnew
    }
    Ieffs[iter] <- Iold
    
    
        
      
  }
  
  plot(1:max_iter, Ieffs, type = "l", ylab = "I-efficiency", xlab = "Iteration")
    
  return(Dold)
  
    
  }
  
### annealing simulation
annealing_exchange <- function(D, B, max_iter = 1e4){
  Dold <- D
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- design_to_PWO(Dold)
  Iold <- I_efficiency(Xold, B)
  # delta <- Inf
  Ieffs <- numeric(max_iter)
  for(iter in 1:max_iter){
    
    print(paste("Iteration", iter, sep = " "))
    
    
    i <- sample(1:n, 1)  # random row
    j<- sample(1:(m-1),1) # random column
    Dnew <- Dold
    # exchange component j with component (j+1)
    Dnew[i,j] <- Dold[i,j+1]
    Dnew[i,j+1] <- Dold[i, j]
    Xnew <- design_to_PWO(Dnew)
    Inew <- I_efficiency(Xnew, B)
    u = runif(1)
    t=0.5
    if(u<exp(-(Inew-Iold)*log(t+1))){
      Iold <- Inew
      Dold <- Dnew
      Xold <- Xnew
    }
    Ieffs[iter] <- Iold
  }
  plot(1:max_iter, Ieffs, type = "l", ylab = "I-efficiency", xlab = "Iteration")
  return(Dold)
}

  
# Implementation of the PWO I-optimality algorithm

m <- 6  # number of components
n <- 180 # target sample size


Df <- generate_full_design(m)
Xf <- design_to_PWO(Df)
Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design

set.seed(1234)
start_found <- FALSE
# we need to find a starting design with n rows that is not singular, just start
while(!start_found){
  D_start <- Df[sample(1:nrow(Df),n),] #assume that n permutations are possible?
  X_start <- design_to_PWO(D_start)
  M_start <- t(X_start)%*%X_start
  if(!is.singular.matrix(M_start)){
    start_found <- TRUE
    print("Found initial matrix")
  }
}


D_out <- random_exchange(D_start, Bf)
X_out <- design_to_PWO(D_out)
I_out <- I_efficiency(X_out, Bf)
print(I_full/I_out) # this is the relative efficiency # larger is better


# A non-singular matrix does not have inverse; There are linear dependent rows?
# When view the I-optimality criteria, trace(X_n^TX_n)^{-1}B
# Start with non-singular random starting design; until find one singular
# Try to find a non-singular so we can do the experiment
# 1. Start with D_0, start I=I(eff*D_0)
# 2. for t=0,1,2,...,n iter
        #2 choose a random row in D_t
        #3 choose a random row in D_{t+1} and switch the order of two adjacent components
        #4 Find I_{t+1}=I_{eff}(D_{t+1})
        #5 if (I_{t+1}<I_t)
        #      D_{best}=D_{t+1}
        #6 Return D_best



# n=1+2p p=choose(m,2) min choice of n is 1+p risky
# in poster m=4,5,6 n=1+2p, boxplots
# for each m, run the algorithm 30/100 times
# plot the relative I eff
# good if above 80%
# comparison of algorithm?
# random pairwise exchange
# delta

### test for multiple n when m=5

set.seed(1234)
tot_iter = 100

m = 5
nlist = seq(20,40,by=5)

I_list = matrix(NA,nrow = tot_iter, ncol = length(mlist)*length(nlist))

I_name = rep(0,length(mlist)*length(nlist))
for (i in 1:length(mlist)){
  for (j in 1:length(nlist)) {
    I_name[(i-1)*length(nlist)+j] <- c(paste("m=", m, ",n=", nlist[j], sep = ""))
  }
}
colnames(I_list) = I_name

for (n in 1:length(nlist)) {
  Df <- generate_full_design(m)
  Xf <- design_to_PWO(Df)
  Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
  I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
  for (iter in 1:tot_iter) {
    start_found <- FALSE
    # we need to find a starting design with n rows that is not singular
    while(!start_found){
      D_start <- Df[sample(1:nrow(Df),nlist[n]),]
      X_start <- design_to_PWO(D_start)
      M_start <- t(X_start)%*%X_start
      if(!is.singular.matrix(M_start)){
        start_found <- TRUE
        print("Found initial matrix")
      }
    }
    D_out <- random_exchange(D_start, Bf)
    X_out <- design_to_PWO(D_out)
    I_out <- I_efficiency(X_out, Bf)
    I_list[iter,(m-1)*length(nlist)+n]=I_full/I_out
  }
}

  


I_mean_mat=matrix(apply(I_list, 2, mean),ncol = length(nlist),byrow = FALSE)

rownames(I_mean_mat)=c(paste("m=", m, sep = ""))
colnames(I_mean_mat)=c(paste("n=", nlist, sep = ""))

preplot = c(I_list[,1],I_list[,2],I_list[,3])
preplot = as.data.frame(cbind(preplot,c(rep(25,100),rep(30,100),rep(35,100))))

colnames(preplot)=c("Ieff","n")
preplot$n=as.factor(preplot$n)

p <- ggplot(preplot, aes(x=n, y=Ieff)) + 
  geom_boxplot()+ggtitle("Box plot for different subset sizes") +
  ylab("I-efficiency")
p


### test for multiple n when m=6

set.seed(1234)
tot_iter = 10

m = 6
nlist = c(60,120,180,240)


I_list = matrix(NA,nrow = tot_iter, ncol = length(nlist))

I_name = rep(0,length(nlist))
for (j in 1:length(nlist)) {
  I_name[j] <- c(paste("m=", m, ",n=", nlist[j], sep = ""))
}
  

colnames(I_list) = I_name

for (n in 1:length(nlist)) {
  Df <- generate_full_design(m)
  Xf <- design_to_PWO(Df)
  Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
  I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
  for (iter in 1:tot_iter) {
    start_found <- FALSE
    # we need to find a starting design with n rows that is not singular
    while(!start_found){
      D_start <- Df[sample(1:nrow(Df),nlist[n]),]
      X_start <- design_to_PWO(D_start)
      M_start <- t(X_start)%*%X_start
      if(!is.singular.matrix(M_start)){
        start_found <- TRUE
        print("Found initial matrix")
      }
    }
    D_out <- random_exchange(D_start, Bf)
    X_out <- design_to_PWO(D_out)
    I_out <- I_efficiency(X_out, Bf)
    I_list[iter,n]=I_full/I_out
  }
}




I_mean_mat=matrix(apply(I_list, 2, mean),ncol = length(nlist),byrow = FALSE)

rownames(I_mean_mat)=c(paste("m=", m, sep = ""))
colnames(I_mean_mat)=c(paste("n=", nlist, sep = ""))

preplot = matrix(I_list,ncol=1)
preplot = as.data.frame(cbind(preplot,rep(nlist,each=tot_iter)))

colnames(preplot)=c("Ieff","n")
preplot$n=as.factor(preplot$n)

p <- ggplot(preplot, aes(x=n, y=Ieff)) + 
  geom_boxplot()+ggtitle("Box plot for different subset sizes when m=6") +
  ylab("Relative I-efficiency")
p

#--------------------------------------------
### Comparison between PWO and Annealing
#--------------------------------------------
set.seed(1234)
tot_iter = 10

m = 5
nlist = c(25,30,35)


I_list_PWO = matrix(NA,nrow = tot_iter, ncol = length(nlist))
I_list_annealing = matrix(NA,nrow = tot_iter, ncol = length(nlist))

I_name = rep(0,length(nlist))
for (j in 1:length(nlist)) {
  I_name[j] <- c(paste("m=", m, ",n=", nlist[j], sep = ""))
}


colnames(I_list) = I_name

for (n in 1:length(nlist)) {
  Df <- generate_full_design(m)
  Xf <- design_to_PWO(Df)
  Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
  I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
  for (iter in 1:tot_iter) {
    start_found <- FALSE
    # we need to find a starting design with n rows that is not singular
    while(!start_found){
      D_start <- Df[sample(1:nrow(Df),nlist[n]),]
      X_start <- design_to_PWO(D_start)
      M_start <- t(X_start)%*%X_start
      if(!is.singular.matrix(M_start)){
        start_found <- TRUE
        print("Found initial matrix")
      }
    }
    D_out_PWO <- random_exchange(D_start, Bf)
    X_out_PWO <- design_to_PWO(D_out_PWO)
    I_out_PWO <- I_efficiency(X_out_PWO, Bf)
    I_list_PWO[iter,n]=I_full/I_out_PWO
    D_out_annealing <- annealing_exchange(D_start, Bf)
    X_out_annealing <- design_to_PWO(D_out_annealing)
    I_out_annealing <- I_efficiency(X_out_annealing, Bf)
    I_list_annealing[iter,n]=I_full/I_out_annealing
  }
}




I_mean_mat=matrix(apply(I_list, 2, mean),ncol = length(nlist),byrow = FALSE)

rownames(I_mean_mat)=c(paste("m=", m, sep = ""))
colnames(I_mean_mat)=c(paste("n=", nlist, sep = ""))

preplot = cbind(matrix(I_list_PWO,ncol = 1),matrix(I_list_annealing,ncol = 1))

preplot = as.data.frame(cbind(preplot,rep(nlist,each=tot_iter)))

colnames(preplot)=c("RandomExchange","AnnealingSimulation","n")
preplot$n=as.factor(preplot$n)

p1 <- ggplot(preplot, aes(x=n, y=Ieff1)) + 
  geom_boxplot()+ggtitle("Box plot for different subset sizes") +
  ylab("I-efficiency")
p1
p2 <- ggplot(preplot, aes(x=n, y=Ieff2)) + 
  geom_boxplot()+ggtitle("Box plot for different subset sizes") +
  ylab("I-efficiency")
p2
preplot_long = melt(preplot,id.vars = c("n"))
p <- ggplot(preplot_long, aes(x=n,y=value,fill=n))+
  geom_boxplot() + labs(title="Algorithm Efficiency Comparison When m=5")+ facet_wrap(~variable) +  ylab("Relative I-efficiency")
p




