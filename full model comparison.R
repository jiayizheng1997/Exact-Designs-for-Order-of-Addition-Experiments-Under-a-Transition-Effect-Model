library(combinat)
library(matrixcalc)
library(stats)

trace <- function(X){
  return(sum(diag(X)))
}

I_efficiency <- function(X, B){
  
  M <- t(X)%*%X/nrow(X)
  if(rcond(M) < 1e-10){
    return(NA)
  }
  
  return(trace( solve(M)%*%B))
  
}

# Takes a model matrix X
# returns D_efficiency
D_efficiency <- function(X){
  k <- ncol(X)
  N <- nrow(X)
  return( (1/N)*det(t(X)%*%X)^(1/k) )
}

design_to_PWO <- function(D){
  m <- ncol(D)
  n <- nrow(D)
  q <- choose(m, 2)
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
  k <- matrix(0, nrow = n, ncol = m)
  for (i in 1:m) {
    k[,i] <- apply(D,1,function(x) which(x==i))
  }
  index <- 1
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      X[, index] <- sign(k[,j]-k[,i])
      index <- index + 1
    }
  }
  # cbind(1, X) adds a column of 1s for the intercept
  X <- cbind(1, X)
  colnames(X) <- cnames
  return(X)
}

# OofA_transition_model_matrix
# Input: D, an n by m matrix of permutations
# Output: X, an n by 2*(m choose 2) matrix of model parameters
OofA_transition_model_matrix <- function(D){
  n <- nrow(D)
  m <- ncol(D)
  X2 <- matrix(0, nrow = n, ncol = m*(m-1))
  x2namevec <- c()
  counter <- 1
  for(i in 1:m){
    for(j in 1:m){
      if(i != j){
        x2namevec <- c(x2namevec, paste("x", i, "|", j, sep = ""))
      }
    }
  }
  
  for(i in 1:n){
    
    for(j in 2:m){
      
      prev <- D[i,j-1]
      curr <- D[i,j]
      col_ind <- which(x2namevec == paste("x", curr, "|", prev, sep = ""))
      X2[i,col_ind] <- 1
      
    }
    
  }
  
  colnames(X2) <- x2namevec
  X <- X2
  return(X)
}

## get baseline CP matrix
# generates a matrix of the CP variables (assuming baseline coding)
# input: a matrix D of permutations of (1,2,...,m)
# z_cj = 1 if component c is in position j, 0 otherwise
design_to_CP <- function(D){
  
  if(is.vector(D)){
    m = length(D)
  } else{
    m = ncol(D)
  }
  
  
  CPnames <- c("Intercept")
  for(j in 1:(m-1)){
    for(c in 2:m){
      CPnames <- c(CPnames, paste("z",c,j,sep=""))
    }
  }
  
  if(is.vector(D) || nrow(D) == 1){
    q <- length(CPnames)
    X <- numeric(q)
    col <- 2
    for(j in 1:(m-1)){
      for(c in 2:m){
        if(D[j] == c){
          X[col] <- 1
        }
        col <- col + 1
      }
    }
    names(X) <- CPnames
    return(X)
    
  }
  
  
  m <- ncol(D)
  n <- nrow(D)
  q <- 1+(m-1)^2
  X <- matrix(0,nrow=n,ncol=q)
  X[,1] <- 1
  
  for(row in 1:n){
    
    col <- 2 
    for(j in 1:(m-1)){ 
      
      for(c in 2:m){ 
        
        if(D[row,j] == c){
          X[row,col] <- 1
        }
        col<- col + 1
        
        
      }
      
    }
    
    
  }
  
  colnames(X) <- CPnames
  
  
  return(X)
  
}

### generates the full design with all m! permutations
generate_full_design <- function(m){
  Plist <- permn(1:m) ## get all permutations of 1:m
  return(matrix(unlist(Plist),byrow=TRUE,nrow=factorial(m)))
}


### simulated annealing
### simulated annealing
simulated_annealing <- function(n, D, B, int1, int2,  max_iter = 100){
  
  start_found <- FALSE
  # we need to find a starting design with vote_n rows that is not singular
  while(!start_found){
    D_start <- Df[sample(1:nrow(Df),n),]
    X_start <- OofA_transition_model_matrix(D_start)
    X_start <- cbind(X_start,X_start[,int1]*X_start[,int2])
    M_start <- t(X_start)%*%X_start
    print("Looking for initial matrix")
    if(!is.singular.matrix(M_start)){
      start_found <- TRUE
      print("Found initial matrix")
    }
  }
  Dold <- D_start
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- design_to_PWO(Dold)
  Xold <- cbind(Xold,Xold[,int1]*Xold[,int2])
  Cold <- solve(t(Xold)%*%Xold/nrow(Xold))
  Iold <- trace(Cold%*%B)
  Ieffs <- numeric(max_iter)
  for (iter in 1:max_iter) {
    print(paste("Iteration", iter, sep = " "))
    pass <- FALSE
    while (!pass) {
      i <- sample(1:n, 1)  # random row
      j<- sample(1:(m-1),1) # random column
      pass<-TRUE
    }
    Dnew <- Dold
    # exchange component j with component (j+1)
    Dnew[i,j] <- Dold[i,j+1]
    Dnew[i,j+1] <- Dold[i, j]
    Xnew <- OofA_transition_model_matrix(Dnew)
    Xnew <- cbind(Xnew,Xnew[,int1]*Xnew[,int2])
    Inew <- I_efficiency(Xnew, B)
    u = runif(1)
    if(u<exp(-(Inew-Iold)*iter)){
      Iold <- Inew
      Dold <- Dnew
    }
    Ieffs[iter] <- Iold
  }
  plot(1:max_iter, Ieffs, type = "l", ylab = "I-efficiency", xlab = "Iteration")
  Ieffs = Ieffs[iter]
  result = list(Dold,Ieffs)
  return(result)
}

# generate_random_design
# m: number of components
# n: sample size
# output: an n times m design matrix where each row is a permutation of 1:m
generate_random_design <- function(m, n){
  
  D0 <- matrix(NA, nrow = n, ncol = m)
  D0 <- t(apply(D0,1,function(x){sample(1:m)}))
  
  
  return(D0)
  
}

# read in data
case1 <- read.csv("case1.csv",header = FALSE,col.names = c('X1','X2','X3','X4','Y'))
case1[,1:4] = case1[,1:4]+1
Df = as.matrix(case1[,1:4])
int1 = 0
int2 = 0
Y <- case1$Y

# TE model

Xf <- OofA_transition_model_matrix(Df)
rank = sort(Y,decreasing = FALSE)
# Xf <- cbind(Xf,Xf[,int1]*Xf[,int2])
Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
data1 = as.data.frame(cbind(Y,Xf))
TE_model <- lm(Y ~ -1 + .,data = data1)

# PWO model

Xf = design_to_PWO(Df)
data1 = as.data.frame(cbind(Y,Xf))
PWO_model = lm(Y ~ 0+., data = data1)

# Triplet Model

Xf = design_to_PWO(Df)
Xf = Xf[,2:7]
Xf = cbind(Xf,Xf[,1]*Xf[,2],Xf[,1]*Xf[,4],Xf[,1]*Xf[,3],Xf[,1]*Xf[,5],Xf[,2]*Xf[,3],Xf[,2]*Xf[,6],Xf[,4]*Xf[,5],Xf[,4]*Xf[,6])
Xf = cbind(Xf,Xf[,1]*Xf[,3]*Xf[,4],Xf[,1]*Xf[,2]*Xf[,5],Xf[,1]*Xf[,2]*Xf[,6],Xf[,1]*Xf[,3]*Xf[,6],Xf[,1]*Xf[,4]*Xf[,6],Xf[,1]*Xf[,5]*Xf[,6],Xf[,2]*Xf[,3]*Xf[,4],Xf[,2]*Xf[,3]*Xf[,5],Xf[,3]*Xf[,5]*Xf[,6])
data1 = as.data.frame(cbind(Y,Xf))
Triplet_model = lm(Y ~ 0+., data = data1)

# CP Model
Xf = design_to_CP(Df)
data1 = as.data.frame(cbind(Y,Xf))
CP_model = lm(Y ~ 0+., data = data1)

AIC(TE_model,PWO_model,Triplet_model,CP_model)
BIC(TE_model,PWO_model,Triplet_model,CP_model)

MSE = c(mean(TE_model$residuals^2),mean(PWO_model$residuals^2),mean(Triplet_model$residuals^2),mean(CP_model$residuals^2))
MSE










