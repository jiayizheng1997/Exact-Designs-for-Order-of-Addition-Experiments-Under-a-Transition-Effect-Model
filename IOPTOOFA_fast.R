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
  for(row in 1:n){
    index <- 1
    temp <- rep(0, q)
    for(i in 1:(m-1)){
      for(j in (i+1):m){
        k <- which(D[row,] == i)
        l <- which(D[row,] == j)
        temp[index] <- sign(l-k)
        index <- index + 1
      }
    }
    X[row,] <- temp
  }
  # cbind(1, X) adds a column of 1s for the intercept
  X <- cbind(1, X)
  colnames(X) <- cnames
  return(X)
}

design_to_PWO_fast <- function(D){
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
        X[, index] <- sign(k[,i]-k[,j])
        index <- index + 1
      }
  }
  # cbind(1, X) adds a column of 1s for the intercept
  X <- cbind(1, X)
  colnames(X) <- cnames
  return(X)
}

D = generate_full_design(8)
D = generate_full_design(10)

t1=Sys.time()
Z=design_to_PWO(D)
t2=Sys.time()
t2-t1

t1=Sys.time()
Z=design_to_PWO_fast(D)
t2=Sys.time()
t2-t1



trace <- function(X){
  return(sum(diag(X)))
}

I_efficiency <- function(X, B){
  
  M <- t(X)%*%X/nrow(X)
  if(is.singular.matrix(M)){
    return(1e10)
  }
  
  return(trace( solve(M)%*%B))
  
}

# random_exchange
# randomly selects adjacent pairs of components to swap
# performs the swap if it improves I optimality
# TODO: Try to improve this. Move of the iterations are being wasted..
# inputs: 
# D: a n by m matrix of permutations of (1:m)
# B: the B matrix for the I opt criterion
# max_iter: max number of iterations, default is 1e4
random_exchange <- function(D, B, max_iter = 1e4){
  Dold <- D
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- design_to_PWO(Dold)
  Cold <- solve(t(Xold)%*%Xold)
  Iold <- trace(Cold%*%B)
  # Iold <- I_efficiency(Xold, B)
  delta <- Inf
  Ieffs <- numeric(max_iter)
  for(iter in 1:max_iter){
    
    print(paste("Iteration", iter, sep = " "))
    
    
    i <- sample(1:n, 1)  # random row
    j<- sample(1:(m-1),1) # random column
    
    Dnew <- Dold
    # exchange component j with component (j+1)
    Dnew[i,j] <- Dold[i,j+1]
    Dnew[i,j+1] <- Dold[i, j]
    XY <- design_to_PWO(rbind(Dold[i,], Dnew[i,]))
    x <- XY[1,]
    y <- XY[2,]
    F1 <- cbind(y,-x)
    F2 <- cbind(y,x)
    I2 <- diag(2)
    Cnew <- Cold - Cold%*%F1%*%solve(I2 + t(F2)%*%Cold%*%F1)%*%t(F2)%*%Cold
    Inew <- trace(Cnew%*%B)
    # Xnew <- design_to_PWO(Dnew)
    # Inew <- I_efficiency(Xnew, B)
    if(Inew < Iold){
      Iold <- Inew
      Dold <- Dnew
      Cold <- Cnew
    }
    Ieffs[iter] <- Iold
    
    
    
    
  }
  
  plot(1:max_iter, Ieffs, type = "l", ylab = "I-efficiency", xlab = "Iteration")
  
  return(Dold)
  
  
}

### simulated annealing
simulated_annealing <- function(D, B, max_iter = 1e4){
  Dold <- D
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- design_to_PWO(Dold)
  Cold <- solve(t(Xold)%*%Xold)
  Iold <- trace(Cold%*%B)
  # Iold <- I_efficiency(Xold, B)
  delta <- Inf
  Ieffs <- numeric(max_iter)
  for(iter in 1:max_iter){
    
    print(paste("Iteration", iter, sep = " "))
    
    
    i <- sample(1:n, 1)  # random row
    j<- sample(1:(m-1),1) # random column
    
    Dnew <- Dold
    # exchange component j with component (j+1)
    Dnew[i,j] <- Dold[i,j+1]
    Dnew[i,j+1] <- Dold[i, j]
    XY <- design_to_PWO(rbind(Dold[i,], Dnew[i,]))
    x <- XY[1,]
    y <- XY[2,]
    F1 <- cbind(y,-x)
    F2 <- cbind(y,x)
    I2 <- diag(2)
    Cnew <- Cold - Cold%*%F1%*%solve(I2 + t(F2)%*%Cold%*%F1)%*%t(F2)%*%Cold
    Inew <- trace(Cnew%*%B)
    u = runif(1)
    if(u<exp(-(Inew-Iold)*log(iter+1))){
      Iold <- Inew
      Dold <- Dnew
      Cold <- Cnew
    }
    Ieffs[iter] <- Iold
  }
  plot(1:max_iter, Ieffs, type = "l", ylab = "I-efficiency", xlab = "Iteration")
  return(Dold)
}


###--------------------------------
### Test of functions when m=8,n=60
###--------------------------------

m <- 8  # number of components
n <- 60 # target sample size


Df <- generate_full_design(m)
Xf <- design_to_PWO(Df)
Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design

set.seed(1234)
start_found <- FALSE
# we need to find a starting design with n rows that is not singular
while(!start_found){
  D_start <- Df[sample(1:nrow(Df),n),]
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
print(I_full/I_out) # this is the relative efficiency



### test for multiple n when m is any positive integer, default to be 5.

set.seed(1234)

# The value of tot_iter, m, nlist needs to be input by hand
tot_iter = 100
m = 5
nlist = c(30,60,90)

###################################################################
###A function to test for m and multiple n for random exchange###
###################################################################

PWO_test = function(tot_iter = 100,m = 5,nlist = c(30,60,90)){
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
    geom_boxplot()+ggtitle(paste("Box plot for different subset sizes when m =",m)) +
    ylab("Relative I-efficiency")
  return(p)
  
}

PWO_test(100,5,c(30,60,90))

test_10 = PWO_test(100,10,c(30,60,90))


#####################################################################
###A function to test for m and multiple n for simulated annealing###
#####################################################################

sa_test = function(tot_iter = 100,m = 5,nlist = c(30,60,90)){
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
      D_out <- simulated_annealing(D_start, Bf)
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
    geom_boxplot()+ggtitle(paste("Box plot for different subset sizes when m =",m)) +
    ylab("Relative I-efficiency")
  return(p)
  
}

sa_test(100,5,c(30,60,90))

test_10 = sa_test(100,10,c(30,60,90))


#--------------------------------------------
### Comparison between PWO and Simulated Annealing
#--------------------------------------------
set.seed(1234)
tot_iter = 100

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
    D_out_annealing <- simulated_annealing(D_start, Bf)
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
  geom_boxplot() + labs(title=paste("Algorithm Efficiency Comparison When m = ",m))+ facet_wrap(~variable) +  ylab("Relative I-efficiency")
p
