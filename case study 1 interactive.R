library(combinat)
library(matrixcalc)

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

#New random_exchange:
#We first randomly select elements, then we see if the swap is legal or not.
#We would randomly select two elements from the same random row in block matrix B

# random_exchange
# randomly selects adjacent pairs of components to swap
# performs the swap if it improves I optimality
# TODO: Try to improve this. Move of the iterations are being wasted..
# inputs: 
# D: a n by m matrix of permutations of (1:m)
# B: the B matrix for the I opt criterion

# max_iter: max number of iterations, default is 1e4
random_exchange <- function(D, B, int1, int2, max_iter = 1e4){
  Dold <- D
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- design_to_PWO(Dold)
  Xold <- cbind(Xold,Xold[,int1]*Xold[,int2])
  Cold <- solve(t(Xold)%*%Xold/nrow(Xold))
  Iold <- trace(Cold%*%B)
  Ieffs <- numeric(max_iter)
  for(iter in 1:max_iter){
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
    Xnew <- design_to_PWO(Dnew)
    Xnew <- cbind(Xnew,Xnew[,int1]*Xnew[,int2])
    Inew <- I_efficiency(Xnew, B)
    if(Inew < Iold){
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

### simulated annealing
simulated_annealing <- function(D, B, int1, int2,  max_iter = 1e4){
  Dold <- D
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
    Xnew <- design_to_PWO(Dnew)
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

################################################################################
####################### Real Data Starts from Here #############################
################################################################################

set.seed(1234)
# read in data
case1 <- read.csv("case1.csv",header = FALSE,col.names = c('X1','X2','X3','X4','Y'))
case1[,1:4] = case1[,1:4]+1
Df = as.matrix(case1[,1:4])

Xf = design_to_PWO(Df)
data1 = as.data.frame(cbind(as.data.frame(case1)$Y,Xf))
colnames(data1)[1] = 'Y'

# Here we are detecting which interaction is strong

t = as.data.frame(matrix(0,nrow=2,ncol=(choose(6,2))))

# fit a full model without constraints
full_model = lm(Y ~ 0+., data = data1)
summary(full_model)
index = 2
for (i in 2:(ncol(Xf)-1)) {
  for (j in (i+1):ncol(Xf)) {
    temp_X = cbind(Xf,Xf[,i]*Xf[,j])
    temp_data = as.data.frame(cbind(data1,Xf[,i]*Xf[,j]))
    int_model = lm(Y ~ 0+., data = temp_data)
    sm = summary(int_model)
    colnames(t)[index-1] = paste(colnames(Xf)[i],colnames(Xf)[j],sep='')
    t[1,index-1] = sm$coefficients[8,3]
    t[2,index-1] = pt(abs(t[1,index-1]),df=24-7-1,lower.tail = FALSE)
    index = index+1
  }
}
t

# from the result we notice that z13z23 is the strongest interaction

# Here is a block of codes studying the I-efficiency with interaction

Y = case1$Y
rank = sort(Y,decreasing = FALSE)
Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
n = ncol(Xf)+2

Int_Ieff <- matrix(0,nrow=3,ncol=choose(6,2)+1)
cnames <- rep(0,ncol(Int_Ieff))
Int_Ieff[1,1] <- I_full
cnames[1] <- 'No interaction'
index = 2
for (i in 2:(ncol(Xf)-1)) {
  for (j in (i+1):ncol(Xf)) {
    temp_X = cbind(Xf,Xf[,i]*Xf[,j])
    temp_B = t(temp_X)%*%temp_X
    Int_Ieff[1,index] = I_efficiency(temp_X,temp_B)
    cnames[index] = paste(colnames(Xf)[i],colnames(Xf)[j]) 
    iter = 100
    RE_eff = rep(0,iter)
    SA_eff = rep(0,iter)
    for (k in 1:iter) {
      start_found <- FALSE
      while(!start_found){
        D_start <- Df[sample(1:nrow(Df),n),]
        X_start <- design_to_PWO(D_start)
        X_start <- cbind(X_start,X_start[,i]*X_start[,j])
        #X_start <- X_start[,-filter$cons_cols]
        M_start <- t(X_start)%*%X_start
        if(!is.singular.matrix(M_start)){
          start_found <- TRUE
          print("Found initial matrix")
        }
        skip_to_next <- FALSE
        tryCatch(RE_eff[k] <- random_exchange(D_start, temp_B, int1=i, int2=j, max_iter = 100)[[2]], error=function(e){skip_to_next<<-TRUE})
        if(skip_to_next) { next }
        skip_to_next <- FALSE
        tryCatch(SA_eff[k] <- simulated_annealing(D_start, temp_B, int1=i, int2=j, max_iter = 100)[[2]], error=function(e){skip_to_next<<-TRUE})
        if(skip_to_next) { next }
      }
      
    }
    Int_Ieff[2,index] = mean(RE_eff)
    Int_Ieff[3,index] = mean(SA_eff)
    index = index+1
  }
}
colnames(Int_Ieff) = cnames
rownames(Int_Ieff) = c('Full','RE','SA')

# now we start calculating the rank with z13z23 interaction added

n = 9
int_X = cbind(Xf,Xf[,3]*Xf[,5])
Bf = t(int_X)%*%int_X
case1_int = as.data.frame(cbind(Y,int_X))
colnames(case1_int)[9] = 'z13z23'
n_sim <- 100
rank_RE <- rep(NA,n_sim)
rank_SA <- rep(NA,n_sim)
rank_rand <- rep(NA,n_sim)
MSE_list <-matrix(NA,nrow = n_sim,ncol = 3)
iter = 1
while (iter<=n_sim){
  start_found <- FALSE
  # we need to find a starting design with n rows that is not singular
  while(!start_found){
    D_start <- Df[sample(1:nrow(Df),n),]
    X_start <- design_to_PWO(D_start)
    X_start <- cbind(X_start,X_start[,2]*X_start[,5])
    #X_start <- X_start[,-filter$cons_cols]
    M_start <- t(X_start)%*%X_start
    if(!is.singular.matrix(M_start)){
      start_found <- TRUE
      print("Found initial matrix")
    }
  }
  #errors are likely to occur, hence we apply a trycatch
  skip_to_next <- FALSE
  tryCatch(D_out_RE<-random_exchange(D_start,Bf,int1=3,int2=5,max_iter=100),error=function(e){skip_to_next<<-TRUE})
  if(skip_to_next) { next }   
  X_out_RE <- design_to_PWO(D_out_RE[[1]])
  X_out_RE <- cbind(X_out_RE,X_out_RE[,3]*X_out_RE[,5])
  skip_to_next <- FALSE
  tryCatch(D_out_SA <- simulated_annealing(D_start,Bf,int1=3,int2=5,max_iter = 100),error=function(e){skip_to_next<<-TRUE})
  if(skip_to_next) { next } 
  X_out_SA <- design_to_PWO(D_out_SA[[1]])
  X_out_SA <- cbind(X_out_SA,X_out_SA[,3]*X_out_SA[,5]) 
  # generate the regression models and calculate ranks of obtained optimal designs
  
  # random exchange
  Y_out_RE <- rep(0,n)
  # separate regressed and predicted data
  for (i in 1:n) {
    Y_out_RE[i] = Y[which(apply(Df, 1, function(x) identical(D_out_RE[[1]][i,],x)))]
  }
  X_out_RE = X_out_RE[,-1]
  case1_RE = as.data.frame(cbind(X_out_RE,Y_out_RE))
  colnames(case1_RE)[7] = 'z13z23'
  cons_RE_model = lm(Y_out_RE ~ ., data = case1_RE)
  RE_pred = predict(cons_RE_model,newdata = case1_int)
  rank_RE[iter] = which(rank == min(Y[which(RE_pred==min(RE_pred))]))
  
  # simulated annealing
  Y_out_SA <- rep(0,n)
  # separate regressed and predicted data
  for (i in 1:n) {
    Y_out_SA[i] = Y[which(apply(Df, 1, function(x) identical(D_out_SA[[1]][i,],x)))]
  }
  X_out_SA = X_out_SA[,-1]
  case1_SA = as.data.frame(cbind(X_out_SA,Y_out_SA))
  colnames(case1_SA)[7] = 'z13z23'
  cons_SA_model = lm(Y_out_SA ~ ., data = case1_SA)
  SA_pred = predict(cons_SA_model,newdata = case1_int)
  rank_SA[iter] = which(rank == min(Y[which(SA_pred==min(SA_pred))]))
  
  #complete random
  row_out_rand = sample(nrow(Xf),n)
  Y_out_rand = Y[row_out_rand]
  X_out_rand = cbind(Xf,Xf[,3]*Xf[,5])
  X_out_rand = X_out_rand[row_out_rand,-1]
  case1_rand = as.data.frame(cbind(X_out_rand,Y_out_rand))
  colnames(case1_rand)[7] = 'z13z23'
  cons_rand_model = lm(Y_out_rand ~ ., data = case1_rand)
  rand_pred = predict(cons_rand_model,newdata = case1_int)
  rank_rand[iter] = which(rank == min(Y[which(rand_pred==min(rand_pred))]))
  MSE = c(mean(sum((Y-RE_pred)^2)),mean(sum((Y-SA_pred)^2)),mean(sum(Y-rand_pred)^2))
  MSE_list[iter,] = MSE
  iter = iter +1
  
}    





result = matrix(c(length(which(rank_rand<=3)),
                  length(which(rank_RE<=3)),
                  length(which(rank_SA<=3)),
                  length(which(rank_rand==1)),
                  length(which(rank_RE==1)),
                  length(which(rank_SA==1))),nrow = 2,byrow = TRUE)
rownames(result) = c("leq 3","equal 1")
colnames(result) = c("rand","RE","SA")
result
data_result = list(rank_rand,rank_RE,rank_SA,result)
hist(data_result[[1]])
hist(data_result[[2]])
hist(data_result[[3]])
length(which(data_result[[1]]<=3))
length(which(data_result[[2]]<=3))
length(which(data_result[[3]]<=3))
length(which(data_result[[1]]==1))
length(which(data_result[[2]]==1))
length(which(data_result[[3]]==1))
colMeans(MSE_list)
hist(MSE_list[,1])
hist(MSE_list[,2])
hist(MSE_list[,3])

#SSE
#summary(full_model)
summary(cons_RE_model)
summary(cons_SA_model)
summary(cons_rand_model)
summary(MSE_list)


# Here is the result without the interaction term
# To begin with, we need to rewrite the algorithms to drop the interaction parameter

random_exchange <- function(D, B, max_iter = 1e4){
  Dold <- D
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- design_to_PWO(Dold)
  Cold <- solve(t(Xold)%*%Xold/nrow(Xold))
  Iold <- trace(Cold%*%B)
  Ieffs <- numeric(max_iter)
  for(iter in 1:max_iter){
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
    Xnew <- design_to_PWO(Dnew)
    Inew <- I_efficiency(Xnew, B)
    if(Inew < Iold){
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

### simulated annealing
simulated_annealing <- function(D, B, max_iter = 1e4){
  Dold <- D
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- design_to_PWO(Dold)
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
    Xnew <- design_to_PWO(Dnew)
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

n = 9
int_X = Xf
Bf = t(int_X)%*%int_X
case1_int = as.data.frame(cbind(Y,int_X))
# colnames(case1_int)[9] = 'z12z23'
n_sim <- 100
rank_RE <- rep(NA,n_sim)
rank_SA <- rep(NA,n_sim)
rank_rand <- rep(NA,n_sim)
MSE_list <-matrix(NA,nrow = n_sim,ncol = 3)
iter = 1
while (iter<=n_sim){
  start_found <- FALSE
  # we need to find a starting design with n rows that is not singular
  while(!start_found){
    D_start <- Df[sample(1:nrow(Df),n),]
    X_start <- design_to_PWO(D_start)
    #X_start <- cbind(X_start,X_start[,2]*X_start[,5])
    #X_start <- X_start[,-filter$cons_cols]
    M_start <- t(X_start)%*%X_start
    if(!is.singular.matrix(M_start)){
      start_found <- TRUE
      print("Found initial matrix")
    }
  }
  #errors are likely to occur, hence we apply a trycatch
  skip_to_next <- FALSE
  tryCatch(D_out_RE<-random_exchange(D_start,Bf,max_iter=100),error=function(e){skip_to_next<<-TRUE})
  if(skip_to_next) { next }   
  X_out_RE <- design_to_PWO(D_out_RE[[1]])
  #X_out_RE <- cbind(X_out_RE,X_out_RE[,2]*X_out_RE[,5])
  skip_to_next <- FALSE
  tryCatch(D_out_SA <- simulated_annealing(D_start,Bf,max_iter = 100),error=function(e){skip_to_next<<-TRUE})
  if(skip_to_next) { next } 
  X_out_SA <- design_to_PWO(D_out_SA[[1]])
  #X_out_SA <- cbind(X_out_SA,X_out_SA[,2]*X_out_SA[,5]) 
  # generate the regression models and calculate ranks of obtained optimal designs
  
  # random exchange
  Y_out_RE <- rep(0,n)
  # separate regressed and predicted data
  for (i in 1:n) {
    Y_out_RE[i] = Y[which(apply(Df, 1, function(x) identical(D_out_RE[[1]][i,],x)))]
  }
  X_out_RE = X_out_RE[,-1]
  case1_RE = as.data.frame(cbind(X_out_RE,Y_out_RE))
  #colnames(case1_RE)[7] = 'z12z23'
  cons_RE_model = lm(Y_out_RE ~ ., data = case1_RE)
  RE_pred = predict(cons_RE_model,newdata = case1_int)
  rank_RE[iter] = which(rank == min(Y[which(RE_pred==min(RE_pred))]))
  
  # simulated annealing
  Y_out_SA <- rep(0,n)
  # separate regressed and predicted data
  for (i in 1:n) {
    Y_out_SA[i] = Y[which(apply(Df, 1, function(x) identical(D_out_SA[[1]][i,],x)))]
  }
  X_out_SA = X_out_SA[,-1]
  case1_SA = as.data.frame(cbind(X_out_SA,Y_out_SA))
  #colnames(case1_SA)[7] = 'z12z23'
  cons_SA_model = lm(Y_out_SA ~ ., data = case1_SA)
  SA_pred = predict(cons_SA_model,newdata = case1_int)
  rank_SA[iter] = which(rank == min(Y[which(SA_pred==min(SA_pred))]))
  
  #complete random
  row_out_rand = sample(nrow(Xf),n)
  Y_out_rand = Y[row_out_rand]
  #X_out_rand = cbind(Xf,Xf[,2]*Xf[,5])
  X_out_rand = Xf[row_out_rand,-1]
  case1_rand = as.data.frame(cbind(X_out_rand,Y_out_rand))
  #colnames(case1_rand)[7] = 'z12z23'
  cons_rand_model = lm(Y_out_rand ~ ., data = case1_rand)
  rand_pred = predict(cons_rand_model,newdata = case1_int)
  rank_rand[iter] = which(rank == min(Y[which(rand_pred==min(rand_pred))]))
  MSE = c(mean(sum((Y-RE_pred)^2)),mean(sum((Y-SA_pred)^2)),mean(sum(Y-rand_pred)^2))
  MSE_list[iter,] = MSE
  iter = iter +1
  
}    





result0 = matrix(c(length(which(rank_rand<=3)),
                  length(which(rank_RE<=3)),
                  length(which(rank_SA<=3)),
                  length(which(rank_rand==1)),
                  length(which(rank_RE==1)),
                  length(which(rank_SA==1))),nrow = 2,byrow = TRUE)
rownames(result0) = c("leq 3","equal 1")
colnames(result0) = c("rand","RE","SA")
result0
data_result0 = list(rank_rand,rank_RE,rank_SA,result0)
hist(data_result0[[1]])
hist(data_result0[[2]])
hist(data_result0[[3]])
length(which(data_result0[[1]]<=3))
length(which(data_result0[[2]]<=3))
length(which(data_result0[[3]]<=3))
length(which(data_result0[[1]]==1))
length(which(data_result0[[2]]==1))
length(which(data_result0[[3]]==1))
colMeans(MSE_list)
hist(MSE_list[,1])
hist(MSE_list[,2])
hist(MSE_list[,3])

#SSE
#summary(full_model)
summary(cons_RE_model)
summary(cons_SA_model)
summary(cons_rand_model)
summary(MSE_list)