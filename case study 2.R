library(combinat)
library(matrixcalc)

#m: number of components
#b: number of blocks
#e: extremeness, it takes value 0 or 1
#If e is 1 then b can only be 2.
#mu: number of free components
cons_block <- function(m,b,e,mu){  
  temp <- sample(m,m-mu)
  if (e==0){ #when extremeness is zero, we want to distribute the constraints equally
    n1=(m-mu)%/%b #in each block, there would be at least n1 components
    n2=(m-mu)%%b #in n2 blocks, there would be n1+1 components
    #we want to write a row of length m as a matrix of n1*(n1+1)
    #Fill the matrix of repeated values with zero then filter the zero
    cons_block = matrix(temp,nrow = b,ncol = n1+((m-mu)%%b>0))
    if (n2>0){for (i in 1:(b-n2)){cons_block[n2+i,n1+1]=0}}
  }
  else if(e==1){
    if (b!=2){
      print("cannot perform extreme divide")
      return()
    }
    else {
      b1 = rep(0,length(temp)-1)
      b1[1] = sample(temp,1)
      b2 = temp[which(temp!=b1[1])]
      cons_block = rbind(b2,b1)
    }
  }
  return(cons_block)
}
#B should be a matrix whose row number are constraint blocks
cons_pair <- function(B){ #input the matrix obtained above 
  cons_pair=matrix(c(0,0),ncol=2)
  #get all binary combinations between two rows
  for (i in 1:(nrow(B)-1)) {
    for (j in 1:(nrow(B)-i)){
      temp=lapply(B[i,],function(X){lapply(B[i+j,],function(Y){c(X,Y)})})
      temp=matrix(unlist(temp),ncol=2,byrow = TRUE)
      cons_pair=rbind(cons_pair,temp)
    }
  }
  cons_pair=cons_pair[apply(cons_pair,1,function(X) all(X!=0)),]
  return(cons_pair) #the final output is a matrix with 2 columns
}

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

#This function is used to filter the Z matrix based on the constraints
#cons is the 2-column constraint matrix
#Z is the full PWO design matrix
PWO_filter <- function(cons,Z){ 
  cons_rows <- c()
  cons_cols <- c()
  if (is.null(nrow(cons))){
    cons <- matrix(cons,nrow = 1)
  }
  for (i in 1:nrow(cons)) {
    if (cons[i,1]<cons[i,2]){
      temp_col=paste('z',cons[i,1],cons[i,2],sep="")
      cons_rows <- c(cons_rows,which(Z[,temp_col]==-1))
      cons_cols <- c(cons_cols,which(colnames(Z)==temp_col))
    }
    else{
      temp_col=paste('z',cons[i,2],cons[i,1],sep="")
      cons_rows <- c(cons_rows,which(Z[,temp_col]==1))
      cons_cols <- c(cons_cols,which(colnames(Z)==temp_col))
    }
    cons_rows=unique(cons_rows)
  }
  output = list(cons_rows=cons_rows,cons_cols=cons_cols)
  return(output)
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
# block: the blockwise constraint matrix
# max_iter: max number of iterations, default is 1e4
random_exchange <- function(D, B, block, max_iter = 1e4){
  Dold <- D
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- design_to_PWO(Dold)
  cons <- cons_pair(block)
  filter <- PWO_filter(cons=cons, Z=Xold)
  Xold <- Xold[,-filter$cons_cols]  #Calculate the constrained I-efficiency
  Cold <- solve(t(Xold)%*%Xold/nrow(Xold))
  Iold <- trace(Cold%*%B)
  # Iold <- I_efficiency(Xold, B)
  Ieffs <- numeric(max_iter)
  for(iter in 1:max_iter){
    print(paste("Iteration", iter, sep = " "))
    pass <- FALSE
    while (!pass) {
      i <- sample(1:n, 1)  # random row
      j<- sample(1:(m-1),1) # random column
      if (Dold[i,j] %in% block & Dold[i,j+1] %in% block){
        if (which(block==Dold[i,j],arr.ind = TRUE)[1]==which(block==Dold[i,j+1],arr.ind = TRUE)[1]){
          pass <- TRUE
        }
      }
      else {pass<-TRUE}
    }
    Dnew <- Dold
    # exchange component j with component (j+1)
    Dnew[i,j] <- Dold[i,j+1]
    Dnew[i,j+1] <- Dold[i, j]
    XY <- design_to_PWO(rbind(Dold[i,], Dnew[i,]))
    XY <- XY[,-filter$cons_cols]
    x <- XY[1,]
    y <- XY[2,]
    F1 <- cbind(y,-x)
    F2 <- cbind(y,x)
    I2 <- diag(2)
    Cnew <- Cold - Cold%*%F1%*%solve(I2 + t(F2)%*%Cold%*%F1)%*%t(F2)%*%Cold
    Inew <- trace(Cnew%*%B)
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
simulated_annealing <- function(D, B, block, max_iter = 1e4){
  Dold <- D
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- design_to_PWO(Dold)
  cons <-cons_pair(block)
  filter <- PWO_filter(cons=cons, Z=Xold)
  Xold <- Xold[,-filter$cons_cols]  #Calculate the constrainted I-efficiency
  Cold <- solve(t(Xold)%*%Xold/nrow(Xold))
  Iold <- trace(Cold%*%B)
  Ieffs <- numeric(max_iter)
  for (iter in 1:max_iter) {
    print(paste("Iteration", iter, sep = " "))
    pass <- FALSE
    while (!pass) {
      i <- sample(1:n, 1)  # random row
      j<- sample(1:(m-1),1) # random column
      if (Dold[i,j] %in% block & Dold[i,j+1] %in% block){
        if (which(block==Dold[i,j],arr.ind = TRUE)[1]==which(block==Dold[i,j+1],arr.ind = TRUE)[1]){
          pass <- TRUE
        }
      }
      else {pass<-TRUE}
    }
    Dnew <- Dold
    # exchange component j with component (j+1)
    Dnew[i,j] <- Dold[i,j+1]
    Dnew[i,j+1] <- Dold[i, j]
    XY <- design_to_PWO(rbind(Dold[i,], Dnew[i,]))
    XY <- XY[,-filter$cons_cols]
    x <- XY[1,]
    y <- XY[2,]
    F1 <- cbind(y,-x)
    F2 <- cbind(y,x)
    I2 <- diag(2)
    Cnew <- Cold - Cold%*%F1%*%solve(I2 + t(F2)%*%Cold%*%F1)%*%t(F2)%*%Cold
    Inew <- trace(Cnew%*%B)
    u = runif(1)
    if(u<exp(-(Inew-Iold)*iter)){
      Iold <- Inew
      Dold <- Dnew
      Cold <- Cnew
    }
    Ieffs[iter] <- Iold
  }
  plot(1:max_iter, Ieffs, type = "l", ylab = "I-efficiency", xlab = "Iteration")
  return(Dold)
}

#############################################################################
#################### From here we start the case study ######################
#############################################################################

#To begin with we do the unconstrained full model

m = 6
Df <- generate_full_design(m) #Calculate the nlist

job = c(1,2,3,4,5,6)
a = c(5,3,2,10,3,1)
c = c(4,5,3,1,6,2)
Y = rep(0,nrow(Df))
for (i in 1:nrow(Df)) {
  y = 0
  for (j in 1:ncol(Df)) {
    t = 0
    for (k in 1:j) {
      t = t+a[Df[i,k]]
    }
    y = y + c[Df[i,j]]*(t^2)
  }
  Y[i] = y 
}
#Y = Y + rnorm(nrow(Df),0,100)
Xf <- design_to_PWO(Df) 
data_full = as.data.frame(cbind(Xf,Y))
full_model = lm(Y~.,data=data_full)
summary(full_model)

sigma = 0


# Then we do the extreme case



m = 6
b = 2
mu = 0
e = 1
set.seed(1234)
block <- matrix(1,rep(0,4),c(2,3,4,5,6),nrow=2,byrow = TRUE) #Generate block matrices
Df <- generate_full_design(m) #Calculate the nlist

n_sim <- 100
rank_RE <- rep(NA,n_sim)
rank_SA <- rep(NA,n_sim)
rank_rand <- rep(NA,n_sim)
MSE_list_e <-matrix(NA,nrow = n_sim,ncol = 3)

# input the dataset, add a white noise of N(0,2)
job = c(1,2,3,4,5,6)
a = c(5,3,2,10,3,1)
c = c(4,5,3,1,6,2)
Y = rep(0,nrow(Df))
for (i in 1:nrow(Df)) {
  y = 0
  for (j in 1:ncol(Df)) {
    t = 0
    for (k in 1:j) {
      t = t+a[Df[i,k]]
    }
    y = y + c[Df[i,j]]*(t^2)
  }
  Y[i] = y 
}
#Y_obs = Y + sample(c(1,-1),nrow(Df),replace = TRUE)*rnorm(nrow(Df),100,1)
Y = (Y-mean(Y))/sd(Y)
Y_obs = Y #+ rnorm(nrow(Df),0,sigma)
Xf <- design_to_PWO(Df) 
cons <- cons_pair(block)
filter <- PWO_filter(cons, Xf)
Df <- Df[-filter$cons_rows,] 
Y <- Y[-filter$cons_rows]
Y_obs <- Y_obs[-filter$cons_rows]
rank = sort(Y,decreasing = FALSE)
Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
case3_cons = as.data.frame(Xf[,-1])
Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
n = 15

iter = 1
while (iter<=n_sim){
  start_found <- FALSE
  # we need to find a starting design with n rows that is not singular
  while(!start_found){
    D_start <- Df[sample(1:nrow(Df),n),]
    X_start <- design_to_PWO(D_start)
    X_start <- X_start[,-filter$cons_cols]
    M_start <- t(X_start)%*%X_start
    if(!is.singular.matrix(M_start)){
      start_found <- TRUE
      print("Found initial matrix")
    }
  }
  #errors are likely to occur, hence we apply a trycatch
  skip_to_next <- FALSE
  tryCatch(D_out_RE<-random_exchange(D_start,Bf,block=block,max_iter=40),error=function(e){skip_to_next<<-TRUE})
  if(skip_to_next) { next }   
  X_out_RE <- design_to_PWO(D_out_RE)
  X_out_RE <- X_out_RE[,-filter$cons_cols]
  skip_to_next <- FALSE
  tryCatch(D_out_SA <- simulated_annealing(D_start,Bf,block = block,max_iter=40),error=function(e){skip_to_next<<-TRUE})
  if(skip_to_next) { next } 
  X_out_SA <- design_to_PWO(D_out_SA)
  X_out_SA <- X_out_SA[,-filter$cons_cols] 
  # generate the regression models and calculate ranks of obtained optimal designs
  # random exchange
  Y_out_RE <- rep(0,n)
  # separate regressed and predicted data
  for (i in 1:n) {
    Y_out_RE[i] = Y_obs[which(apply(Df, 1, function(x) identical(D_out_RE[i,],x)))]
  }
  X_out_RE = X_out_RE[,-1]
  case3_RE = as.data.frame(cbind(X_out_RE,Y_out_RE))
  cons_RE_model = lm(Y_out_RE ~ ., data = case3_RE)
  RE_pred = predict(cons_RE_model,newdata = case3_cons)
  rank_RE[iter] = which(rank == min(Y[which(RE_pred==min(RE_pred))]))
  # simulated annealing
  Y_out_SA <- rep(0,n)
  # separate regressed and predicted data
  for (i in 1:n) {
    Y_out_SA[i] = Y_obs[which(apply(Df, 1, function(x) identical(D_out_SA[i,],x)))]
  }
  X_out_SA = X_out_SA[,-1]
  case3_SA = as.data.frame(cbind(X_out_SA,Y_out_SA))
  cons_SA_model = lm(Y_out_SA ~ ., data = case3_SA)
  SA_pred = predict(cons_SA_model,newdata = case3_cons)
  rank_SA[iter] = which(rank == min(Y[which(SA_pred==min(SA_pred))]))
  #complete random
  row_out_rand = sample(nrow(Xf),n)
  Y_out_rand = Y_obs[row_out_rand]
  X_out_rand = Xf[row_out_rand,-1]
  case3_rand = as.data.frame(cbind(X_out_rand,Y_out_rand))
  cons_rand_model = lm(Y_out_rand ~ ., data = case3_rand)
  rand_pred = predict(cons_rand_model,newdata = case3_cons)
  rank_rand[iter] = which(rank == min(Y[which(rand_pred==min(rand_pred))]))
  MSE = c(mean(sum((Y-RE_pred)^2)),mean(sum((Y-SA_pred)^2)),mean(sum(Y-rand_pred)^2))
  MSE_list_e[iter,] = MSE
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
data_extreme = list(rank_rand,rank_RE,rank_SA,result)
hist(data_extreme[[1]])
hist(data_extreme[[2]])
hist(data_extreme[[3]])
length(which(data_extreme[[1]]<=3))
length(which(data_extreme[[2]]<=3))
length(which(data_extreme[[3]]<=3))
length(which(data_extreme[[1]]==1))
length(which(data_extreme[[2]]==1))
length(which(data_extreme[[3]]==1))




#Then we do the average case study




m = 6
b = 2
mu = 3
e = 0
block <- matrix(c(1,0,4,6),nrow=2,byrow = TRUE) #Generate block matrices
Df <- generate_full_design(m) #Calculate the nlist
n_sim <- 100
rank_RE <- rep(NA,n_sim)
rank_SA <- rep(NA,n_sim)
rank_rand <- rep(NA,n_sim)
MSE_list_ne <-matrix(NA,nrow = n_sim,ncol = 3)

set.seed(1234)
# input the dataset, add a white noise of N(0,2)
job = c(1,2,3,4,5,6)
a = c(5,3,2,10,3,1)
c = c(4,5,3,1,6,2)
Y = rep(0,nrow(Df))
for (i in 1:nrow(Df)) {
  y = 0
  for (j in 1:ncol(Df)) {
    t = 0
    for (k in 1:j) {
      t = t+a[Df[i,k]]
    }
    y = y + c[Df[i,j]]*(t^2)
  }
  Y[i] = y
}
#Y_obs = Y + sample(c(1,-1),nrow(Df),replace = TRUE)*rnorm(nrow(Df),500,1)
Y = (Y-mean(Y))/sd(Y)
Y_obs = Y + rnorm(nrow(Df),0,sigma)
Xf <- design_to_PWO(Df) 
cons <- cons_pair(block)
filter <- PWO_filter(cons, Xf)
Df <- Df[-filter$cons_rows,] 
Y <- Y[-filter$cons_rows]
Y_obs <- Y_obs[-filter$cons_rows]
rank = sort(Y,decreasing = FALSE)
Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
case3_cons = as.data.frame(Xf[,-1])
Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
n = 15

iter = 1
while (iter<=n_sim){
  start_found <- FALSE
  # we need to find a starting design with n rows that is not singular
  while(!start_found){
    D_start <- Df[sample(1:nrow(Df),n),]
    X_start <- design_to_PWO(D_start)
    X_start <- X_start[,-filter$cons_cols]
    M_start <- t(X_start)%*%X_start
    if(!is.singular.matrix(M_start)){
      start_found <- TRUE
      print("Found initial matrix")
    }
  }
  #errors are likely to occur, hence we apply a trycatch
  skip_to_next <- FALSE
  tryCatch(D_out_RE<-random_exchange(D_start,Bf,block=block,max_iter=40),error=function(e){skip_to_next<<-TRUE})
  if(skip_to_next) { next }   
  X_out_RE <- design_to_PWO(D_out_RE)
  X_out_RE <- X_out_RE[,-filter$cons_cols]
  skip_to_next <- FALSE
  tryCatch(D_out_SA <- simulated_annealing(D_start,Bf,block=block,max_iter=40),error=function(e){skip_to_next<<-TRUE})
  if(skip_to_next) { next } 
  X_out_SA <- design_to_PWO(D_out_SA)
  X_out_SA <- X_out_SA[,-filter$cons_cols] 
  # generate the regression models and calculate ranks of obtained optimal designs
  # random exchange
  Y_out_RE <- rep(0,n)
  # separate regressed and predicted data
  for (i in 1:n) {
    Y_out_RE[i] = Y_obs[which(apply(Df, 1, function(x) identical(D_out_RE[i,],x)))]
  }
  X_out_RE = X_out_RE[,-1]
  case3_RE = as.data.frame(cbind(X_out_RE,Y_out_RE))
  cons_RE_model = lm(Y_out_RE ~ ., data = case3_RE)
  RE_pred = predict(cons_RE_model,newdata = case3_cons)
  rank_RE[iter] = which(rank == min(Y[which(RE_pred==min(RE_pred))]))
  # simulated annealing
  Y_out_SA <- rep(0,n)
  # separate regressed and predicted data
  for (i in 1:n) {
    Y_out_SA[i] = Y_obs[which(apply(Df, 1, function(x) identical(D_out_SA[i,],x)))]
  }
  X_out_SA = X_out_SA[,-1]
  case3_SA = as.data.frame(cbind(X_out_SA,Y_out_SA))
  cons_SA_model = lm(Y_out_SA ~ ., data = case3_SA)
  SA_pred = predict(cons_SA_model,newdata = case3_cons)
  rank_SA[iter] = which(rank == min(Y[which(SA_pred==min(SA_pred))]))
  #complete random
  row_out_rand = sample(nrow(Xf),n)
  Y_out_rand = Y_obs[row_out_rand]
  X_out_rand = Xf[row_out_rand,-1]
  case3_rand = as.data.frame(cbind(X_out_rand,Y_out_rand))
  cons_rand_model = lm(Y_out_rand ~ ., data = case3_rand)
  rand_pred = predict(cons_rand_model,newdata = case3_cons)
  rank_rand[iter] = which(rank == min(Y[which(rand_pred==min(rand_pred))]))
  MSE = c(mean(sum((Y-RE_pred)^2)),mean(sum((Y-SA_pred)^2)),mean(sum(Y-rand_pred)^2))
  MSE_list_ne[iter,] = MSE
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
data_nonextreme = list(rank_rand,rank_RE,rank_SA,result)
hist(data_nonextreme[[1]])
hist(data_nonextreme[[2]])
hist(data_nonextreme[[3]])
length(which(data_nonextreme[[1]]<=3))
length(which(data_nonextreme[[2]]<=3))
length(which(data_nonextreme[[3]]<=3))
length(which(data_nonextreme[[1]]==1))
length(which(data_nonextreme[[2]]==1))
length(which(data_nonextreme[[3]]==1))

summary(MSE_list_e)
summary(MSE_list_ne)