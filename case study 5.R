library(combinat)
library(matrixcalc)

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
random_exchange <- function(D, D_full, B, block, max_iter = 1e4){
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
    Dnew <- Dold
    pass <- FALSE
    while (!pass) {
      i <- sample(1:n, 1)  # random row
      j<- sample(1:(m-1),1) # random column
      # exchange component j with component (j+1)
      Dnew[i,j] <- Dold[i,j+1]
      Dnew[i,j+1] <- Dold[i, j]
      if (sum(which(apply(D_full, 1, function(x) all.equal(x, Dnew[i,])) == "TRUE"))>0){
        pass <-TRUE
      }
    }
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
simulated_annealing <- function(D, D_full, B, block, max_iter = 1e4){
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
    Dnew <- Dold
    pass <- FALSE
    while (!pass) {
      i <- sample(1:n, 1)  # random row
      j<- sample(1:(m-1),1) # random column
      # exchange component j with component (j+1)
      Dnew[i,j] <- Dold[i,j+1]
      Dnew[i,j+1] <- Dold[i, j]
      if (sum(which(apply(D_full, 1, function(x) all.equal(x, Dnew[i,])) == "TRUE"))>0){
        pass <-TRUE
      }
    }
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


# read in data
case2 <- read.csv("case2.csv",header = TRUE)
case2$s6 = case2$block+5
case2$s7 = -case2$block+8

D_cons = as.matrix(cbind(case2[,2:6],case2[,8:9]))

X_cons = design_to_PWO(D_cons)

block = rbind(c(1,2,3,4,5),c(6,7,0,0,0))
cons = cons_pair(block)

# constrain the data, then fit a full model
filter = PWO_filter(Z=X_cons,cons=cons)
X_cons = X_cons[,-filter$cons_cols]

# fit a full model with constraints
Y = case2$Y
case2 = as.data.frame(cbind(X_cons,Y))
case2 = case2[,-1]
full_model = lm(Y ~ ., data = case2)
summary(full_model)

# generate initial matrix
B = t(X_cons)%*%X_cons
#nlist = c(ncol(X_cons)+1,ncol(X_cons)+2,ncol(X_cons)+3)
set.seed(1234)
n = 15
start_found <- FALSE
# we need to find a starting design with n rows that is not singular
while(!start_found){
  D_start <- D_cons[sample(1:nrow(D_cons),n),]
  X_start <- design_to_PWO(D_start)
  X_start <- X_start[,-filter$cons_cols]
  M_start <- t(X_start)%*%X_start
  if(!is.singular.matrix(M_start)){
    start_found <- TRUE
    print("Found initial matrix")
  }
}
#random exchange algorithm
D_out_RE <- random_exchange(D_start, D_full=D_cons, B, block = block,max_iter = 100)
X_out_RE <- design_to_PWO(D_out_RE)
X_out_RE <- X_out_RE[,-filter$cons_cols]

#simulated annealing algorithm
D_out_SA <- simulated_annealing(D_start, D_full=D_cons, B, block = block,max_iter = 100)
X_out_SA <- design_to_PWO(D_out_SA)
X_out_SA <- X_out_SA[,-filter$cons_cols]


# now we fit the models on the subsets

Y_out_RE <- rep(0,n)
for (i in 1:n) {
  Y_out_RE[i] = case1[case1$X1==D_out_RE[i,1] & case1$X2==D_out_RE[i,2] & case1$X3==D_out_RE[i,3],5] 
}
case1_RE = as.data.frame(cbind(X_out_RE,Y_out_RE))

# we did not notice collinearity here
cons_RE_model = lm(Y_out_RE ~ z12 + z14 + z23 + z24 + z34, data = case1_RE)
RE_pred = predict(cons_RE_model,newdata = case1_cons)



Y_out_SA <- rep(0,n)
for (i in 1:n) {
  Y_out_SA[i] = case1[case1$X1==D_out_SA[i,1] & case1$X2==D_out_SA[i,2] & case1$X3==D_out_SA[i,3],5] 
}
case1_SA = as.data.frame(cbind(X_out_SA,Y_out_SA))
# 
cons_SA_model = lm(Y_out_SA ~ z12 + z14 +z23 + z24 + z34, data = case1_SA)
SA_pred = predict(cons_SA_model,newdata = case1_cons)

# We use a total random subset to compare

case1_rand = case1[sample(nrow(case1_cons),n),]
cons_rand_model = lm(Y_out_SA ~ z12 + z14 +z23 + z24 + z34, data = case1_rand)
rand_pred = predict(cons_rand_model,newdata = case1_cons)

SSE = c(sum((Y-RE_pred)^2),sum((Y-SA_pred)^2),sum(Y-rand_pred)^2)
SSE
summary(cons_full_model)
summary(cons_RE_model)
summary(cons_SA_model)
summary(cons_rand_model)

# from the result we can see that SA and RE model are both much better than a totally random model.

# The real constrained optimal order is 3 1 2 4 
# The optimal order by full constrained model is 3 4 1 2 ranked 3rd of the 12
# The optimal order by total random is 4 3 1 2 ranked 7th of the 12
# The optimal order by RE is 3 1 4 2 ranked 5th of the 12
# The optimal order by SA is 3 4 1 2 or 4 3 1 2 ranked 3rd and 7th respectively


