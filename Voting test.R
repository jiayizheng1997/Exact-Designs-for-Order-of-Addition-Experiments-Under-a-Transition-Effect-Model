library(combinat)
library(matrixcalc)
library(doParallel)
library(foreach)

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


# random_exchange_v2
# randomly selects adjacent pairs of components to swap
# performs the swap if it improves I optimality
# inputs: 
# n: target sample size
# Df: full design
# B: the B matrix for the I opt criterion
# block: the blockwise constraint matrix
# max_iter: max number of iterations, default is 1e4
# filter: PWO_filter
random_exchange_v2 <- function(n, Df, B, block, filter, max_iter = 1e4){
  start_found <- FALSE
  # we need to find a starting design with vote_n rows that is not singular
  
  
  if(!missing(filter)){
    Df <- Df[-filter$cons_rows,] 
  }
  
  
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
  
  Dold <- D_start
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



### simulated annealing
simulated_annealing_v2 <- function(n, Df, B, block, filter, max_iter = 1e4){
  
  
  if(!missing(filter)){
    Df <- Df[-filter$cons_rows,] 
  }
  
  
  
  start_found <- FALSE
  # we need to find a starting design with vote_n rows that is not singular
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
  
  Dold <- D_start
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


# n: target sample size
# Df: full design (if filtered, use the filter argument)
# vote_n: design size that voters vote on
# vote_num: number of voters
# filter: contains pairwise constraints
# Bf: B-matrix to evaluate I-optimality criterion
# max_iter: max number of iterations for RE
RE_voting <- function(n,Df, vote_n, vote_num, filter, Bf, max_iter){
  
  
  vote_RE = matrix(0,nrow = vote_n,ncol = vote_num)
  for (voter in 1:vote_num) {
    
    skip_to_next <- FALSE
    tryCatch(D_out_RE<-random_exchange_v2(n, Df,Bf,block=block, filter = filter,max_iter=max_iter),error=function(e){skip_to_next<<-TRUE})
    if(skip_to_next) { next }   
    for (i in 1:nrow(D_out_RE)){
      vote_RE[i,voter] = which(apply(Df, 1, function(x) isTRUE(all.equal(D_out_SA[i,],x))))
    }
    
    
    
  }
  D_out_RE <- Df[t(as.numeric(names(sort(table(vote_RE),decreasing=TRUE)[1:n]))),]
  
  return(D_out_RE)
  
}


SA_voting <- function(n,Df, vote_n, vote_num, filter, Bf, max_iter){
  
  
  vote_SA = matrix(0,nrow = vote_n,ncol = vote_num)
  for (voter in 1:vote_num) {
    
    skip_to_next <- FALSE
    tryCatch(D_out_SA<-simulated_annealing_v2(n, Df,Bf,block=block, filter = filter,max_iter=max_iter),error=function(e){skip_to_next<<-TRUE})
    if(skip_to_next) { next }   
    for (i in 1:nrow(D_out_SA)){
      vote_SA[i,voter] = which(apply(Df, 1, function(x) isTRUE(all.equal(D_out_SA[i,],x))))
    }
    
    
    
  }
  D_out_SA <- Df[t(as.numeric(names(sort(table(vote_SA),decreasing=TRUE)[1:n]))),]
  
  return(D_out_SA)
  
}

########## Test Code Here ##########
set.seed(1234)
m = 6
n = 20
vote_n = n
vote_num = 100
n_sim = 10
Df = generate_full_design(m) 
Xf <- design_to_PWO(Df) 
block = matrix(c(c(2,3,4,5,6),1,rep(0,4)),nrow=2,byrow = TRUE) 
cons = cons_pair(block)
filter = PWO_filter(cons, Xf)

Xf2 <- Xf[-filter$cons_rows,-filter$cons_cols]
Bf <- t(Xf2)%*%Xf2

I_full <- I_efficiency(Xf2, Bf)

rel_I_voting = numeric(n_sim)
rel_I_RE = numeric(n_sim)

my_cluster = parallel::makeCluster(1, n.cores = 11)

out1 = foreach(i = 1:n_sim, .packages = c("combinat","matrixcalc"),
               .errorhandling = "stop", .combine = "rbind") %dopar% {
  
  D_RE_voting = RE_voting(n = n, Df = Df, vote_n = vote_n, vote_num = vote_num, filter = filter, Bf = Bf,
                          max_iter = 100)
  X_RE_voting = design_to_PWO(D_RE_voting)
  X_RE_voting = X_RE_voting[,-filter$cons_cols]
  I_RE_voting = I_efficiency(X_RE_voting, Bf)
  
  
  # compare to just random exchange
  D_RE = random_exchange_v2(n, Df, Bf, block, filter, max_iter = 100)
  X_RE = design_to_PWO(D_RE)
  X_RE = X_RE[,-filter$cons_cols]
  I_RE = I_efficiency(X_RE, Bf)
  
  c(I_full/I_RE_voting, I_full/I_RE)

  
               }

out2 = foreach(i = 1:n_sim, .packages = c("combinat","matrixcalc"),
               .errorhandling = "stop", .combine = "rbind") %dopar% {
                 
                 D_SA_voting = SA_voting(n = n, Df = Df, vote_n = vote_n, vote_num = vote_num, filter = filter, Bf = Bf,
                                         max_iter = 100)
                 X_SA_voting = design_to_PWO(D_SA_voting)
                 X_SA_voting = X_SA_voting[,-filter$cons_cols]
                 I_SA_voting = I_efficiency(X_SA_voting, Bf)
                 
                 
                 # compare to just random exchange
                 D_SA = simulated_annealing_v2(n, Df, Bf, block, filter, max_iter = 100)
                 X_SA = design_to_PWO(D_SA)
                 X_SA = X_SA[,-filter$cons_cols]
                 I_SA = I_efficiency(X_SA, Bf)
                 
                 c(I_full/I_SA_voting, I_full/I_SA)
                 
                 
               }

# 
# for(i in 1:n_sim){
#   print(paste("SIMULATION ", i, sep = ""))
#   D_RE_voting = RE_voting(n = n, Df = Df, vote_n = vote_n, vote_num = vote_num, filter = filter, Bf = Bf,
#                           max_iter = 100)
#   X_RE_voting = design_to_PWO(D_RE_voting)
#   X_RE_voting = X_RE_voting[,-filter$cons_cols]
#   I_RE_voting = I_efficiency(X_RE_voting, Bf)
#   
#   
#   # compare to just random exchange
#   D_RE = random_exchange_v2(n, Df, Bf, block, filter, max_iter = 100)
#   X_RE = design_to_PWO(D_RE)
#   X_RE = X_RE[,-filter$cons_cols]
#   I_RE = I_efficiency(X_RE, Bf)
#   
#   
#   rel_I_voting[i] = I_full/I_RE_voting
#   rel_I_RE[i] = I_full/I_RE
#   
#   
# }

parallel::stopCluster(my_cluster)

rel_I_voting = out1[,1]
rel_I_RE = out1[,2]

summary(rel_I_voting)
summary(rel_I_RE)
df_out_RE = data.frame(rel_I_voting, rel_I_RE)

rel_I_voting = out2[,1]
rel_I_SA = out2[,2]

summary(rel_I_voting)
summary(rel_I_SA)
df_out_SA = data.frame(rel_I_voting, rel_I_SA)
df_out = data.frame(df_out_RE,df_out_SA)
write.csv(df_out,"Voting_Sim_Results.csv")


