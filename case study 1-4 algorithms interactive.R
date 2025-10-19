library(combinat)
library(matrixcalc)
library(doParallel)
library(foreach)
library(AlgDesign)
library(ggplot2)
library(reshape2)
library(retry)

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

### simulated annealing
simulated_annealing <- function(n, D, B, int1, int2,  max_iter = 100){
  
  start_found <- FALSE
  # we need to find a starting design with vote_n rows that is not singular
  while(!start_found){
    D_start <- Df[sample(1:nrow(Df),n),]
    X_start <- design_to_PWO(D_start)
    X_start <- cbind(X_start,X_start[,int1]*X_start[,int2])
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

SA_voting_v3 <- function(n, Df, vote_n, vote_num, Bf, int1, int2, max_iter=1e4){
  id = matrix(apply(Df,1,function(x) as.numeric(paste(x,sep='',collapse=''))),ncol=1)
  vote_SA = matrix(0,nrow = vote_n,ncol = vote_num)
  for (voter in 1:vote_num) {
    skip_to_next <- FALSE
    tryCatch(D_out_SA<-simulated_annealing(n=n, D=Df, B=Bf, int1=int1, int2=int2,  max_iter = 1000)[[1]],error=function(e){skip_to_next<<-TRUE})
    if(skip_to_next) { next }   
    candidate_id = matrix(apply(D_out_SA,1,function(x) as.numeric(paste(x,sep='',collapse=''))),ncol=1)
    vote_SA[,voter] = apply(candidate_id, 1, function(x) which(id==x))
  }
  D_out_SA <- Df[as.numeric(names(sort(table(vote_SA),decreasing=TRUE)[1:n])),]
  return(D_out_SA)
}

insertion_sort <- function(n, Df, B, int1, int2, max_iter = 5){
  
  start_found <- FALSE
  # we need to find a starting design with vote_n rows that is not singular
  while(!start_found){
    D_start <- Df[sample(1:nrow(Df),n),]
    X_start <- design_to_PWO(D_start)
    X_start <- cbind(X_start,X_start[,int1]*X_start[,int2])
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
  Xold <- cbind(Xold,Xold[,int1]*Xold[,int2])
  Cold <- solve(t(Xold)%*%Xold)
  Iold <- trace(Cold%*%B)
  for(iter in 1:max_iter){
    print(paste("iteration", iter, sep = " "))
    for(row in 1:n){
      current_row = Dold[row,]
      print(current_row)
      for(i in 2:m){
        for(j in i:2){
          # try swapping elements j and j-1
          # if it is better, perform the swap
          # if the swap is invalid, skip it
          Dtmp <- Dold
          Dtmp[row,j] <- Dold[row,j-1]
          Dtmp[row,j-1] <- Dold[row, j]
          XY <- design_to_PWO(rbind(Dold[row,], Dtmp[row,]))
          XY <- cbind(XY,XY[,int1]*XY[,int2])
          x <- XY[1,]
          y <- XY[2,]
          F1 <- cbind(y,-x)
          F2 <- cbind(y,x)
          I2 <- diag(2)
          # Ctmp <- Cold - Cold%*%F1%*%solve(I2 + t(F2)%*%Cold%*%F1)%*%t(F2)%*%Cold
          tmpmat <- I2 + t(F2)%*%Cold%*%F1
          invmat <- 1/(tmpmat[1,1]*tmpmat[2,2] - tmpmat[1,2]*tmpmat[2,1])*matrix(c(tmpmat[2,2], -tmpmat[1,2],-tmpmat[2,1],tmpmat[1,1]),
                                                                                 byrow = TRUE, nrow = 2, ncol = 2)
          Ctmp <- Cold - Cold%*%F1%*%invmat%*%t(F2)%*%Cold
          Itmp <- trace(Ctmp%*%B)
          if(Itmp < Iold){
            # swap
            Dold <- Dtmp
            Cold <- Ctmp
            Iold <- Itmp
            print(Iold)
          }
        }
      }
    }
  }
  return(Dold)
}

comparison <- function(Y,int1,int2,n_sim,n){
  
  Y_obs = Y #+ rnorm(nrow(Df),0,sigma)
  Xf <- design_to_PWO(Df) 
  rank = sort(Y,decreasing = FALSE)
  Xf <- cbind(Xf,Xf[,int1]*Xf[,int2])
  case1 = as.data.frame(Xf[,-1])
  Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
  I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
  
  # simulated annealing
  D_out_SA <- simulated_annealing(n=n, Df, Bf, int1, int2, max_iter = 1000)[[1]]
  X_out_SA <- design_to_PWO(D_out_SA)
  X_out_SA <- cbind(X_out_SA,X_out_SA[,int1]*X_out_SA[,int2]) 
  Y_out_SA <- rep(0,nrow(D_out_SA))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_SA)) {
    Y_out_SA[i] = Y_obs[which(apply(Df, 1, function(x) identical(D_out_SA[i,],x)))]
  }
  X_out_SA = X_out_SA[,-1]
  case1_SA = as.data.frame(cbind(X_out_SA,Y_out_SA))
  cons_SA_model = lm(Y_out_SA ~ ., data = case1_SA)
  SA_pred = predict(cons_SA_model,newdata = case1)
  rank_SA = which(rank == min(Y[which(SA_pred==min(SA_pred))]))[1]
  
  # SA voting
  vote_n = n
  vote_num = 200
  D_out_voting <- SA_voting_v3(n, Df, vote_n, vote_num, Bf, int1, int2, max_iter=1000)
  X_out_voting <- design_to_PWO(D_out_voting)
  X_out_voting <- cbind(X_out_voting,X_out_voting[,int1]*X_out_voting[,int2])
  Y_out_voting <- rep(0,nrow(D_out_voting))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_voting)) {
    Y_out_voting[i] = Y_obs[which(apply(Df, 1, function(x) identical(D_out_voting[i,],x)))]
  }
  X_out_voting = X_out_voting[,-1]
  case1_voting = as.data.frame(cbind(X_out_voting,Y_out_voting))
  cons_voting_model = lm(Y_out_voting ~ ., data = case1_voting)
  voting_pred = predict(cons_voting_model,newdata = case1)
  rank_voting = which(rank == min(Y[which(voting_pred==min(voting_pred))]))[1]
  
  # insertion sorting
  D_out_IS <- insertion_sort(n=n, Df, Bf, int1, int2, max_iter = 3)
  X_out_IS <- design_to_PWO(D_out_IS)
  X_out_IS <- cbind(X_out_IS,X_out_IS[,int1]*X_out_IS[,int2])
  Y_out_IS <- rep(0,nrow(D_out_IS))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_IS)) {
    Y_out_IS[i] = Y_obs[which(apply(Df, 1, function(x) identical(D_out_IS[i,],x)))]
  }
  X_out_IS = X_out_IS[,-1]
  case1_IS = as.data.frame(cbind(X_out_IS,Y_out_IS))
  cons_IS_model = lm(Y_out_IS ~ ., data = case1_IS)
  IS_pred = predict(cons_IS_model,newdata = case1)
  rank_IS = which(rank == min(Y[which(IS_pred==min(IS_pred))]))[1]
  
  # optFederov
  out=retry(optFederov(data=Xf[,-1],nTrials=n,criterion="I"),when = 'Singular design')
  D_out_Fed = Df[out$rows,]
  X_out_Fed <- design_to_PWO(D_out_Fed)
  X_out_Fed <- cbind(X_out_Fed,X_out_Fed[,int1]*X_out_Fed[,int2])
  Y_out_Fed <- rep(0,nrow(D_out_Fed))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_Fed)) {
    Y_out_Fed[i] = Y_obs[which(apply(Df, 1, function(x) identical(D_out_Fed[i,],x)))]
  }
  X_out_Fed = X_out_Fed[,-1]
  case1_Fed = as.data.frame(cbind(X_out_Fed,Y_out_Fed))
  cons_Fed_model = lm(Y_out_Fed ~ ., data = case1_Fed)
  Fed_pred = predict(cons_Fed_model,newdata = case1)
  rank_Fed = which(rank == min(Y[which(Fed_pred==min(Fed_pred))]))[1]
  
  #complete random
  row_out_rand = sample(nrow(Xf),n)
  Y_out_rand = Y_obs[row_out_rand]
  X_out_rand = Xf[row_out_rand,-1]
  case1_rand = as.data.frame(cbind(X_out_rand,Y_out_rand))
  cons_rand_model = lm(Y_out_rand ~ ., data = case1_rand)
  rand_pred = predict(cons_rand_model,newdata = case1)
  rank_rand = which(rank == min(Y[which(rand_pred==min(rand_pred))]))[1]
  
  MSE = c(mean((Y-SA_pred)^2),mean((Y-voting_pred)^2),mean((Y-IS_pred)^2),mean((Y-Fed_pred)^2),mean((Y-rand_pred)^2))
  # MSE = c(mean(Y-SA_pred)^2,mean((Y-voting_pred)^2),mean((Y-IS_pred)^2),mean((Y-Fed_pred)^2),mean(Y-rand_pred)^2)
  return(c(rank_SA,rank_voting,rank_IS,rank_Fed,rank_rand,MSE))
}


################################################################################
####################### Real Data Starts from Here #############################
################################################################################


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

# now we start calculating the rank with z13z23 interaction added

n = 9
int_X = cbind(Xf,Xf[,3]*Xf[,5])
Bf = t(int_X)%*%int_X
case1_int = as.data.frame(cbind(Y,int_X))
colnames(case1_int)[9] = 'z13z23'
n_sim <- 100
int1 = 3
int2 = 5
Y = data1[,1]

set.seed(1234)

n.cores <- 10
my.cluster <- parallel::makeCluster(n.cores)
registerDoParallel(my.cluster)

# Then we do the extreme case

data_extreme <- foreach(l = 1:100,.errorhandling = 'remove',
                        .packages = c("combinat", "matrixcalc", "reshape2",  "AlgDesign",
                                      "ggplot2", "retry")) %dopar% {
                                        comparison(Y,int1,int2,n_sim,n)
                                      }


success = length(data_extreme)
data_extreme = matrix(unlist(data_extreme),nrow = success,byrow = TRUE)
result = matrix(c(length(which(data_extreme[,1]<=3)),length(which(data_extreme[,2]<=3)),length(which(data_extreme[,3]<=3)),length(which(data_extreme[,4]<=3)),length(which(data_extreme[,5]<=3)),length(which(data_extreme[,1]==1)),length(which(data_extreme[,2]==1)),length(which(data_extreme[,3]==1)),length(which(data_extreme[,4]==1)),length(which(data_extreme[,5]==1))),nrow = 2,byrow = TRUE)
rownames(result) = c("leq 3","equal 1")
colnames(result) = c("SA","SA_voting","InsertionSorting","optFedrov","Random")
result

summary(data_extreme[,6:10])

write.csv(data_extreme,"case1studyinteractive.csv")


# read in data
case1 <- read.csv("case1.csv",header = FALSE,col.names = c('X1','X2','X3','X4','Y'))
case1[,1:4] = case1[,1:4]+1
Df = as.matrix(case1[,1:4])

Xf = design_to_PWO(Df)
data1 = as.data.frame(cbind(as.data.frame(case1)$Y,Xf))
colnames(data1)[1] = 'Y'
n = 9
Bf = t(Xf)%*%Xf
case1_int = as.data.frame(cbind(Y,Xf))
n_sim <- 100
Y = data1[,1]
int1 = 0
int2 = 0
# Then we do the extreme case
data_extreme <- foreach(l = 1:100,.errorhandling = 'remove',
                        .packages = c("combinat", "matrixcalc", "reshape2",  "AlgDesign",
                                      "ggplot2", "retry")) %dopar% {
                                        comparison(Y,int1,int2,n_sim,n)
                                      }


success = length(data_extreme)
data_extreme = matrix(unlist(data_extreme),nrow = success,byrow = TRUE)
result = matrix(c(length(which(data_extreme[,1]<=3)),length(which(data_extreme[,2]<=3)),length(which(data_extreme[,3]<=3)),length(which(data_extreme[,4]<=3)),length(which(data_extreme[,5]<=3)),length(which(data_extreme[,1]==1)),length(which(data_extreme[,2]==1)),length(which(data_extreme[,3]==1)),length(which(data_extreme[,4]==1)),length(which(data_extreme[,5]==1))),nrow = 2,byrow = TRUE)
rownames(result) = c("leq 3","equal 1")
colnames(result) = c("SA","SA_voting","InsertionSorting","optFedrov","Random")
result

summary(data_extreme[,6:10])

write.csv(data_extreme,"case1studynoninteractive.csv")

parallel::stopCluster(my.cluster)

