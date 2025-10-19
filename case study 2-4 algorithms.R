library(combinat)
library(matrixcalc)
library(doParallel)
library(foreach)
library(AlgDesign)
library(ggplot2)
library(reshape2)
library(retry)

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


### simulated annealing
simulated_annealing_v2 <- function(n, Df, B, max_iter = 1e4){
  start_found <- FALSE
  while(!start_found){
    D_start <- Df[sample(1:nrow(Df),n),]
    X_start <- design_to_PWO(D_start)
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
  # Cold <- solve(t(Xold)%*%Xold/nrow(Xold))
  Cold <- solve(t(Xold)%*%Xold)
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
    tmpmat <- I2 + t(F2)%*%Cold%*%F1
    invmat <- 1/(tmpmat[1,1]*tmpmat[2,2] - tmpmat[1,2]*tmpmat[2,1])*matrix(c(tmpmat[2,2], -tmpmat[1,2],-tmpmat[2,1],tmpmat[1,1]),
                                                                           byrow = TRUE, nrow = 2, ncol = 2)
    # Cnew <- Cold - Cold%*%F1%*%solve(I2 + t(F2)%*%Cold%*%F1)%*%t(F2)%*%Cold
    Cnew <- Cold - Cold%*%F1%*%invmat%*%t(F2)%*%Cold
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


SA_voting <- function(n,Df, vote_n, vote_num, Bf, max_iter=1e4){
  
  id = matrix(apply(Df,1,function(x) as.numeric(paste(x,sep='',collapse=''))),ncol=1)
  vote_SA = matrix(0,nrow = vote_n,ncol = vote_num)
  
  for (voter in 1:vote_num) {
    
    skip_to_next <- FALSE
    tryCatch(D_out_SA<-simulated_annealing_v2(n=n,Df=Df,B=Bf,max_iter=max_iter),error=function(e){skip_to_next<<-TRUE})
    if(skip_to_next) { next }   
    
    candidate_id = matrix(apply(D_out_SA,1,function(x) as.numeric(paste(x,sep='',collapse=''))),ncol=1)
    #print(candidate_matrix)
    
    vote_SA[,voter] = apply(candidate_id, 1, function(x) which(id==x))
    
  }
  D_out_SA <- Df[as.numeric(names(sort(table(vote_SA),decreasing=TRUE)[1:n])),]
  
  return(D_out_SA)
  
}




insertion_sort <- function(n, Df, B, block, filter, max_iter = 5){
  
  
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
  Xold <- Xold[,-filter$cons_cols]  #Calculate the constrained I-efficiency
  # Cold <- solve(t(Xold)%*%Xold/nrow(Xold))
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
          pass <- FALSE
          if (Dold[row,j] %in% block & Dold[row,j-1] %in% block){
            if (which(block==Dold[row,j],arr.ind = TRUE)[1]==which(block==Dold[row,j-1],arr.ind = TRUE)[1]){
              pass <- TRUE
            }
          }  else {pass<-TRUE}
          if(!pass){
            next
          }
          Dtmp <- Dold
          Dtmp[row,j] <- Dold[row,j-1]
          Dtmp[row,j-1] <- Dold[row, j]
          XY <- design_to_PWO(rbind(Dold[row,], Dtmp[row,]))
          XY <- XY[,-filter$cons_cols]
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

comparison <- function(Y,block,sigma,n_sim,n){
  
  Y_obs = Y + rnorm(nrow(Df),0,sigma)
  Xf <- design_to_PWO(Df) 
  cons <- cons_pair(block)
  filter <- PWO_filter(cons, Xf)
  Y <- Y[-filter$cons_rows]
  Y_obs <- Y_obs[-filter$cons_rows]
  rank = sort(Y,decreasing = FALSE)
  Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
  # fix: we need to filter Df when finding the matching Y's
  Df_filtered <- Df[-filter$cons_rows,]
  case3_cons = as.data.frame(Xf[,-1])
  Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
  I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
  
  # simulated annealing
  # skip_to_next <- FALSE
  # tryCatch(D_out_SA <- simulated_annealing_v2(n=n, Df, Bf, block=block, filter, max_iter = 40),error=function(e){skip_to_next<<-TRUE})
  # if(skip_to_next) { next } 
  # D_out_SA <- simulated_annealing_v2(n=n, Df, Bf, block=block, filter, max_iter = 40)
  D_out_SA <- simulated_annealing_v2(n=n, Df, Bf, block=block, filter, max_iter = 1000)
  X_out_SA <- design_to_PWO(D_out_SA)
  X_out_SA <- X_out_SA[,-filter$cons_cols] 
  Y_out_SA <- rep(0,nrow(D_out_SA))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_SA)) {
    Y_out_SA[i] = Y_obs[which(apply(Df_filtered, 1, function(x) identical(D_out_SA[i,],x)))]
  }
  X_out_SA = X_out_SA[,-1]
  case3_SA = as.data.frame(cbind(X_out_SA,Y_out_SA))
  cons_SA_model = lm(Y_out_SA ~ ., data = case3_SA)
  SA_pred = predict(cons_SA_model,newdata = case3_cons)
  rank_SA = which(rank == min(Y[which(SA_pred==min(SA_pred))]))[1]
  
  # SA voting
  vote_n = n
  vote_num = 200
  #skip_to_next <- FALSE
  #tryCatch(D_out_voting <- SA_voting_v3(n=n,Df, vote_n, vote_num, Bf, block=block, filter, max_iter=1e4),error=function(e){skip_to_next<<-TRUE})
  #if(skip_to_next) { next } 
  # D_out_voting <- SA_voting_v3(n=n,Df, vote_n, vote_num, Bf, block=block, filter, max_iter= 40)
  D_out_voting <- SA_voting_v3(n=n,Df, vote_n, vote_num, Bf, block=block, filter, max_iter= 100)
  X_out_voting <- design_to_PWO(D_out_voting)
  X_out_voting <- X_out_voting[,-filter$cons_cols] 
  Y_out_voting <- rep(0,nrow(D_out_voting))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_voting)) {
    Y_out_voting[i] = Y_obs[which(apply(Df_filtered, 1, function(x) identical(D_out_voting[i,],x)))]
  }
  X_out_voting = X_out_voting[,-1]
  case3_voting = as.data.frame(cbind(X_out_voting,Y_out_voting))
  cons_voting_model = lm(Y_out_voting ~ ., data = case3_voting)
  voting_pred = predict(cons_voting_model,newdata = case3_cons)
  rank_voting = which(rank == min(Y[which(voting_pred==min(voting_pred))]))[1]
  
  # insertion sorting
  #skip_to_next <- FALSE
  #tryCatch(D_out_IS <- insertion_sort(n=n, Df, Bf, block=block, filter, max_iter = 5),error=function(e){skip_to_next<<-TRUE})
  #if(skip_to_next) { next } 
  D_out_IS <- insertion_sort(n=n, Df, Bf, block=block, filter, max_iter = 5)
  # D_out_IS <- insertion_sort(n=n, Df, Bf, block=block, filter, max_iter = 10)
  X_out_IS <- design_to_PWO(D_out_IS)
  X_out_IS <- X_out_IS[,-filter$cons_cols] 
  Y_out_IS <- rep(0,nrow(D_out_IS))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_IS)) {
    Y_out_IS[i] = Y_obs[which(apply(Df_filtered, 1, function(x) identical(D_out_IS[i,],x)))]
  }
  X_out_IS = X_out_IS[,-1]
  case3_IS = as.data.frame(cbind(X_out_IS,Y_out_IS))
  cons_IS_model = lm(Y_out_IS ~ ., data = case3_IS)
  IS_pred = predict(cons_IS_model,newdata = case3_cons)
  rank_IS = which(rank == min(Y[which(IS_pred==min(IS_pred))]))[1]
  
  # optFederov
  out=retry(optFederov(data=Xf[,-1],nTrials=n,criterion="I"),when = 'Singular design')
  Df2 <- Df[-filter$cons_rows,] 
  D_out_Fed = Df2[out$rows,]
  X_out_Fed <- design_to_PWO(D_out_Fed)
  X_out_Fed <- X_out_Fed[,-filter$cons_cols] 
  Y_out_Fed <- rep(0,nrow(D_out_Fed))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_Fed)) {
    Y_out_Fed[i] = Y_obs[which(apply(Df_filtered, 1, function(x) identical(D_out_Fed[i,],x)))]
  }
  X_out_Fed = X_out_Fed[,-1]
  case3_Fed = as.data.frame(cbind(X_out_Fed,Y_out_Fed))
  cons_Fed_model = lm(Y_out_Fed ~ ., data = case3_Fed)
  Fed_pred = predict(cons_Fed_model,newdata = case3_cons)
  rank_Fed = which(rank == min(Y[which(Fed_pred==min(Fed_pred))]))[1]
  
  #complete random
  row_out_rand = sample(nrow(Xf),n)
  Y_out_rand = Y_obs[row_out_rand]
  X_out_rand = Xf[row_out_rand,-1]
  case3_rand = as.data.frame(cbind(X_out_rand,Y_out_rand))
  cons_rand_model = lm(Y_out_rand ~ ., data = case3_rand)
  rand_pred = predict(cons_rand_model,newdata = case3_cons)
  rank_rand = which(rank == min(Y[which(rand_pred==min(rand_pred))]))[1]
  
  MSE = c(mean(sum((Y-SA_pred)^2)),mean(sum((Y-voting_pred)^2)),mean(sum((Y-IS_pred)^2)),mean(sum((Y-Fed_pred)^2)),mean(sum(Y-rand_pred)^2))
  # MSE = c(mean(Y-SA_pred)^2,mean((Y-voting_pred)^2),mean((Y-IS_pred)^2),mean((Y-Fed_pred)^2),mean(Y-rand_pred)^2)
  return(c(rank_SA,rank_voting,rank_IS,rank_Fed,rank_rand,MSE))
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

n_sim <- 100
Y = (Y-mean(Y))/sd(Y)
n = 15
block <- matrix(c(1,rep(0,4),c(2:6)),nrow=2,byrow = TRUE) #Generate block matrices

set.seed(1234)

n.cores <- 10
my.cluster <- parallel::makeCluster(n.cores)
registerDoParallel(my.cluster)

# Then we do the extreme case

data_extreme <- foreach(l = 1:100,.errorhandling = 'remove',
                            .packages = c("combinat", "matrixcalc", "reshape2",  "AlgDesign",
                                          "ggplot2", "retry")) %dopar% {
                                            comparison(Y,block,sigma,n_sim,n)
                                          }

   
success = length(data_extreme)
data_extreme = matrix(unlist(data_extreme),nrow = success,byrow = TRUE)
result = matrix(c(length(which(data_extreme[,1]<=3)),length(which(data_extreme[,2]<=3)),length(which(data_extreme[,3]<=3)),length(which(data_extreme[,4]<=3)),length(which(data_extreme[,5]<=3)),length(which(data_extreme[,1]==1)),length(which(data_extreme[,2]==1)),length(which(data_extreme[,3]==1)),length(which(data_extreme[,4]==1)),length(which(data_extreme[,5]==1))),nrow = 2,byrow = TRUE)
rownames(result) = c("leq 3","equal 1")
colnames(result) = c("SA","SA_voting","InsertionSorting","optFedrov","Random")
result

summary(data_extreme[,6:10])

write.csv(data_extreme,paste(c("case2studyextremesigma",sigma,".csv"),sep="",collapse=""))



#Then we do the average case study
block <- matrix(c(1,0,4,6),nrow=2,byrow = TRUE) #Generate block matrices
data_nonextreme <- foreach(l = 1:100,.errorhandling = 'remove',
                        .packages = c("combinat", "matrixcalc", "reshape2",  "AlgDesign",
                                      "ggplot2", "retry")) %dopar% {
                                        comparison(Y,block,sigma,n_sim,n)
                                      }

success = length(data_nonextreme)
data_nonextreme = matrix(unlist(data_nonextreme),nrow = success,byrow = TRUE)
result = matrix(c(length(which(data_nonextreme[,1]<=3)),length(which(data_nonextreme[,2]<=3)),length(which(data_nonextreme[,3]<=3)),length(which(data_nonextreme[,4]<=3)),length(which(data_nonextreme[,5]<=3)),length(which(data_nonextreme[,1]==1)),length(which(data_nonextreme[,2]==1)),length(which(data_nonextreme[,3]==1)),length(which(data_nonextreme[,4]==1)),length(which(data_nonextreme[,5]==1))),nrow = 2,byrow = TRUE)
rownames(result) = c("leq 3","equal 1")
colnames(result) = c("SA","SA_voting","InsertionSorting","optFedrov","Random")
result

summary(data_nonextreme[,6:10])

write.csv(data_nonextreme, paste(c("case2studynonextremesigma",sigma,".csv"),sep="",collapse=""))

sigma = 0.01
# Then we do the extreme case
block <- matrix(c(1,rep(0,4),c(2:6)),nrow=2,byrow = TRUE) #Generate block matrices
data_extreme <- foreach(l = 1:100,.errorhandling = 'remove',
                        .packages = c("combinat", "matrixcalc", "reshape2",  "AlgDesign",
                                      "ggplot2", "retry")) %dopar% {
                                        comparison(Y,block,sigma,n_sim,n)
                                      }


success = length(data_extreme)
data_extreme = matrix(unlist(data_extreme),nrow = success,byrow = TRUE)
result = matrix(c(length(which(data_extreme[,1]<=3)),length(which(data_extreme[,2]<=3)),length(which(data_extreme[,3]<=3)),length(which(data_extreme[,4]<=3)),length(which(data_extreme[,5]<=3)),length(which(data_extreme[,1]==1)),length(which(data_extreme[,2]==1)),length(which(data_extreme[,3]==1)),length(which(data_extreme[,4]==1)),length(which(data_extreme[,5]==1))),nrow = 2,byrow = TRUE)
rownames(result) = c("leq 3","equal 1")
colnames(result) = c("SA","SA_voting","InsertionSorting","optFedrov","Random")
result

summary(data_extreme[,6:10])

write.csv(data_extreme,paste(c("case2studyextremesigma",sigma,".csv"),sep="",collapse=""))

#Then we do the average case study
block <- matrix(c(1,0,4,6),nrow=2,byrow = TRUE) #Generate block matrices
data_nonextreme <- foreach(l = 1:100,.errorhandling = 'remove',
                           .packages = c("combinat", "matrixcalc", "reshape2",  "AlgDesign",
                                         "ggplot2", "retry")) %dopar% {
                                           comparison(Y,block,sigma,n_sim,n)
                                         }

success = length(data_nonextreme)
data_nonextreme = matrix(unlist(data_nonextreme),nrow = success,byrow = TRUE)
result = matrix(c(length(which(data_nonextreme[,1]<=3)),length(which(data_nonextreme[,2]<=3)),length(which(data_nonextreme[,3]<=3)),length(which(data_nonextreme[,4]<=3)),length(which(data_nonextreme[,5]<=3)),length(which(data_nonextreme[,1]==1)),length(which(data_nonextreme[,2]==1)),length(which(data_nonextreme[,3]==1)),length(which(data_nonextreme[,4]==1)),length(which(data_nonextreme[,5]==1))),nrow = 2,byrow = TRUE)
rownames(result) = c("leq 3","equal 1")
colnames(result) = c("SA","SA_voting","InsertionSorting","optFedrov","Random")
result

summary(data_nonextreme[,6:10])

write.csv(data_nonextreme, paste(c("case2studynonextremesigma",sigma,".csv"),sep="",collapse=""))

sigma = 0.1
# Then we do the extreme case
block <- matrix(c(1,rep(0,4),c(2:6)),nrow=2,byrow = TRUE) #Generate block matrices
data_extreme <- foreach(l = 1:100,.errorhandling = 'remove',
                        .packages = c("combinat", "matrixcalc", "reshape2",  "AlgDesign",
                                      "ggplot2", "retry")) %dopar% {
                                        comparison(Y,block,sigma,n_sim,n)
                                      }


success = length(data_extreme)
data_extreme = matrix(unlist(data_extreme),nrow = success,byrow = TRUE)
result = matrix(c(length(which(data_extreme[,1]<=3)),length(which(data_extreme[,2]<=3)),length(which(data_extreme[,3]<=3)),length(which(data_extreme[,4]<=3)),length(which(data_extreme[,5]<=3)),length(which(data_extreme[,1]==1)),length(which(data_extreme[,2]==1)),length(which(data_extreme[,3]==1)),length(which(data_extreme[,4]==1)),length(which(data_extreme[,5]==1))),nrow = 2,byrow = TRUE)
rownames(result) = c("leq 3","equal 1")
colnames(result) = c("SA","SA_voting","InsertionSorting","optFedrov","Random")
result

summary(data_extreme[,6:10])

write.csv(data_extreme,paste(c("case2studyextremesigma",sigma,".csv"),sep="",collapse=""))

#Then we do the average case study
block <- matrix(c(1,0,4,6),nrow=2,byrow = TRUE) #Generate block matrices
data_nonextreme <- foreach(l = 1:100,.errorhandling = 'remove',
                           .packages = c("combinat", "matrixcalc", "reshape2",  "AlgDesign",
                                         "ggplot2", "retry")) %dopar% {
                                           comparison(Y,block,sigma,n_sim,n)
                                         }

success = length(data_nonextreme)
data_nonextreme = matrix(unlist(data_nonextreme),nrow = success,byrow = TRUE)

result = matrix(c(length(which(data_nonextreme[,1]<=3)),length(which(data_nonextreme[,2]<=3)),length(which(data_nonextreme[,3]<=3)),length(which(data_nonextreme[,4]<=3)),length(which(data_nonextreme[,5]<=3)),length(which(data_nonextreme[,1]==1)),length(which(data_nonextreme[,2]==1)),length(which(data_nonextreme[,3]==1)),length(which(data_nonextreme[,4]==1)),length(which(data_nonextreme[,5]==1))),nrow = 2,byrow = TRUE)
rownames(result) = c("leq 3","equal 1")
colnames(result) = c("SA","SA_voting","InsertionSorting","optFedrov","Random")
result


summary(data_nonextreme[,6:10])

write.csv(data_nonextreme, paste(c("case2studynonextremesigma",sigma,".csv"),sep="",collapse=""))

parallel::stopCluster(my.cluster)






