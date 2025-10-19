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

# random_exchange
# randomly selects adjacent pairs of components to swap
# performs the swap if it improves I optimality
# inputs: 
# D: a n by m matrix of permutations of (1:m)
# B: the B matrix for the I opt criterion

# max_iter: max number of iterations, default is 1e4
random_exchange_v2 <- function(n, Df, B, inter, D_start, max_iter = 1e4){
  name <- c('z12','z23','z34','z45','z56','z67','z78')
  # start_found <- FALSE
  # while(!start_found){
  #    D_start <- Df[sample(1:nrow(Df),n),]
  #  X_start <- design_to_PWO(D_start)
  #  for (i in 1:inter) {
  #    X_start = cbind(X_start,X_start[,which(colnames(X_start)==name[i])]*X_start[,which(colnames(X_start)==name[i+1])])
  #  }
  #  M_start <- t(X_start)%*%X_start
  #  if(!is.singular.matrix(M_start)){
  #    start_found <- TRUE
  #    print("Found initial matrix")
  #  }
  #}
  
  
  Dold <- D_start
  X_start <- design_to_PWO(D_start)
  for(i in 1:inter){
      X_start <- cbind(X_start,X_start[,which(colnames(X_start)==name[i])]*X_start[,which(colnames(X_start)==name[i+1])])
  }
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- X_start
  Iold <- I_efficiency(Xold,B)
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
    for (i in 1:inter) {
      Xnew = cbind(Xnew,Xnew[,which(colnames(Xnew)==name[i])]*Xnew[,which(colnames(Xnew)==name[i+1])])
    }
    Inew <- trace(Xnew%*%B)
    if(Inew < Iold){
      Iold <- Inew
      Dold <- Dnew
    }
    Ieffs[iter] <- Iold
  }
  plot(1:max_iter, Ieffs, type = "l", ylab = "I-efficiency", xlab = "Iteration")
  Ieffs = Ieffs[iter]
  result = Dold
  return(result)
}

### simulated annealing
simulated_annealing_v2 <- function(n, Df, B, inter, D_start, max_iter = 1e4){
  name <- c('z12','z23','z34','z45','z56','z67','z78')
  # start_found <- FALSE
  # while(!start_found){
  #    D_start <- Df[sample(1:nrow(Df),n),]
  #  X_start <- design_to_PWO(D_start)
  #  for (i in 1:inter) {
  #    X_start = cbind(X_start,X_start[,which(colnames(X_start)==name[i])]*X_start[,which(colnames(X_start)==name[i+1])])
  #  }
  #  M_start <- t(X_start)%*%X_start
  #  if(!is.singular.matrix(M_start)){
  #    start_found <- TRUE
  #    print("Found initial matrix")
  #  }
  # }
 
  
  Dold <- D_start
  X_start <- design_to_PWO(D_start)
  for(i in 1:inter){
      X_start <- cbind(X_start,X_start[,which(colnames(X_start)==name[i])]*X_start[,which(colnames(X_start)==name[i+1])])
  }
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- X_start
  Iold <- I_efficiency(Xold,B)
  Ieffs <- numeric(max_iter)
  for (iter in 1:max_iter) {
    print(paste("Iteration", iter, sep = " "))
    i <- sample(1:n, 1)  # random row
    j<- sample(1:(m-1),1) # random column
    Dnew <- Dold
    # exchange component j with component (j+1)
    Dnew[i,j] <- Dold[i,j+1]
    Dnew[i,j+1] <- Dold[i, j]
    Xnew <- design_to_PWO(Dnew)
    for (i in 1:inter) {
      Xnew = cbind(Xnew,Xnew[,which(colnames(Xnew)==name[i])]*Xnew[,which(colnames(Xnew)==name[i+1])])
    }
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
  result = Dold
  return(result)
}




# n: target sample size
# Df: full design (if filtered, use the filter argument)
# vote_n: design size that voters vote on
# vote_num: number of voters
# filter: contains pairwise constraints
# Bf: B-matrix to evaluate I-optimality criterion
# max_iter: max number of iterations for RE
# v2 update: Instead of using the full design and calling identical() to count votes,
# we create a set of unique rows that is updated as more "candidates" for "votes" come in
RE_voting_v3 <- function(n,Df, vote_n, vote_num, inter, Bf, max_iter=1e4){
  
  id = matrix(apply(Df,1,function(x) as.numeric(paste(x,sep='',collapse=''))),ncol=1)
  vote_RE = matrix(0,nrow = vote_n,ncol = vote_num)
  
  for (voter in 1:vote_num) {
    
    skip_to_next <- FALSE
    tryCatch(D_out_RE<-random_exchange_v2(n=n,Df=Df,B=Bf,inter=inter,max_iter),error=function(e){skip_to_next<<-TRUE})
    if(skip_to_next) { next }   
    
    candidate_id = matrix(apply(D_out_RE,1,function(x) as.numeric(paste(x,sep='',collapse=''))),ncol=1)
    #print(candidate_matrix)
    
    vote_RE[,voter] = apply(candidate_id, 1, function(x) which(id==x))
    
  }
  D_out_RE <- Df[as.numeric(names(sort(table(vote_RE),decreasing=TRUE)[1:n])),]
  
  return(D_out_RE)
  
}

# n: target sample size
# Df: full design (if filtered, use the filter argument)
# vote_n: design size that voters vote on
# vote_num: number of voters
# filter: contains pairwise constraints
# Bf: B-matrix to evaluate I-optimality criterion
# max_iter: max number of iterations for RE
# v2 update: Instead of using the full design and calling identical() to count votes,
# we create a set of unique rows that is updated as more "candidates" for "votes" come in
SA_voting_v3 <- function(n,Df, vote_n, vote_num, inter, Bf, D_start, max_iter=1e4){
  
  id = matrix(apply(Df,1,function(x) as.numeric(paste(x,sep='',collapse=''))),ncol=1)
  vote_SA = matrix(0,nrow = vote_n,ncol = vote_num)
  
  for (voter in 1:vote_num) {
    
    skip_to_next <- FALSE
    tryCatch(D_out_SA<-simulated_annealing_v2(n=n,Df=Df,B=Bf,inter=inter, D_start = D_start, max_iter=max_iter),error=function(e){skip_to_next<<-TRUE})
    if(skip_to_next) { next }   
    
    candidate_id = matrix(apply(D_out_SA,1,function(x) as.numeric(paste(x,sep='',collapse=''))),ncol=1)
    #print(candidate_matrix)
    
    vote_SA[,voter] = apply(candidate_id, 1, function(x) which(id==x))
    
  }
  D_out_SA <- Df[as.numeric(names(sort(table(vote_SA),decreasing=TRUE)[1:n])),]
  
  return(D_out_SA)
  
}


insertion_sort <- function(n, Df, B, max_iter = 5){
  
  
  name <- c('z12','z23','z34','z45','z56','z67','z78')
  start_found <- FALSE
  # we need to find a starting design with vote_n rows that is not singular
  while(!start_found){
    D_start <- Df[sample(1:nrow(Df),n),]
    X_start <- design_to_PWO(D_start)
    for(i in 1:inter){
      X_start <- cbind(X_start,X_start[,which(colnames(X_start)==name[i])]*X_start[,which(colnames(X_start)==name[i+1])])
    }
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
  for(i in 1:inter){
      Xold <- cbind(Xold,Xold[,which(colnames(Xold)==name[i])]*Xold[,which(colnames(Xold)==name[i+1])])
    }
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
          Dtmp <- Dold
          Dtmp[row,j] <- Dold[row,j-1]
          Dtmp[row,j-1] <- Dold[row, j]
          XY <- design_to_PWO(rbind(Dold[row,], Dtmp[row,]))
          for(i in 1:inter){
            XY <- cbind(XY,XY[,which(colnames(XY)==name[i])]*XY[,which(colnames(XY)==name[i+1])])
        }
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




#--------------------------------------------
### Comparison between Random Exchange and Simulated Annealing
#--------------------------------------------
comparison <- function(m=6,inter=2,vote_num=50,tot_iter=100,
                       max_iter = 2000){
  
  Df <- generate_full_design(m) #Calculate the nlist
  Xf <- design_to_PWO(Df)
  name <- c('z12','z23','z34','z45','z56','z67','z78')
  for (i in 1:inter) {
    Xf = cbind(Xf,Xf[,which(colnames(Xf)==name[i])]*Xf[,which(colnames(Xf)==name[i+1])])
  }
  
  nlist = pmin(c(90,100,110),round(c(nrow(Xf)/4,nrow(Xf)/2,3*nrow(Xf)/4))) #list of sample sizes
  if (nlist[1]<=ncol(Xf)){
    nlist=c(ncol(Xf)+1,ncol(Xf)+2,ncol(Xf)+3)
  }
  
  I_list_SA = matrix(NA,nrow = tot_iter, ncol = length(nlist))
  I_list_SA_vote = matrix(NA,nrow = tot_iter, ncol = length(nlist))
  I_list_Federov = matrix(NA,nrow = tot_iter, ncol = length(nlist))
  I_list_insertion_sort = matrix(NA,nrow = tot_iter, ncol = length(nlist))
  Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
  I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
  
  for (i in 1:length(nlist)) {
  
    n.cores <- 30
    my.cluster <- parallel::makeCluster(n.cores)
    registerDoParallel(my.cluster)
  
    I_effs <- foreach(iter = 1:tot_iter, .packages = c("combinat", "matrixcalc", "reshape2", "AlgDesign",
                                            "ggplot2", "retry"), .combine = 'rbind',
                                            .export=c("design_to_PWO", "SA_voting_v3", "simulated_annealing_v2", "I_efficiency",
                                            "random_exchange_v2", "trace", "generate_full_design","insertion_sort")) %dopar% { 
      print(paste(c('n=',nlist[i],',iteration',iter),sep='',collapse=''))
      n = nlist[i]
      #optFederov
      out=retry(optFederov(data=Xf[,-1],nTrials=n,criterion="I"),when = 'Singular design')
      D_out_Federov = Df[out$rows,]
      X_out <- design_to_PWO(D_out_Federov)
      for (j in 1:inter) {
        X_out = cbind(X_out,X_out[,which(colnames(X_out)==name[j])]*X_out[,which(colnames(X_out)==name[j+1])])
      }
      I_out <- I_efficiency(X_out, Bf)
      # I_list_Federov[iter,i]=I_full/I_out
      I_Federov=I_full/I_out
      # insertion sort
      D_out_insertion_sort = insertion_sort(n=n, Df = Df, B = Bf)
      X_out <- design_to_PWO(D_out_insertion_sort)
      for (j in 1:inter) {
        X_out = cbind(X_out,X_out[,which(colnames(X_out)==name[j])]*X_out[,which(colnames(X_out)==name[j+1])])
      }
      I_out <- I_efficiency(X_out, Bf)
      # I_list_RE[iter,i]=I_full/I_out
      I_insertion = I_full/I_out
      #voting
      vote_n = n
      D_out_SA_vote <- SA_voting_v3(n=n,Df=Df, vote_n=vote_n, vote_num=vote_num, inter=inter, Bf=Bf, max_iter=1000, D_start = D_out_Federov)
      X_out <- design_to_PWO(D_out_SA_vote)
      for (j in 1:inter) {
        X_out = cbind(X_out,X_out[,which(colnames(X_out)==name[j])]*X_out[,which(colnames(X_out)==name[j+1])])
      }
      I_out <- I_efficiency(X_out, Bf)
      I_SA_vote = I_full/I_out
      #simulated annealing
      D_out_SA = simulated_annealing_v2(n=n, D=Df, B=Bf, inter=inter,max_iter=max_iter, D_start = D_out_Federov)
      X_out <- design_to_PWO(D_out_SA)
      for (j in 1:inter) {
        X_out = cbind(X_out,X_out[,which(colnames(X_out)==name[j])]*X_out[,which(colnames(X_out)==name[j+1])])
      }
      I_out <- I_efficiency(X_out, Bf)
      I_SA = I_full/I_out
      # c(I_Federov, I_RE, I_RE_vote, I_SA)
      c(I_Federov, I_insertion, I_SA_vote, I_SA)
    }
    
    I_list_Federov[,i] = I_effs[,1]
    I_list_insertion_sort[,i] = I_effs[,2]
    I_list_SA_vote[,i] = I_effs[,3]
    I_list_SA[,i] = I_effs[,4]
    
    parallel::stopCluster(my.cluster)
    
  }
  preplot = cbind(matrix(I_list_Federov,ncol = 1),matrix(I_list_insertion_sort,ncol = 1),matrix(I_list_SA_vote,ncol = 1),matrix(I_list_SA,ncol = 1))
  preplot = as.data.frame(cbind(preplot,rep(nlist,each=tot_iter)))
  colnames(preplot)=c("optFederov","InsertionSort","SA_voting","SimulatedAnnealing","n")
  preplot$n=as.factor(preplot$n)
  preplot_long = melt(preplot,id.vars = c("n"))
  return(preplot_long)
} 




# Generate the full data frame for all possible combination of characteristics

M <- c(8) #Number of components
inter <- c(4)

set.seed(1234)
# n.cores <- 30
vote_num <- 30
# max_iter <- 5000
# my.cluster <- parallel::makeCluster(n.cores)
# registerDoParallel(my.cluster)

# nsim = length(M)*length(inter)
# DF <- as.matrix(expand.grid(inter, M))
# out_list <- foreach(iter = 1:nsim, .packages = c("combinat", "matrixcalc", "reshape2", "AlgDesign",
#                                            "ggplot2", "retry")) %dopar% {
#     m = DF[iter,2]
#     i = DF[iter,1]
#     comparison_list <- comparison(m=m,inter=i,vote_num=vote_num,tot_iter=100, max_iter = max_iter)
#     namestr <- paste("comparison_list_m",m,"interactions_",i,".csv",sep = "")
#     write.csv(comparison_list, namestr)
#     comparison_list 
 #}

# saveRDS(out_list, file = "interaction_sim_results.R")

#for(m in M){

    #miter = 10000
    #for(i in inter){
    #    comparison_list <- comparison(m=m,inter=i,vote_num=vote_num,tot_iter=30, max_iter = miter)
    #    namestr <- paste("comparison_list_m",m,"interactions_",i,".csv",sep = "")
    #    write.csv(comparison_list, namestr)
    #}
    #print(paste("m = ",m," DONE", sep = ""))
#}


comparison_list1 <- comparison(m=6,inter=2,vote_num=vote_num,tot_iter=30, max_iter = 10000)
write.csv(comparison_list1, "final_comparison_list_m6interactions2.csv")
print("m = 6 DONE")

comparison_list2 <- comparison(m=7,inter=3,vote_num=vote_num,tot_iter=30, max_iter = 10000)
write.csv(comparison_list2, "final_comparison_list_m7interactions3.csv")
print("m = 7 DONE")

comparison_list3 <- comparison(m=8,inter=4,vote_num=vote_num,tot_iter=30, max_iter = 10000)
write.csv(comparison_list3, "final_comparison_list_m8interactions4.csv")
print("m = 8 DONE")

# parallel::stopCluster(my.cluster)






