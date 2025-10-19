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

# random_exchange
# randomly selects adjacent pairs of components to swap
# performs the swap if it improves I optimality
# TODO: Try to improve this. Move of the iterations are being wasted..
# inputs: 
# D: a n by m matrix of permutations of (1:m)
# B: the B matrix for the I opt criterion

# max_iter: max number of iterations, default is 1e4
random_exchange_v2 <- function(n, Df, B, inter, max_iter = 1e4){
  name <- c('z12','z23','z34','z45','z56','z67','z78')
  start_found <- FALSE
  while(!start_found){
    D_start <- Df[sample(1:nrow(Df),n),]
    X_start <- design_to_PWO(D_start)
    for (i in 1:inter) {
      X_start = cbind(X_start,X_start[,which(colnames(X_start)==name[i])]*X_start[,which(colnames(X_start)==name[i+1])])
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
simulated_annealing_v2 <- function(n, Df, B, inter, max_iter = 1e4){
  name <- c('z12','z23','z34','z45','z56','z67','z78')
  start_found <- FALSE
  while(!start_found){
    D_start <- Df[sample(1:nrow(Df),n),]
    X_start <- design_to_PWO(D_start)
    for (i in 1:inter) {
      X_start = cbind(X_start,X_start[,which(colnames(X_start)==name[i])]*X_start[,which(colnames(X_start)==name[i+1])])
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
SA_voting_v3 <- function(n,Df, vote_n, vote_num, inter, Bf, max_iter=1e4){
  
  id = matrix(apply(Df,1,function(x) as.numeric(paste(x,sep='',collapse=''))),ncol=1)
  vote_SA = matrix(0,nrow = vote_n,ncol = vote_num)
  
  for (voter in 1:vote_num) {
    
    skip_to_next <- FALSE
    tryCatch(D_out_SA<-random_exchange_v2(n=n,Df=Df,B=Bf,inter=inter,max_iter),error=function(e){skip_to_next<<-TRUE})
    if(skip_to_next) { next }   
    
    candidate_id = matrix(apply(D_out_SA,1,function(x) as.numeric(paste(x,sep='',collapse=''))),ncol=1)
    #print(candidate_matrix)
    
    vote_SA[,voter] = apply(candidate_id, 1, function(x) which(id==x))
    
  }
  D_out_SA <- Df[as.numeric(names(sort(table(vote_SA),decreasing=TRUE)[1:n])),]
  
  return(D_out_SA)
  
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
  
  I_list_RE_vote = matrix(NA,nrow = tot_iter, ncol = length(nlist))
  I_list_RE = matrix(NA,nrow = tot_iter, ncol = length(nlist))
  I_list_SA = matrix(NA,nrow = tot_iter, ncol = length(nlist))
  I_list_Federov = matrix(NA,nrow = tot_iter, ncol = length(nlist))
  Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
  I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
  
  for (i in 1:length(nlist)) {
    for (iter in 1:tot_iter) { 
      print(paste(c('n=',nlist[i],',iteration',iter),sep='',collapse=''))
      n = nlist[i]
      #optFederov
      out=retry(optFederov(data=Xf[,-1],nTrials=n,criterion="I"),when = 'Singular design')
      D_out = Df[out$rows,]
      X_out <- design_to_PWO(D_out)
      for (j in 1:inter) {
        X_out = cbind(X_out,X_out[,which(colnames(X_out)==name[j])]*X_out[,which(colnames(X_out)==name[j+1])])
      }
      I_out <- I_efficiency(X_out, Bf)
      I_list_Federov[iter,i]=I_full/I_out
      #random exchange
      D_out_RE = random_exchange_v2(n=n, Df=Df, B=Bf, inter=inter,max_iter=max_iter)
      X_out <- design_to_PWO(D_out_RE)
      for (j in 1:inter) {
        X_out = cbind(X_out,X_out[,which(colnames(X_out)==name[j])]*X_out[,which(colnames(X_out)==name[j+1])])
      }
      I_out <- I_efficiency(X_out, Bf)
      I_list_RE[iter,i]=I_full/I_out
      #random exchange voting
      vote_n = n
      D_out_RE_vote <- RE_voting_v3(n=n,Df=Df, vote_n=vote_n, vote_num=vote_num, inter=inter, Bf=Bf, max_iter=50)
      X_out <- design_to_PWO(D_out_RE_vote)
      for (j in 1:inter) {
        X_out = cbind(X_out,X_out[,which(colnames(X_out)==name[j])]*X_out[,which(colnames(X_out)==name[j+1])])
      }
      I_out <- I_efficiency(X_out, Bf)
      I_list_RE_vote[iter,i]=I_full/I_out
      #simulated annealing
      D_out_SA = simulated_annealing_v2(n=n, D=Df, B=Bf, inter=inter,max_iter=max_iter)
      X_out <- design_to_PWO(D_out_SA)
      for (j in 1:inter) {
        X_out = cbind(X_out,X_out[,which(colnames(X_out)==name[j])]*X_out[,which(colnames(X_out)==name[j+1])])
      }
      I_out <- I_efficiency(X_out, Bf)
      I_list_SA[iter,i]=I_full/I_out
    }
  }
  preplot = cbind(matrix(I_list_Federov,ncol = 1),matrix(I_list_RE,ncol = 1),matrix(I_list_RE_vote,ncol = 1),matrix(I_list_SA,ncol = 1))
  preplot = as.data.frame(cbind(preplot,rep(nlist,each=tot_iter)))
  colnames(preplot)=c("optFederov","RandomExchange","RE_voting","SimulatedAnnealing","n")
  preplot$n=as.factor(preplot$n)
  preplot_long = melt(preplot,id.vars = c("n"))
  return(preplot_long)
} 




# Generate the full data frame for all possible combination of characteristics

M <- c(6,7,8) #Number of components
inter <- c(2,3,4)

set.seed(1234)
t1 = Sys.time()


comparison_list1 <- comparison(m=6,inter=2,vote_num=50,tot_iter=2, max_iter = 1000)
write.csv(comparison_list1, "df1_inter_out.csv")
# comparison_list2 <- comparison(m=7,inter=3,vote_num=50,tot_iter=100)
# write.csv(comparison_list2, "df2_inter_out.csv")
# comparison_list3 <- comparison(m=8,inter=4,vote_num=50,tot_iter=100)
# write.csv(comparison_list3, "df3_inter_out.csv")




print("m = 6, 7, 8 DONE")

t2 = Sys.time()
t2-t1

# n.cores <- 30
# vote_num <- 50
# my.cluster <- parallel::makeCluster(n.cores)
# registerDoParallel(my.cluster)
# parallel::stopCluster(my.cluster)






