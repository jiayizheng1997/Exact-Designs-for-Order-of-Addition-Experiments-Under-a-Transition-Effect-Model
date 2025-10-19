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
simulated_annealing_v2 <- function(n, Df, B, block, filter, max_iter = 1e4){
  
  
  start_found <- FALSE
  # we need to find a starting design with vote_n rows that is not singular
  
  
  if(!missing(filter)){
    Df <- Df[-filter$cons_rows,] 
    Xf <- design_to_PWO(Df)
    Xf <- Xf[,-filter$cons_cols]
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



# n: target sample size
# Df: full design (if filtered, use the filter argument)
# vote_n: design size that voters vote on
# vote_num: number of voters
# filter: contains pairwise constraints
# Bf: B-matrix to evaluate I-optimality criterion
# max_iter: max number of iterations for RE
# v2 update: Instead of using the full design and calling identical() to count votes,
# we create a set of unique rows that is updated as more "candidates" for "votes" come in
SA_voting_v3 <- function(n,Df, vote_n, vote_num, block, filter, Bf, max_iter=1e4){
  
  id = matrix(apply(Df,1,function(x) as.numeric(paste(x,sep='',collapse=''))),ncol=1)
  vote_SA = matrix(0,nrow = vote_n,ncol = vote_num)
  
  for (voter in 1:vote_num) {
    
    skip_to_next <- FALSE
    tryCatch(D_out_SA<-simulated_annealing_v2(n, Df,Bf,block=block, filter = filter,max_iter),error=function(e){skip_to_next<<-TRUE})
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
comparison <- function(data_row,tot_iter=100, max_iter = 5000){
  m = data_row[1] #number of components
  b = data_row[2] #number of blocks
  mu = data_row[3] #number of free components
  e = data_row[4] #extremeness
  alg = data_row[5] #random exchange or simulated annealing
  block <- cons_block(m,b,e,mu) #Generate block matrices
  Df <- generate_full_design(m) #Calculate the nlist
  Xf <- design_to_PWO(Df)
  cons <- cons_pair(block)
  filter <- PWO_filter(cons, Xf)
  Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
  nlist = pmin(c(90,100,110),round(c(nrow(Xf)/4,nrow(Xf)/2,3*nrow(Xf)/4))) #list of sample sizes
  # nlist = pmin(c(90,100,110),round(c(nrow(Xf)/4,nrow(Xf)/2,3*nrow(Xf)/4)) + 2) #list of sample sizes
  if (nlist[1]<=ncol(Xf)){
    nlist=c(ncol(Xf)+1,ncol(Xf)+2,ncol(Xf)+3)
  }
  I_list = matrix(NA,nrow = tot_iter, ncol = length(nlist))
  Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
  I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
  I_name <- c(paste("m=",m,",n=",nlist,"b=",b,",mu=",mu,",e=",e,c("RE","SA")[alg],sep=""))
  colnames(I_list) = I_name
  for (i in 1:length(nlist)) {
    for (iter in 1:tot_iter) { 
      print(paste(c('n=',nlist[i],',iteration',iter),sep=''))
      n = nlist[i]
      skip_to_next <- FALSE
      if (alg==1){
        tryCatch(D_out <- insertion_sort(n, Df, B = Bf, block, filter, max_iter = 10), error = function(e){skip_to_next<<-TRUE})
      }
      if (alg==2){
        tryCatch(D_out<-simulated_annealing_v2(n=n, Df=Df, B=Bf, block=block, filter=filter,max_iter = max_iter),error=function(e){skip_to_next<<-TRUE})
      }
      if (alg==3){
          # out=optFederov(data=Xf[,-1],nTrials=n,criterion="I",nRepeats = 25)
          out=retry(optFederov(data=Xf[,-1],nTrials=n,criterion="I", nRepeats = 2),when = 'Singular design')
          Df2 <- Df[-filter$cons_rows,] 
          D_out = Df2[out$rows,]
      }
      if (alg==4){
          D_out <- SA_voting_v3(n=n,Df=Df,Bf=Bf,block=block,filter=filter,vote_n=n,vote_num=30,max_iter=1000)
      }
      if(skip_to_next) { next }   
      X_out <- design_to_PWO(D_out)
      X_out <- X_out[,-filter$cons_cols]
      I_out <- I_efficiency(X_out, Bf)
      I_list[iter,i]=I_full/I_out
    }
  }
  preplot = cbind(matrix(I_list,ncol = 1),rep(nlist,each=tot_iter))
  if (alg==1){colnames(preplot)=c("InsertionSort","n")}
  if (alg==2){colnames(preplot)=c("SimulatedAnnealing","n")}
  if (alg==3){colnames(preplot)=c("optFederov","n")}
  if (alg==4){colnames(preplot)=c("SA_voting","n")}
  return(preplot)
  
  # preplot = cbind(matrix(I_list,ncol = 1),matrix(I_list_vote,ncol = 1))
  # preplot = as.data.frame(cbind(preplot,rep(nlist,each=tot_iter)))
  # if (alg==1){colnames(preplot)=c("RandomExchange","RandomExchangeVote","n")}
  # else {colnames(preplot)=c("SimulatedAnnealing","SimulatedAnnealingVote","n")}
  # preplot$n=as.factor(preplot$n)
  # preplot_long = melt(preplot,id.vars = c("n"))
  # return(preplot_long)
  # p <- ggplot(preplot_long,aes(x=n,y=value,fill=n),plot=FALSE)+geom_boxplot()+facet_wrap(~variable)+ylim(0.5,1)+theme(axis.title.y=element_blank(),axis.title.x=element_blank())+theme(legend.position="none")
  # png(paste("m",m,"b",b,"mu",mu,"e",e,c("RE.png","SA.png")[alg],sep = ""))
  # print(p)
  # dev.off()
  # return(p)
} 
#p=comparison(DF[4,])



# Generate the full data frame for all possible combination of characteristics

M <- c(6,7,8) #Number of components
B <- c(2,3) #Number of blocks
mu <- c(0,1,2,3) #Number of free components
E <- c(0,1) #extreme/non-extreme division
alg <- c(1,2,3,4) 
DF <- matrix(NA,ncol = 5, nrow = length(B)*length(M)*length(mu)*length(E)*length(alg)) #The experimental matrix
i = 1
for (i1 in 1:length(M)) {
  for (i2 in 1:length(B)) {
    if (B[i2]==2){
      for (i3 in 1:length(mu)) {
        for (i4 in 1:length(E)) {
          for (i5 in 1:length(alg)) {
            DF[i,1] <- M[i1]
            DF[i,2] <- B[i2]
            DF[i,3] <- mu[i3]
            DF[i,4] <- E[i4]
            DF[i,5] <- alg[i5]
            i = i+1
          }
        }
      }
    }
    else {
      for (i3 in 1:length(mu)) {
        for (i5 in 1:length(alg)) {
          DF[i,1] <- M[i1]
          DF[i,2] <- B[i2]
          DF[i,3] <- mu[i3]
          DF[i,4] <- 0
          DF[i,5] <- alg[i5]
          i = i+1
        }
      }
    }
  }
}
DF <- DF[which(DF[,1]!=0),]
colnames(DF) <- c("components","blocks","free components","extremeness","algorithm")

set.seed(1234)
n.cores <- 30
max_iter <- 10000
my.cluster <- parallel::makeCluster(n.cores)
registerDoParallel(my.cluster)

# comparison_list1 <- foreach(l = 1:24,
#                             .packages = c("combinat", "matrixcalc", "reshape2", "AlgDesign",
#                                           "ggplot2", "retry")) %dopar% {
#                                             comparison(DF[l,], max_iter = max_iter)
#                                           }

#write.csv(comparison_list1, paste("df1_out", max_iter, "max_iter.csv", sep = ""))
#print("m = 6 DONE")


comparison_list2 <- foreach(l = 25:48,
                             .packages = c("combinat", "matrixcalc", "reshape2",  "AlgDesign",
                                           "ggplot2", "retry")) %dopar% {
                                             cbind(l,comparison(DF[l,], max_iter = max_iter))
                                           }


write.csv(comparison_list2, paste("df2_out", max_iter, "max_iter.csv", sep = ""))
print("m = 7 DONE")


#comparison_list3 <- foreach(l = 49:72,
#                             .packages = c("combinat", "matrixcalc", "reshape2",  "AlgDesign",
#                                           "ggplot2", "retry")) %dopar% {
#                                             comparison(DF[l,], max_iter = max_iter)
#                                           }

#write.csv(comparison_list3, paste("df3_out", max_iter, "max_iter.csv", sep = ""))
#print("m = 8 DONE")


parallel::stopCluster(my.cluster)