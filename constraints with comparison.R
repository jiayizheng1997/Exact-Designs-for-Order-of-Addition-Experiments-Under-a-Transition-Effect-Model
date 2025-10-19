# libraries
library(combinat)
library(matrixcalc)
library(ggplot2)
library(reshape2)

#Now we want to input constraints to our current model.
#The constraints should be a matrix of t rows, filled with numbers between 1 and m.
#m is the number of components, and fill the empty entries with zero if necessary.


#This is a function to generate equal sized blocks, I don't recommend using it.
#Instead, we can manually input blocks.
cons_block <- function(m,t){  
  n1=m%/%t #in each block, there would be at least n1 components
  n2=m%%t #in n2 blocks, there would be n1+1 components
  #we want to write a row of length m as a matrix of n1*(n1+1)
  #Fill the matrix of repeated values with zero then filter the zero
  cons_block = matrix(sample(m,size=m,replace=FALSE),nrow = t,ncol = n1+1)
  if (n2>0){for (i in 1:(t-n2)){cons_block[n2+i,n1+1]=0}}
  return(cons_block)
}
cons_block(10,3)
cons_block(8,3)
cons_block(8,2)


#After the generation of block, either by hand or automatically,
#we need to convert them into pairwise indications and therefore filter the constraint rows

#The result would be an n by 2 matrix, where a row (a,b) means a must go in front of b
#B should be a matrix whose row number are constraint blocks
cons_pair <- function(B){ #input the matrix obtained above 
  cons_pair=matrix(c(0,0),ncol=2)
  #get all binary combinations between two rows
  for (i in 1:(nrow(B)-1)) {
    temp=lapply(B[i,],function(X){lapply(B[i+1,],function(Y){c(X,Y)})})
    temp=matrix(unlist(temp),ncol=2,byrow = TRUE)
    cons_pair=rbind(cons_pair,temp)
  }
  cons_pair=cons_pair[apply(cons_pair,1,function(X) all(X!=0)),]
  return(cons_pair) #the final output is a matrix with 2 columns
}
cons_pair(matrix(c(2,3),byrow = TRUE,nrow = 2))
cons_pair(matrix(c(1,0,2,3),byrow = TRUE,nrow = 2))
cons_pair(matrix(c(5,3,2,1,4,6,7,9,8),byrow = TRUE,nrow = 3))


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

# m=8
# block=matrix(c(2,4,3,0),byrow = TRUE,nrow = 2)
# cons=cons_pair(block)
# 
# D1=generate_full_design(m)
# z1=design_to_PWO(D) #full PWO matrix
# filter <- PWO_filter(cons = cons,Z = z1)
# z2=z1[-filter$cons_rows,-filter$cons_cols] #reduced PWO matrix
# D2=D1[-filter$cons_rows,]
# nrow(z2)


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
  Xold <- Xold[,-filter$cons_cols]  #Calculate the constrainted I-efficiency
  Cold <- solve(t(Xold)%*%Xold/nrow(Xold))
  Iold <- trace(Cold%*%B)
  # Iold <- I_efficiency(Xold, B)
  Ieffs <- numeric(max_iter)
  for(iter in 1:max_iter){
    print(paste("Iteration", iter, sep = " "))
    pass <- FALSE
    i <- sample(1:n, 1)  # random row
    j<- sample(1:(m-1),1) # random column
    while (!pass) {
      if (Dold[i,j] %in% block & Dold[i,j+1] %in% block){
        if (which(block==Dold[i,j],arr.ind = TRUE)[1]==which(block==Dold[i,j+1],arr.ind = TRUE)[1]){
        pass <- TRUE
        }
      }
      else {pass<-TRUE}
      i <- sample(1:n, 1)
      j<- sample(1:(m-1),1)
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
    i <- sample(1:n, 1)  # random row
    j<- sample(1:(m-1),1) # random column
    while (!pass) {
      if (Dold[i,j] %in% block & Dold[i,j+1] %in% block){
        if (which(block==Dold[i,j],arr.ind = TRUE)[1]==which(block==Dold[i,j+1],arr.ind = TRUE)[1]){
          pass <- TRUE
        }
      }
      else {pass<-TRUE}
      i <- sample(1:n, 1)
      j<- sample(1:(m-1),1)
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

# m <- 7  # number of components
# n <- 60 # target sample size
# #t <- 2  # number of constraint blocks, only useful if automatic block used
# 
# 
# Df <- generate_full_design(m)
# Xf <- design_to_PWO(Df)
# block <- matrix(c(2,4,3,0),byrow = TRUE,nrow = 2)
# cons <- cons_pair(block)
# filter <- PWO_filter(cons,Xf)
# Df <- Df[-filter$cons_rows,]
# Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
# Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
# I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
# 
# set.seed(1234)
# start_found <- FALSE
# # we need to find a starting design with n rows that is not singular
# while(!start_found){
#   D_start <- Df[sample(c(1:nrow(Df)),n),]
#   X_start <- design_to_PWO(D_start)
#   X_start <- X_start[,-filter$cons_cols]
#   M_start <- t(X_start)%*%X_start
#   if(!is.singular.matrix(M_start)){
#     start_found <- TRUE
#     print("Found initial matrix")
#   }
# }
# 
# 
# D_out <- random_exchange(D_start, Bf, block)
# X_out <- design_to_PWO(D_out)
# X_out <- X_out[,-filter$cons_cols]
# I_out <- I_efficiency(X_out, Bf)
# print(I_full/I_out) # this is the relative efficiency


### test for multiple n when m is any positive integer, default to be 5.
### block is default to be matrix(c(2,4,3,0),byrow = TRUE,nrow = 2)

#set.seed(1234)

# The value of tot_iter, m, nlist needs to be input by hand
#tot_iter = 100
#m = 7
#nlist = c(60,90,120)
#block = matrix(c(2,4,3,0),byrow = TRUE,nrow = 2)

###################################################################
###A function to test for m and multiple n for random exchange###
###################################################################

PWO_test = function(tot_iter,m,nlist,block){
  cons = cons_pair(block)
  I_list = matrix(NA,nrow = tot_iter, ncol = length(nlist))
  
  I_name = rep(0,length(nlist))
  for (j in 1:length(nlist)) {
    I_name[j] <- c(paste("m=", m, ",n=", nlist[j], sep = ""))
  }
  
  
  colnames(I_list) = I_name
  
  for (n in 1:length(nlist)) {
    Df <- generate_full_design(m)
    Xf <- design_to_PWO(Df)
    filter <- PWO_filter(cons, Xf)
    Df <- Df[-filter$cons_rows,]
    Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
    Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
    I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
    for (iter in 1:tot_iter) {
      start_found <- FALSE
      # we need to find a starting design with n rows that is not singular
      while(!start_found){
        D_start <- Df[sample(1:nrow(Df),nlist[n]),]
        X_start <- design_to_PWO(D_start)
        X_start <- X_start[,-filter$cons_cols]
        M_start <- t(X_start)%*%X_start
        if(!is.singular.matrix(M_start)){
          start_found <- TRUE
          print("Found initial matrix")
        }
      }
      D_out <- random_exchange(D_start, Bf, block = block)
      X_out <- design_to_PWO(D_out)
      X_out <- X_out[,-filter$cons_cols]
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

#PWO_test(100,m=7,c(30,60,90),block = block)

#test_7 = PWO_test(100,m=7,c(30,60,90),block = block)


#####################################################################
###A function to test for m and multiple n for simulated annealing###
#####################################################################

sa_test = function(tot_iter = 100,m = 5,nlist = c(30,60,90),block = matrix(c(2,4,3,0),byrow = TRUE,nrow = 2)){
  cons = cons_pair(block)
  I_list = matrix(NA,nrow = tot_iter, ncol = length(nlist))
  
  I_name = rep(0,length(nlist))
  for (j in 1:length(nlist)) {
    I_name[j] <- c(paste("m=", m, ",n=", nlist[j], sep = ""))
  }
  
  
  colnames(I_list) = I_name
  
  for (n in 1:length(nlist)) {
    Df <- generate_full_design(m)
    Xf <- design_to_PWO(Df)
    filter <- PWO_filter(cons, Xf)
    Df <- Df[-filter$cons_rows,]
    Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
    Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
    I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
    for (iter in 1:tot_iter) {
      start_found <- FALSE
      # we need to find a starting design with n rows that is not singular
      while(!start_found){
        D_start <- Df[sample(1:nrow(Df),nlist[n]),]
        X_start <- design_to_PWO(D_start)
        X_start <- X_start[,-filter$cons_cols]
        M_start <- t(X_start)%*%X_start
        if(!is.singular.matrix(M_start)){
          start_found <- TRUE
          print("Found initial matrix")
        }
      }
      D_out <- simulated_annealing(D_start, Bf, block = block)
      X_out <- design_to_PWO(D_out)
      X_out <- X_out[,-filter$cons_cols]
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

#sa_test(100,5,c(30,60,90),block = block)

#test_7 = sa_test(100,7,c(30,60,90),block = block)


#--------------------------------------------
### Comparison between PWO and Simulated Annealing
#--------------------------------------------
set.seed(1234)
tot_iter = 100

m = 7
nlist = c(60,90,120)
block = matrix(c(1,2,4,3,5,6,7,0,0,0,0,0),byrow = TRUE,nrow = 2)

I_list_PWO = matrix(NA,nrow = tot_iter, ncol = length(nlist))
I_list_annealing = matrix(NA,nrow = tot_iter, ncol = length(nlist))

I_name = rep(0,length(nlist))
for (j in 1:length(nlist)) {
  I_name[j] <- c(paste("m=", m, ",n=", nlist[j], sep = ""))
}


colnames(I_list_PWO) = I_name
colnames(I_list_annealing) = I_name

for (n in 1:length(nlist)) {
  Df <- generate_full_design(m)
  Xf <- design_to_PWO(Df)
  cons <- cons_pair(block)
  filter <- PWO_filter(cons, Xf)
  Df <- Df[-filter$cons_rows,]
  Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
  Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
  I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
  for (iter in 1:tot_iter) {
    start_found <- FALSE
    # we need to find a starting design with n rows that is not singular
    while(!start_found){
      D_start <- Df[sample(1:nrow(Df),nlist[n]),]
      X_start <- design_to_PWO(D_start)
      X_start <- X_start[,-filter$cons_cols]
      M_start <- t(X_start)%*%X_start
      if(!is.singular.matrix(M_start)){
        start_found <- TRUE
        print("Found initial matrix")
      }
    }
    D_out_PWO <- random_exchange(D_start, Bf, block = block)
    X_out_PWO <- design_to_PWO(D_out_PWO)
    X_out_PWO <- X_out_PWO[,-filter$cons_cols]
    I_out_PWO <- I_efficiency(X_out_PWO, Bf)
    I_list_PWO[iter,n]=I_full/I_out_PWO
    D_out_annealing <- simulated_annealing(D_start, Bf, block = block)
    X_out_annealing <- design_to_PWO(D_out_annealing)
    X_out_annealing <- X_out_annealing[,-filter$cons_cols]
    I_out_annealing <- I_efficiency(X_out_annealing, Bf)
    I_list_annealing[iter,n]=I_full/I_out_annealing
  }
}



# #calculate mean for each n
# matrix(apply(I_list_PWO, 2, mean),ncol = length(nlist),byrow = FALSE)
# matrix(apply(I_list_annealing, 2, mean),ncol = length(nlist),byrow = FALSE)
# 
# #Box plots
# rownames(I_mean_mat)=c(paste("m=", m, sep = ""))
# colnames(I_mean_mat)=c(paste("n=", nlist, sep = ""))
# 
# preplot = cbind(matrix(I_list_PWO,ncol = 1),matrix(I_list_annealing,ncol = 1))
# 
# preplot = as.data.frame(cbind(preplot,rep(nlist,each=tot_iter)))
# 
# colnames(preplot)=c("RandomExchange","AnnealingSimulation","n")
# preplot$n=as.factor(preplot$n)
# 
# p1 <- ggplot(preplot, aes(x=n, y=RandomExchange)) + 
#   geom_boxplot()+ggtitle("Box plot for different subset sizes") +
#   ylab("I-efficiency")
# p1
# p2 <- ggplot(preplot, aes(x=n, y=AnnealingSimulation)) + 
#   geom_boxplot()+ggtitle("Box plot for different subset sizes") +
#   ylab("I-efficiency")
# p2
# preplot_long = melt(preplot,id.vars = c("n"))
# p <- ggplot(preplot_long, aes(x=n,y=value,fill=n))+
#   geom_boxplot() + labs(title=paste("Algorithm Efficiency Comparison When m = ",m))
#                     + facet_wrap(~variable) +  ylab("Relative I-efficiency")
# p


#block complete balanced design
#5 treatments, 1 car, 4 tyre
#one treatment each tyre
#can't fit the 5th 
#block1 1-4, block2 1235, block3 1245, block4 1345, block5 2345
#pairwise balance
#Is full design optimal if constraints from BIBD?

#m=7, n=30 failed of singular

# # of blocks b
# # of components m
# sample size n
# # unconstrainted components m_u
# for m = 6,7,8,9
# mu = 0,1,2,3; b = 2, 3
# arrange m-mu components into b blocks balanced/imbalanced(1,mc-1)
# n = 1+2choose(p,2),1+3choose(p,2)+... < nrow(full constrained design)
# compare balance/imbalanced

# 2^k factorial design
# k factors; each has 2 levels; each has 2 levels (-1,+1)
# response: y
# full design: all combinations (2^k) of the levels 

# 4*3*3*2*2




