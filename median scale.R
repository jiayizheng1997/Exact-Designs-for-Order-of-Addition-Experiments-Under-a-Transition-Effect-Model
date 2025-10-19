library(combinat)
library(matrixcalc)
library(doParallel)
library(foreach)
library(AlgDesign)
library(ggplot2)
library(reshape2)
library(retry)
library(readr)
library(xtable)

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

# m = 6
data = read.csv('df1_out10000max_iter.csv')
data[data<0.00001] <- NA


m1=as.data.frame(lapply(data[1:100,], median, na.rm = TRUE)) #remove the na and small values
m2=as.data.frame(lapply(data[101:200,], median, na.rm = TRUE))
m3=as.data.frame(lapply(data[201:300,], median, na.rm = TRUE))
mlist = rbind(m1,m2,m3)
mlist = as.data.frame(mlist)

for (i in 1:12) {
  data_row = DF[i*4,]
  b = data_row[2] #number of blocks
  mu = data_row[3] #number of free components
  e = data_row[4] #extremeness
  colnames(mlist)[12*i-9] = paste(c(colnames(mlist)[12*i-9],'b',b,'mu',mu,'e',e),sep='',collapse='')
  colnames(mlist)[12*i-6] = paste(c(colnames(mlist)[12*i-6],'b',b,'mu',mu,'e',e),sep='',collapse='')
  colnames(mlist)[12*i-3] = paste(c(colnames(mlist)[12*i-3],'b',b,'mu',mu,'e',e),sep='',collapse='')
  colnames(mlist)[12*i] = paste(c(colnames(mlist)[12*i],'b',b,'mu',mu,'e',e),sep='',collapse='')
}

write.csv(mlist,'scaled median for m=6.csv')


# m = 7
data = read.csv('df2_out10000max_iter.csv')
data[data<0.00001] <- NA


m1=as.data.frame(lapply(data[1:100,], median, na.rm = TRUE)) #remove the na and small values
m2=as.data.frame(lapply(data[101:200,], median, na.rm = TRUE))
m3=as.data.frame(lapply(data[201:300,], median, na.rm = TRUE))
mlist = rbind(m1,m2,m3)
mlist = as.data.frame(mlist)


for (i in 1:12) {
  data_row = DF[i*4,]
  b = data_row[2] #number of blocks
  mu = data_row[3] #number of free components
  e = data_row[4] #extremeness
  colnames(mlist)[12*i-9] = paste(c(colnames(mlist)[12*i-9],'b',b,'mu',mu,'e',e),sep='',collapse='')
  colnames(mlist)[12*i-6] = paste(c(colnames(mlist)[12*i-6],'b',b,'mu',mu,'e',e),sep='',collapse='')
  colnames(mlist)[12*i-3] = paste(c(colnames(mlist)[12*i-3],'b',b,'mu',mu,'e',e),sep='',collapse='')
  colnames(mlist)[12*i] = paste(c(colnames(mlist)[12*i],'b',b,'mu',mu,'e',e),sep='',collapse='')
}

write.csv(mlist,'scaled median for m=7.csv')

# m = 8
data = read.csv('df3_out10000max_iter.csv')
data[data<0.00001] <- NA


m1=as.data.frame(lapply(data[1:100,], median, na.rm = TRUE)) #remove the na and small values
m2=as.data.frame(lapply(data[101:200,], median, na.rm = TRUE))
m3=as.data.frame(lapply(data[201:300,], median, na.rm = TRUE))
mlist = rbind(m1,m2,m3)
mlist = as.data.frame(mlist)



for (i in 1:12) {
  data_row = DF[i*4,]
  b = data_row[2] #number of blocks
  mu = data_row[3] #number of free components
  e = data_row[4] #extremeness
  colnames(mlist)[12*i-9] = paste(c(colnames(mlist)[12*i-9],'b',b,'mu',mu,'e',e),sep='',collapse='')
  colnames(mlist)[12*i-6] = paste(c(colnames(mlist)[12*i-6],'b',b,'mu',mu,'e',e),sep='',collapse='')
  colnames(mlist)[12*i-3] = paste(c(colnames(mlist)[12*i-3],'b',b,'mu',mu,'e',e),sep='',collapse='')
  colnames(mlist)[12*i] = paste(c(colnames(mlist)[12*i],'b',b,'mu',mu,'e',e),sep='',collapse='')
}

write.csv(mlist,'scaled median for m=8.csv')


s <- read_csv("scaled median for m=6.csv")
data = rbind(cbind(s$SimulatedAnnealingb2mu0e0,s$SA_votingb2mu0e0,s$optFederovb2mu0e0,s$InsertionSortb2mu0e0),cbind(s$SimulatedAnnealing.8b3mu0e0,s$SA_voting.8b3mu0e0,s$optFederov.8b3mu0e0,s$SimulatedAnnealing.8b3mu0e0),cbind(s$SimulatedAnnealing.1b2mu0e1,s$SA_voting.1b2mu0e1,s$optFederov.1b2mu0e1,s$InsertionSort.1b2mu0e1))
print(xtable(data,digits = 3))
data = rbind(cbind(s$SimulatedAnnealing.2b2mu1e0,s$SA_voting.2b2mu1e0,s$optFederov.2b2mu1e0,s$InsertionSort.2b2mu1e0),cbind(s$SimulatedAnnealing.9b3mu1e0,s$SA_voting.9b3mu1e0,s$optFederov.9b3mu1e0,s$SimulatedAnnealing.9b3mu1e0),cbind(s$SimulatedAnnealing.3b2mu1e1,s$SA_voting.3b2mu1e1,s$optFederov.3b2mu1e1,s$InsertionSort.3b2mu1e1))
print(xtable(data,digits = 3))


s <- read_csv("scaled median for m=7.csv")
data = rbind(cbind(s$SimulatedAnnealingb2mu0e0,s$SA_votingb2mu0e0,s$optFederovb2mu0e0,s$InsertionSortb2mu0e0),cbind(s$SimulatedAnnealing.8b3mu0e0,s$SA_voting.8b3mu0e0,s$optFederov.8b3mu0e0,s$SimulatedAnnealing.8b3mu0e0),cbind(s$SimulatedAnnealing.1b2mu0e1,s$SA_voting.1b2mu0e1,s$optFederov.1b2mu0e1,s$InsertionSort.1b2mu0e1))
print(xtable(data,digits = 3))
data = rbind(cbind(s$SimulatedAnnealing.2b2mu1e0,s$SA_voting.2b2mu1e0,s$optFederov.2b2mu1e0,s$InsertionSort.2b2mu1e0),cbind(s$SimulatedAnnealing.9b3mu1e0,s$SA_voting.9b3mu1e0,s$optFederov.9b3mu1e0,s$SimulatedAnnealing.9b3mu1e0),cbind(s$SimulatedAnnealing.3b2mu1e1,s$SA_voting.3b2mu1e1,s$optFederov.3b2mu1e1,s$InsertionSort.3b2mu1e1))
print(xtable(data,digits = 3))

s <- read_csv("scaled median for m=8.csv")
data = rbind(cbind(s$SimulatedAnnealingb2mu0e0,s$SA_votingb2mu0e0,s$optFederovb2mu0e0,s$InsertionSortb2mu0e0),cbind(s$SimulatedAnnealing.8b3mu0e0,s$SA_voting.8b3mu0e0,s$optFederov.8b3mu0e0,s$SimulatedAnnealing.8b3mu0e0),cbind(s$SimulatedAnnealing.1b2mu0e1,s$SA_voting.1b2mu0e1,s$optFederov.1b2mu0e1,s$InsertionSort.1b2mu0e1))
print(xtable(data,digits = 3))
data = rbind(cbind(s$SimulatedAnnealing.2b2mu1e0,s$SA_voting.2b2mu1e0,s$optFederov.2b2mu1e0,s$InsertionSort.2b2mu1e0),cbind(s$SimulatedAnnealing.9b3mu1e0,s$SA_voting.9b3mu1e0,s$optFederov.9b3mu1e0,s$SimulatedAnnealing.9b3mu1e0),cbind(s$SimulatedAnnealing.3b2mu1e1,s$SA_voting.3b2mu1e1,s$optFederov.3b2mu1e1,s$InsertionSort.3b2mu1e1))
print(xtable(data,digits = 3))
