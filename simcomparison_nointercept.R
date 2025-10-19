# Try simulating 1 dataset and compare the MSPE, AIC, BIC, and show the parameter estimates for TE1
# Please run TE_PWO_CP_FP_modelsource.R before running this file.

# setwd("C:/Users/riosn/Downloads/Example51_nointercept")

source("TE_PWO_CP_FP_modelsource.R")

set.seed(1234)

# OofA_transition_model_matrix
# Input: D, an n by m matrix of permutations
# Output: X, an n by 2*(m choose 2) matrix of model parameters
OofA_transition_model_matrix <- function(D){
  
  if(is.vector(D)){
    
    m <- length(D)
    # X2 <- matrix(0, nrow = 1, ncol = m*(m-1))
    X2 <- numeric(m*(m-1))
    x2namevec <- c()
    counter <- 1
    for(i in 1:m){
      for(j in 1:m){
        if(i != j){
          x2namevec <- c(x2namevec, paste("x", i, "|", j, sep = ""))
        }
      }
    }
    
    for(j in 2:m){
      
      prev <- D[j-1]
      curr <- D[j]
      col_ind <- which(x2namevec == paste("x", curr, "|", prev, sep = ""))
      X2[col_ind] <- 1
      
    }
    #X2 <- c(1,X2) # add the intercept
    names(X2) <- x2namevec
    return(X2)
  }
  n <- nrow(D)
  m <- ncol(D)
  X2 <- matrix(0, nrow = n, ncol = m*(m-1))
  x2namevec <- c()
  counter <- 1
  for(i in 1:m){
    for(j in 1:m){
      if(i != j){
        x2namevec <- c(x2namevec, paste("x", i, "|", j, sep = ""))
      }
    }
  }
  
  for(i in 1:n){
    for(j in 2:m){
      prev <- D[i,j-1]
      curr <- D[i,j]
      col_ind <- which(x2namevec == paste("x", curr, "|", prev, sep = ""))
      X2[i,col_ind] <- 1
    }
  }
  #X <- cbind(1,X2)
  colnames(X2) <- x2namevec
  #X <- X[,-ncol(X)]
  return(X2)
}

SD <- read.table("SD.txt")
ExpectedTerminalTravelTime <- read.table("ExpectedTerminalTravelTime.txt")

m=7
D_full <- generate_full_design(m)
block = matrix(c(1,rep(0,5),seq(2,7)),nrow=2,byrow = TRUE)
cons <-cons_pair(block)

Xf = OofA_transition_model_matrix(D_full)
filter <- TE_filter1(cons=cons, D=D_full, Z=Xf)
D_full <- D_full[-filter$cons_rows,]
Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
B = t(Xf)%*%Xf/nrow(Xf)
id <- apply(D_full,1, function(x) paste(x,sep='',collapse=''))

n = 200
n_iter = 100

AIC_BB_mat = matrix(nrow = n_iter, ncol = 5)
AIC_GRASP_mat = matrix(nrow = n_iter, ncol = 5)
MSPE_BB_mat = matrix(nrow = n_iter, ncol = 5)
MSPE_GRASP_mat = matrix(nrow = n_iter, ncol = 5)

# # BB
# Dtest = D_full[sample(nrow(D_full),n),]
# D_out_BB = bubble_sort_design_B1(init_runs = Dtest,B,criterion="I",block=block)
# 
# # GRASP
# D_out_GRASP = grasp_v3_B(D = Dtest,B,block=block,criterion="I", max_iter = 10, n_rand = 250)

for(iter in 1:n_iter){
  
  
  # Generate Y for each iteration
  Y = rep(0,nrow(D_full))
  for (i in 1:nrow(D_full)) {
    t = rnorm(1, mean = ExpectedTerminalTravelTime[D_full[i,1],D_full[i,2]-1], sd = SD[D_full[i,1],D_full[i,2]-1])
    for (j in (2:m-1)) {
      t = t + rnorm(1, mean = ExpectedTerminalTravelTime[D_full[i,j],D_full[i,j+1]-1], sd = SD[D_full[i,j],D_full[i,j+1]-1])
    }
    Y[i] = t
  }
  rank = sort(Y,decreasing = FALSE)
  
  # BB
  Dtest = D_full[sample(nrow(D_full),n),]
  D_out_BB = bubble_sort_design_B1(init_runs = Dtest,B,criterion="I",block=block)

  # GRASP
  D_out_GRASP = grasp_v3_B(D = Dtest,B,block=block,criterion="I", max_iter = 10, n_rand = 250)
  
  # TE1*BB
  
  Xf = OofA_transition_model_matrix(D_full)
  filter <- TE_filter1(cons=cons, D=D_full, Z=Xf)
  Xf <- Xf[,-filter$cons_cols]
  case5_cons = as.data.frame(Xf)
  
  TE_out_BB <- OofA_transition_model_matrix(D_out_BB)
  TE_out_BB <- TE_out_BB[,-filter$cons_cols] 
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  }
  #TE_out_BB = TE_out_BB[,-1]
  esc_BB = as.data.frame(cbind(TE_out_BB,Y_out_BB))
  cons_BB_TE_model = lm(Y_out_BB ~-1+., data = esc_BB)
  BB_TE_pred = predict(cons_BB_TE_model,newdata = case5_cons)
  rank_BB_TE = which(rank == min(Y[which(BB_TE_pred==min(BB_TE_pred))]))[1]
  
  
  # TE1*GRASP
  TE_out_GRASP <- OofA_transition_model_matrix(D_out_GRASP)
  TE_out_GRASP <- TE_out_GRASP[,-filter$cons_cols] 
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  #TE_out_GRASP = TE_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(TE_out_GRASP,Y_out_GRASP))
  cons_GRASP_TE_model = lm(Y_out_GRASP ~-1+., data = esc_GRASP)
  GRASP_TE_pred = predict(cons_GRASP_TE_model,newdata = case5_cons)
  rank_GRASP_TE = which(rank == min(Y[which(GRASP_TE_pred==min(GRASP_TE_pred))]))[1]
  
  
  # TE2*BB
  
  Xf = OofA_transition_model_matrix_2(D_full)
  filter <- TE_filter1(cons=cons, D=D_full, Z=Xf)
  Xf <- Xf[,-filter$cons_cols]
  case5_cons = as.data.frame(Xf)
  
  TE2_out_BB <- OofA_transition_model_matrix_2(D_out_BB)
  TE2_out_BB <- TE2_out_BB[,-filter$cons_cols] 
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  }
  #TE2_out_BB = TE2_out_BB[,-1]
  esc_BB = as.data.frame(cbind(TE2_out_BB,Y_out_BB))
  cons_BB_TE2_model = lm(Y_out_BB ~-1+., data = esc_BB)
  BB_TE2_pred = predict(cons_BB_TE2_model,newdata = case5_cons)
  rank_BB_TE2 = which(rank == min(Y[which(BB_TE2_pred==min(BB_TE2_pred))]))[1]
  
  
  # TE2*GRASP
  TE2_out_GRASP <- OofA_transition_model_matrix_2(D_out_GRASP)
  TE2_out_GRASP <- TE2_out_GRASP[,-filter$cons_cols] 
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  #TE2_out_GRASP = TE2_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(TE2_out_GRASP,Y_out_GRASP))
  cons_GRASP_TE2_model = lm(Y_out_GRASP ~-1+., data = esc_GRASP)
  GRASP_TE2_pred = predict(cons_GRASP_TE2_model,newdata = case5_cons)
  rank_GRASP_TE2 = which(rank == min(Y[which(GRASP_TE2_pred==min(GRASP_TE2_pred))]))[1]
  
  
  # PWO*BB
  Xf = design_to_PWO(D_full)
  filter <- PWO_filter(cons=cons, Z=Xf)
  Xf <- Xf[,-filter$cons_cols]
  case5_cons = as.data.frame(Xf)
  
  PWO_out_BB <- design_to_PWO(D_out_BB)
  PWO_out_BB <- PWO_out_BB[,-filter$cons_cols] 
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separaPWO regressed and predicPWOd data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
    #Y_out_BB[i] = Y[which(apply(D_out_BB, 1, function(x) identical(id,pasPWO(x,sep='',collapse=''))))]
  }
  #PWO_out_BB = PWO_out_BB[,-1]
  esc_BB = as.data.frame(cbind(PWO_out_BB,Y_out_BB))
  cons_BB_PWO_model = lm(Y_out_BB ~., data = esc_BB)
  BB_PWO_pred = predict(cons_BB_PWO_model,newdata = case5_cons)
  rank_BB_PWO = which(rank == min(Y[which(BB_PWO_pred==min(BB_PWO_pred))]))[1]
  
  # PWO*GRASP
  PWO_out_GRASP <- design_to_PWO(D_out_GRASP)
  PWO_out_GRASP <- PWO_out_GRASP[,-filter$cons_cols] 
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separaPWO regressed and predicPWOd data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
    #Y_out_GRASP[i] = Y[which(apply(D_out_GRASP, 1, function(x) identical(id,pasPWO(x,sep='',collapse=''))))]
  }
  #PWO_out_GRASP = PWO_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(PWO_out_GRASP,Y_out_GRASP))
  cons_GRASP_PWO_model = lm(Y_out_GRASP ~., data = esc_GRASP)
  GRASP_PWO_pred = predict(cons_GRASP_PWO_model,newdata = case5_cons)
  rank_GRASP_PWO = which(rank == min(Y[which(GRASP_PWO_pred==min(GRASP_PWO_pred))]))[1]
  
  # CP*BB
  case5_cons = as.data.frame(design_to_CP(D_full))
  CP_out_BB <- design_to_CP(D_out_BB)
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  }
  
  esc_BB = as.data.frame(cbind(CP_out_BB,Y_out_BB))
  
  cons_BB_CP_model = lm(Y_out_BB ~., data = esc_BB)
  BB_CP_pred = predict(cons_BB_CP_model,newdata = case5_cons)
  rank_BB_CP = which(rank == min(Y[which(BB_CP_pred==min(BB_CP_pred))]))[1]
  
  
  # CP*GRASP
  case5_cons = as.data.frame(design_to_CP(D_full))
  CP_out_GRASP <- design_to_CP(D_out_GRASP)
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  
  esc_GRASP = as.data.frame(cbind(CP_out_GRASP,Y_out_GRASP))
  cons_GRASP_CP_model = lm(Y_out_GRASP ~., data = esc_GRASP)
  GRASP_CP_pred = predict(cons_GRASP_CP_model,newdata = case5_cons)
  rank_GRASP_CP = which(rank == min(Y[which(GRASP_CP_pred==min(GRASP_CP_pred))]))[1]
  
  
  # FP*BB
  case5_cons = as.data.frame(design_to_FP(D_full))
  Xf = OofA_transition_model_matrix_2(D_full)
  filter <- TE_filter1(cons=cons, D=D_full, Z=Xf)
  Xf <- Xf[,-filter$cons_cols]
  
  FP_out_BB <- design_to_FP(D_out_BB)
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  }
  #FP_out_BB = FP_out_BB[,-1]
  esc_BB = as.data.frame(cbind(FP_out_BB,Y_out_BB))
  cons_BB_FP_model = lm(Y_out_BB ~., data = esc_BB)
  BB_FP_pred = predict(cons_BB_FP_model,newdata = case5_cons)
  rank_BB_FP = which(rank == min(Y[which(BB_FP_pred==min(BB_FP_pred))]))[1]
  
  
  # FP*GRASP
  FP_out_GRASP <- design_to_FP(D_out_GRASP)
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  #FP_out_GRASP = FP_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(FP_out_GRASP,Y_out_GRASP))
  cons_GRASP_FP_model = lm(Y_out_GRASP ~., data = esc_GRASP)
  GRASP_FP_pred = predict(cons_GRASP_FP_model,newdata = case5_cons)
  rank_GRASP_FP = which(rank == min(Y[which(GRASP_FP_pred==min(GRASP_FP_pred))]))[1]
  
  
  AIC_BB_mat[iter,] = c(AIC(cons_BB_PWO_model),AIC(cons_BB_TE_model),AIC(cons_BB_TE2_model),AIC(cons_BB_CP_model),AIC(cons_BB_FP_model))
  AIC_GRASP_mat[iter,] = c(AIC(cons_GRASP_PWO_model),AIC(cons_GRASP_TE_model),AIC(cons_GRASP_TE2_model),AIC(cons_GRASP_CP_model),AIC(cons_GRASP_FP_model))
  MSPE_BB_mat[iter,] = c(mean((Y-BB_PWO_pred)^2),mean((Y-BB_TE_pred)^2),mean((Y-BB_TE2_pred)^2),mean((Y-BB_CP_pred)^2),mean((Y-BB_FP_pred)^2))
  MSPE_GRASP_mat[iter,] = c(mean((Y-GRASP_PWO_pred)^2),mean((Y-GRASP_TE_pred)^2),mean((Y-GRASP_TE2_pred)^2),mean((Y-GRASP_CP_pred)^2),mean((Y-GRASP_FP_pred)^2))
  
}


res_table = cbind(round(colMeans(AIC_BB_mat),2), round(colMeans(AIC_GRASP_mat),2), round(colMeans(MSPE_BB_mat),2), round(colMeans(MSPE_GRASP_mat),2))
colnames(res_table) = c("AIC_BB","AIC_GRASP","MSPE_BB","MSPE_GRASP")
res_table

