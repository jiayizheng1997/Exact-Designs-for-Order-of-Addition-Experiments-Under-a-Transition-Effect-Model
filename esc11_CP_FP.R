# This file is a parallel of esc11_PWO_TE.R
# esc11_PWO_TE.R was named esc11.R when we wrote the first version of the article.



#setwd("C:/Users/riosn/OneDrive/Desktop")
set.seed(1234)

library(igraph)
esc11_raw <- read.table("ESC11.txt", header = F)

D_full <- as.matrix(read.csv("esc11_mod_full_design.csv")[,-1])

# returns the total cost for a single permutation of 1:m
get_cost <- function(perm, cost_matrix){
  m <- length(perm)
  cost <- 0
  for(i in 1:(m-1)){
    cost <- cost + cost_matrix[perm[i],perm[i+1]]
  }
  return(cost)
}

esc11 <- esc11_raw[-c(1,13),-c(1,13)]

# Modified so that:
# Group 1: {1,4,7,9,10}
# Group 2: {2,6,8}
# Group 3: {3,5,11}

esc11[2,1] <- esc11[2,4] <- esc11[2,7] <- esc11[2,9] <- esc11[2,10] <- -1
esc11[6,1] <- esc11[6,4] <- esc11[6,7] <- esc11[6,9] <- esc11[6,10] <- -1
esc11[8,1] <- esc11[8,4] <- esc11[8,7] <- esc11[8,9] <- esc11[8,10] <- -1
esc11[3,1] <- esc11[3,4] <- esc11[3,7] <- esc11[3,9] <- esc11[3,10] <- esc11[3,2] <- esc11[3,6] <- esc11[3,8] <- -1
esc11[5,1] <- esc11[5,4] <- esc11[5,7] <- esc11[5,9] <- esc11[5,10] <- esc11[5,2] <- esc11[5,6] <- esc11[5,8] <- -1
esc11[11,1] <- esc11[11,4] <- esc11[11,7] <- esc11[11,9] <- esc11[11,10] <- esc11[11,2] <- esc11[11,6] <- esc11[11,8] <- -1




adjmat <- matrix(0, nrow = nrow(esc11), ncol = ncol(esc11))
for(i in 1:ncol(esc11)){
  for(j in 1:ncol(esc11)){
    if(esc11[i,j] == -1){
      adjmat[j,i] = 1
    }
  }
}


cost_matrix <- pmax(esc11,0)

G <- graph_from_adjacency_matrix(adjmat, mode = "directed")

layout_mat = cbind(c(1,2,3,1,3,2,1,2,1,1,3),c(1,1,1,2,2,2,3,3,4,5,3))
plot(G, layout = layout_mat)

Y <- apply(D_full, 1, get_cost, cost_matrix = cost_matrix)  # these are all of the true responses

rank = sort(Y,decreasing = FALSE)

min_cost <- min(Y)  # 2330
opt_order <- D_full[which.min(Y),]  # 10,7,4,1,9,2,6,8,5,3,11 is the optimal order

##################################################
# Two by two comparison of CP, FC
# and bubble sorting, GRASP.
##################################################

n = 200

block = matrix(c(1,4,7,9,10,2,6,8,0,0,3,5,11,0,0),nrow=3,byrow = TRUE)
Xf = OofA_transition_model_matrix(D_full)
cons <-cons_pair(block)
filter <- TE_filter1(cons=cons, D=D_full, Z=Xf)
id <- apply(D_full,1, function(x) paste(x,sep='',collapse=''))
Y <- apply(D_full, 1, get_cost, cost_matrix = cost_matrix)


Xf <- Xf[,-filter$cons_cols]
B = t(Xf)%*%Xf/nrow(Xf)

result1 = data.frame(matrix(nrow=100,ncol = 4))


colnames(result1) = c('CP*BB', 'CP*GRASP','FP*BB', 'FP*GRASP')


N = 100 # the number of runs we want to run for each combination

for (iter in 1:N) {  # make 100 runs for each combination
  
  # BB
  Dtest = D_full[sample(nrow(D_full),n),]
  D_out_BB = bubble_sort_design_B1(init_runs = Dtest,B,criterion="I",block=block)
  
  # GRASP
  D_out_GRASP = grasp_v3_B(D = Dtest,B,block=block,criterion="I", max_iter = 10, n_rand = 250)
  
  # CP*BB
  case5_cons = as.data.frame(design_to_CP(D_full))
  CP_out_BB <- design_to_CP(D_out_BB)
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  }
  
  esc_BB = as.data.frame(cbind(CP_out_BB,Y_out_BB))
  
  # Stepwise selection, as CP model cannot be filtered
  # full_model <- lm(Y_out_BB ~ ., data = esc_BB)
  # null_model <- lm(Y_out_BB ~ 1, data = esc_BB) # Only intercept
  # cons_BB_CP_model <- step(null_model, 
  #                       scope = list(lower = null_model, upper = full_model), 
  #                       direction = "both", 
  #                       trace = 1) # Set to 1 for step-by-step output
  
  cons_BB_CP_model = lm(Y_out_BB ~ ., data = esc_BB)
  BB_CP_pred = predict(cons_BB_CP_model,newdata = case5_cons)
  rank_BB_CP = which(rank == min(Y[which(BB_CP_pred==min(BB_CP_pred))]))[1]
  result1[iter,1] = rank_BB_CP
  
  # CP*GRASP
  case5_cons = as.data.frame(design_to_CP(D_full))
  CP_out_GRASP <- design_to_CP(D_out_GRASP)
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  
  esc_GRASP = as.data.frame(cbind(CP_out_GRASP,Y_out_GRASP))
  cons_GRASP_CP_model = lm(Y_out_GRASP ~ ., data = esc_GRASP)
  GRASP_CP_pred = predict(cons_GRASP_CP_model,newdata = case5_cons)
  rank_GRASP_CP = which(rank == min(Y[which(GRASP_CP_pred==min(GRASP_CP_pred))]))[1]
  result1[iter,2] = rank_GRASP_CP
  
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
  FP_out_BB = FP_out_BB[,-1]
  esc_BB = as.data.frame(cbind(FP_out_BB,Y_out_BB))
  cons_BB_FP_model = lm(Y_out_BB ~ ., data = esc_BB)
  BB_FP_pred = predict(cons_BB_FP_model,newdata = case5_cons)
  rank_BB_FP = which(rank == min(Y[which(BB_FP_pred==min(BB_FP_pred))]))[1]
  result1[iter,3] = rank_BB_FP
  
  # FP*GRASP
  FP_out_GRASP <- design_to_FP(D_out_GRASP)
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  FP_out_GRASP = FP_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(FP_out_GRASP,Y_out_GRASP))
  cons_GRASP_FP_model = lm(Y_out_GRASP ~ ., data = esc_GRASP)
  GRASP_FP_pred = predict(cons_GRASP_FP_model,newdata = case5_cons)
  rank_GRASP_FP = which(rank == min(Y[which(GRASP_FP_pred==min(GRASP_FP_pred))]))[1]
  result1[iter,4] = rank_GRASP_FP
  
}


n = 150

block = matrix(c(1,4,7,9,10,2,6,8,0,0,3,5,11,0,0),nrow=3,byrow = TRUE)
Xf = OofA_transition_model_matrix(D_full)
cons <-cons_pair(block)
filter <- TE_filter1(cons=cons, D=D_full, Z=Xf)
id <- apply(D_full,1, function(x) paste(x,sep='',collapse=''))
Y <- apply(D_full, 1, get_cost, cost_matrix = cost_matrix)


Xf <- Xf[,-filter$cons_cols]
B = t(Xf)%*%Xf/nrow(Xf)

result2 = data.frame(matrix(nrow=100,ncol = 4))


colnames(result2) = c('CP*BB', 'CP*GRASP','FP*BB', 'FP*GRASP')


N = 100 # the number of runs we want to run for each combination

for (iter in 1:N) {  # make 100 runs for each combination
  
  # BB
  Dtest = D_full[sample(nrow(D_full),n),]
  D_out_BB = bubble_sort_design_B1(init_runs = Dtest,B,criterion="I",block=block)
  
  # GRASP
  D_out_GRASP = grasp_v3_B(D = Dtest,B,block=block,criterion="I", max_iter = 10, n_rand = 250)
  
  # CP*BB
  case5_cons = as.data.frame(design_to_CP(D_full))
  CP_out_BB <- design_to_CP(D_out_BB)
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  }
  
  esc_BB = as.data.frame(cbind(CP_out_BB,Y_out_BB))
  
  cons_BB_CP_model = lm(Y_out_BB ~ ., data = esc_BB)
  BB_CP_pred = predict(cons_BB_CP_model,newdata = case5_cons)
  rank_BB_CP = which(rank == min(Y[which(BB_CP_pred==min(BB_CP_pred))]))[1]
  result2[iter,1] = rank_BB_CP
  
  # CP*GRASP
  case5_cons = as.data.frame(design_to_CP(D_full))
  CP_out_GRASP <- design_to_CP(D_out_GRASP)
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  
  esc_GRASP = as.data.frame(cbind(CP_out_GRASP,Y_out_GRASP))
  cons_GRASP_CP_model = lm(Y_out_GRASP ~ ., data = esc_GRASP)
  GRASP_CP_pred = predict(cons_GRASP_CP_model,newdata = case5_cons)
  rank_GRASP_CP = which(rank == min(Y[which(GRASP_CP_pred==min(GRASP_CP_pred))]))[1]
  result2[iter,2] = rank_GRASP_CP
  
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
  FP_out_BB = FP_out_BB[,-1]
  esc_BB = as.data.frame(cbind(FP_out_BB,Y_out_BB))
  cons_BB_FP_model = lm(Y_out_BB ~ ., data = esc_BB)
  BB_FP_pred = predict(cons_BB_FP_model,newdata = case5_cons)
  rank_BB_FP = which(rank == min(Y[which(BB_FP_pred==min(BB_FP_pred))]))[1]
  result2[iter,3] = rank_BB_FP
  
  # FP*GRASP
  FP_out_GRASP <- design_to_FP(D_out_GRASP)
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  FP_out_GRASP = FP_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(FP_out_GRASP,Y_out_GRASP))
  cons_GRASP_FP_model = lm(Y_out_GRASP ~ ., data = esc_GRASP)
  GRASP_FP_pred = predict(cons_GRASP_FP_model,newdata = case5_cons)
  rank_GRASP_FP = which(rank == min(Y[which(GRASP_FP_pred==min(GRASP_FP_pred))]))[1]
  result2[iter,4] = rank_GRASP_FP
  
}


colMeans(result1)
colMeans(result2)









