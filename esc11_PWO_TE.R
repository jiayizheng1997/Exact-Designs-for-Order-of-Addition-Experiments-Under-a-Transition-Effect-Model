# This file used to be esc11.R
# Name was changed due to addition of more algorithms.
# Nothing but name was changed.

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

min_cost <- min(Y)  # 2330
opt_order <- D_full[which.min(Y),]  # 10,7,4,1,9,2,6,8,5,3,11 is the optimal order

########################################
# Step 1: Find a design of size n = 200
########################################

n = 200

block = matrix(c(1,4,7,9,10,2,6,8,0,0,3,5,11,0,0),nrow=3,byrow = TRUE)

Xf = OofA_transition_model_matrix(D_full)
cons <-cons_pair(block)
filter <- TE_filter1(cons=cons, D=D_full, Z=Xf)

Xf <- Xf[,-filter$cons_cols]
B = t(Xf)%*%Xf/nrow(Xf)
Dtest = D_full[sample(nrow(D_full),n),]
#D_out_SA = simulated_annealing_B1(Dtest,B,max_iter = 100, block, temp_init = 1, criterion = "I")
D_out_BB = bubble_sort_design_B1(init_runs = Dtest,B,block=block)



########################################
# Step 2: Fit the TE and PWO models
########################################

id <- apply(D_full,1, function(x) paste(x,sep='',collapse=''))
Y <- apply(D_full, 1, get_cost, cost_matrix = cost_matrix)  # these are all of the true responses
cons <-cons_pair(block)

min_cost <- min(Y)  # 2330
opt_order <- D_full[which.min(Y),]  # 10,7,4,1,9,2,6,8,5,3,11 is the optimal order

rank = sort(Y,decreasing = FALSE)


### TE model

Xf = OofA_transition_model_matrix(D_full)
filter <- TE_filter1(cons=cons, D=D_full, Z=Xf)
Xf <- Xf[,-filter$cons_cols]
case5_cons = as.data.frame(Xf[,-1])
TE_out_BB <- OofA_transition_model_matrix(D_out_BB)
TE_out_BB <- TE_out_BB[,-filter$cons_cols] 
Y_out_BB <- rep(0,nrow(D_out_BB))
# separate regressed and predicted data
for (i in 1:nrow(D_out_BB)) {
  Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  #Y_out_BB[i] = Y[which(apply(D_out_BB, 1, function(x) identical(id,paste(x,sep='',collapse=''))))]
}
TE_out_BB = TE_out_BB[,-1]
esc_BB = as.data.frame(cbind(TE_out_BB,Y_out_BB))
cons_BB_TE_model = lm(Y_out_BB ~ ., data = esc_BB)
BB_TE_pred = predict(cons_BB_TE_model,newdata = case5_cons)
rank_BB_TE = which(rank == min(Y[which(BB_TE_pred==min(BB_TE_pred))]))[1]


### PWO model

Xf = design_to_PWO(D_full)
filter <- PWO_filter(cons=cons, Z=Xf)
Xf <- Xf[,-filter$cons_cols]
case5_cons = as.data.frame(Xf[,-1])
PWO_out_BB <- design_to_PWO(D_out_BB)
PWO_out_BB <- PWO_out_BB[,-filter$cons_cols] 
Y_out_BB <- rep(0,nrow(D_out_BB))
# separaPWO regressed and predicPWOd data
for (i in 1:nrow(D_out_BB)) {
  Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  #Y_out_BB[i] = Y[which(apply(D_out_BB, 1, function(x) identical(id,pasPWO(x,sep='',collapse=''))))]
}
PWO_out_BB = PWO_out_BB[,-1]
esc_BB = as.data.frame(cbind(PWO_out_BB,Y_out_BB))
cons_BB_PWO_model = lm(Y_out_BB ~ ., data = esc_BB)
BB_PWO_pred = predict(cons_BB_PWO_model,newdata = case5_cons)
rank_BB_PWO = which(rank == min(Y[which(BB_PWO_pred==min(BB_PWO_pred))]))[1]


########################################
# Step 3: Estimate the optimal order.
# Compare to the true optimal order.
########################################

c(rank_BB_TE,rank_BB_PWO)





##################################################
# Step 4: two by two comparison of TE-1, TE-2, PWO
# and bubble sorting, GRASP.
##################################################

n = 200

result1 = data.frame(matrix(nrow=100,ncol = 6))
resultAIC = data.frame(matrix(nrow=100,ncol = 6))
resultBIC = data.frame(matrix(nrow=100,ncol = 6))
resultR2 = data.frame(matrix(nrow=100,ncol = 6))
resultadjR2 = data.frame(matrix(nrow=100,ncol = 6))

colnames(result1) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')
colnames(resultAIC) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')
colnames(resultBIC) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')
colnames(resultR2) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')
colnames(resultadjR2) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')

N = 100 # the number of runs we want to run for each combination

for (iter in 1:N) {  # make 100 runs for each combination
  
  # BB
  Dtest = D_full[sample(nrow(D_full),n),]
  D_out_BB = bubble_sort_design_B1(init_runs = Dtest,B,criterion="I",block=block)
  
  # GRASP
  D_out_GRASP = grasp_v3_B(D = Dtest,B,block=block,criterion="I", max_iter = 10, n_rand = 250)
  
  # TE1*BB
  
  Xf = OofA_transition_model_matrix(D_full)
  filter <- TE_filter1(cons=cons, D=D_full, Z=Xf)
  Xf <- Xf[,-filter$cons_cols]
  case5_cons = as.data.frame(Xf[,-1])
  
  TE_out_BB <- OofA_transition_model_matrix(D_out_BB)
  TE_out_BB <- TE_out_BB[,-filter$cons_cols] 
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  }
  TE_out_BB = TE_out_BB[,-1]
  esc_BB = as.data.frame(cbind(TE_out_BB,Y_out_BB))
  cons_BB_TE_model = lm(Y_out_BB ~ ., data = esc_BB)
  BB_TE_pred = predict(cons_BB_TE_model,newdata = case5_cons)
  rank_BB_TE = which(rank == min(Y[which(BB_TE_pred==min(BB_TE_pred))]))[1]
  result1[iter,1] = rank_BB_TE
  resultAIC[iter,1] = AIC(cons_BB_TE_model)
  resultBIC[iter,1] = BIC(cons_BB_TE_model)
  resultR2[iter,1] = summary(cons_BB_TE_model)$r.squared
  resultadjR2[iter,1] = summary(cons_BB_TE_model)$adj.r.squared
  
  # TE1*GRASP
  TE_out_GRASP <- OofA_transition_model_matrix(D_out_GRASP)
  TE_out_GRASP <- TE_out_GRASP[,-filter$cons_cols] 
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  TE_out_GRASP = TE_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(TE_out_GRASP,Y_out_GRASP))
  cons_GRASP_TE_model = lm(Y_out_GRASP ~ ., data = esc_GRASP)
  GRASP_TE_pred = predict(cons_GRASP_TE_model,newdata = case5_cons)
  rank_GRASP_TE = which(rank == min(Y[which(GRASP_TE_pred==min(GRASP_TE_pred))]))[1]
  result1[iter,2] = rank_GRASP_TE
  resultAIC[iter,2] = AIC(cons_GRASP_TE_model)
  resultBIC[iter,2] = BIC(cons_GRASP_TE_model)
  resultR2[iter,2] = summary(cons_GRASP_TE_model)$r.squared
  resultadjR2[iter,2] = summary(cons_GRASP_TE_model)$adj.r.squared
  
  # TE2*BB
  
  Xf = OofA_transition_model_matrix_2(D_full)
  filter <- TE_filter2(cons=cons, D=D_full, Z=Xf)
  Xf <- Xf[,-filter$cons_cols]
  case5_cons = as.data.frame(Xf[,-1])
  
  TE_out_BB <- OofA_transition_model_matrix_2(D_out_BB)
  TE_out_BB <- TE_out_BB[,-filter$cons_cols] 
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  }
  TE_out_BB = TE_out_BB[,-1]
  esc_BB = as.data.frame(cbind(TE_out_BB,Y_out_BB))
  cons_BB_TE_model = lm(Y_out_BB ~ ., data = esc_BB)
  BB_TE_pred = predict(cons_BB_TE_model,newdata = case5_cons)
  rank_BB_TE = which(rank == min(Y[which(BB_TE_pred==min(BB_TE_pred))]))[1]
  result1[iter,3] = rank_BB_TE
  resultAIC[iter,3] = AIC(cons_BB_TE_model)
  resultBIC[iter,3] = BIC(cons_BB_TE_model)
  resultR2[iter,3] = summary(cons_BB_TE_model)$r.squared
  resultadjR2[iter,3] = summary(cons_BB_TE_model)$adj.r.squared
  
  # TE2*GRASP
  TE_out_GRASP <- OofA_transition_model_matrix_2(D_out_GRASP)
  TE_out_GRASP <- TE_out_GRASP[,-filter$cons_cols] 
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  TE_out_GRASP = TE_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(TE_out_GRASP,Y_out_GRASP))
  cons_GRASP_TE_model = lm(Y_out_GRASP ~ ., data = esc_GRASP)
  GRASP_TE_pred = predict(cons_GRASP_TE_model,newdata = case5_cons)
  rank_GRASP_TE = which(rank == min(Y[which(GRASP_TE_pred==min(GRASP_TE_pred))]))[1]
  result1[iter,4] = rank_GRASP_TE
  resultAIC[iter,4] = AIC(cons_GRASP_TE_model)
  resultBIC[iter,4] = BIC(cons_GRASP_TE_model)
  resultR2[iter,4] = summary(cons_GRASP_TE_model)$r.squared
  resultadjR2[iter,4] = summary(cons_GRASP_TE_model)$adj.r.squared
  
  # PWO*BB
  Xf = design_to_PWO(D_full)
  filter <- PWO_filter(cons=cons, Z=Xf)
  Xf <- Xf[,-filter$cons_cols]
  case5_cons = as.data.frame(Xf[,-1])
  
  PWO_out_BB <- design_to_PWO(D_out_BB)
  PWO_out_BB <- PWO_out_BB[,-filter$cons_cols] 
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separaPWO regressed and predicPWOd data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
    #Y_out_BB[i] = Y[which(apply(D_out_BB, 1, function(x) identical(id,pasPWO(x,sep='',collapse=''))))]
  }
  PWO_out_BB = PWO_out_BB[,-1]
  esc_BB = as.data.frame(cbind(PWO_out_BB,Y_out_BB))
  cons_BB_PWO_model = lm(Y_out_BB ~ ., data = esc_BB)
  BB_PWO_pred = predict(cons_BB_PWO_model,newdata = case5_cons)
  rank_BB_PWO = which(rank == min(Y[which(BB_PWO_pred==min(BB_PWO_pred))]))[1]
  result1[iter,5] = rank_BB_PWO
  resultAIC[iter,5] = AIC(cons_BB_PWO_model)
  resultBIC[iter,5] = BIC(cons_BB_PWO_model)
  resultR2[iter,5] = summary(cons_BB_PWO_model)$r.squared
  resultadjR2[iter,5] = summary(cons_BB_PWO_model)$adj.r.squared
  
  
  # PWO*GRASP
  PWO_out_GRASP <- design_to_PWO(D_out_GRASP)
  PWO_out_GRASP <- PWO_out_GRASP[,-filter$cons_cols] 
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separaPWO regressed and predicPWOd data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
    #Y_out_GRASP[i] = Y[which(apply(D_out_GRASP, 1, function(x) identical(id,pasPWO(x,sep='',collapse=''))))]
  }
  PWO_out_GRASP = PWO_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(PWO_out_GRASP,Y_out_GRASP))
  cons_GRASP_PWO_model = lm(Y_out_GRASP ~ ., data = esc_GRASP)
  GRASP_PWO_pred = predict(cons_GRASP_PWO_model,newdata = case5_cons)
  rank_GRASP_PWO = which(rank == min(Y[which(GRASP_PWO_pred==min(GRASP_PWO_pred))]))[1]
  result1[iter,6] = rank_GRASP_PWO
  resultAIC[iter,6] = AIC(cons_GRASP_PWO_model)
  resultBIC[iter,6] = BIC(cons_GRASP_PWO_model)
  resultR2[iter,6] = summary(cons_GRASP_PWO_model)$r.squared
  resultadjR2[iter,6] = summary(cons_GRASP_PWO_model)$adj.r.squared
  
}


n = 150

result2 = data.frame(matrix(nrow=100,ncol = 6))
colnames(result2) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')

N = 100 # the number of runs we want to run for each combination

for (iter in 1:N) {  # make 100 runs for each combination
  
  # BB
  Dtest = D_full[sample(nrow(D_full),n),]
  D_out_BB = bubble_sort_design_B1(init_runs = Dtest,B,criterion="I",block=block)
  
  # GRASP
  D_out_GRASP = grasp_v3_B(D = Dtest,B,block=block,criterion="I", max_iter = 10, n_rand = 250)
  
  # TE1*BB
  
  Xf = OofA_transition_model_matrix(D_full)
  filter <- TE_filter1(cons=cons, D=D_full, Z=Xf)
  Xf <- Xf[,-filter$cons_cols]
  case5_cons = as.data.frame(Xf[,-1])
  
  TE_out_BB <- OofA_transition_model_matrix(D_out_BB)
  TE_out_BB <- TE_out_BB[,-filter$cons_cols] 
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  }
  TE_out_BB = TE_out_BB[,-1]
  esc_BB = as.data.frame(cbind(TE_out_BB,Y_out_BB))
  cons_BB_TE_model = lm(Y_out_BB ~ ., data = esc_BB)
  BB_TE_pred = predict(cons_BB_TE_model,newdata = case5_cons)
  rank_BB_TE = which(rank == min(Y[which(BB_TE_pred==min(BB_TE_pred))]))[1]
  result2[iter,1] = rank_BB_TE
  
  # TE1*GRASP
  TE_out_GRASP <- OofA_transition_model_matrix(D_out_GRASP)
  TE_out_GRASP <- TE_out_GRASP[,-filter$cons_cols] 
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  TE_out_GRASP = TE_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(TE_out_GRASP,Y_out_GRASP))
  cons_GRASP_TE_model = lm(Y_out_GRASP ~ ., data = esc_GRASP)
  GRASP_TE_pred = predict(cons_GRASP_TE_model,newdata = case5_cons)
  rank_GRASP_TE = which(rank == min(Y[which(GRASP_TE_pred==min(GRASP_TE_pred))]))[1]
  result2[iter,2] = rank_GRASP_TE
  
  # TE2*BB
  
  Xf = OofA_transition_model_matrix_2(D_full)
  filter <- TE_filter2(cons=cons, D=D_full, Z=Xf)
  Xf <- Xf[,-filter$cons_cols]
  case5_cons = as.data.frame(Xf[,-1])
  
  TE_out_BB <- OofA_transition_model_matrix_2(D_out_BB)
  TE_out_BB <- TE_out_BB[,-filter$cons_cols] 
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  }
  TE_out_BB = TE_out_BB[,-1]
  esc_BB = as.data.frame(cbind(TE_out_BB,Y_out_BB))
  cons_BB_TE_model = lm(Y_out_BB ~ ., data = esc_BB)
  BB_TE_pred = predict(cons_BB_TE_model,newdata = case5_cons)
  rank_BB_TE = which(rank == min(Y[which(BB_TE_pred==min(BB_TE_pred))]))[1]
  result2[iter,3] = rank_BB_TE
  
  # TE2*GRASP
  TE_out_GRASP <- OofA_transition_model_matrix_2(D_out_GRASP)
  TE_out_GRASP <- TE_out_GRASP[,-filter$cons_cols] 
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  TE_out_GRASP = TE_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(TE_out_GRASP,Y_out_GRASP))
  cons_GRASP_TE_model = lm(Y_out_GRASP ~ ., data = esc_GRASP)
  GRASP_TE_pred = predict(cons_GRASP_TE_model,newdata = case5_cons)
  rank_GRASP_TE = which(rank == min(Y[which(GRASP_TE_pred==min(GRASP_TE_pred))]))[1]
  result2[iter,4] = rank_GRASP_TE
  
  # PWO*BB
  Xf = design_to_PWO(D_full)
  filter <- PWO_filter(cons=cons, Z=Xf)
  Xf <- Xf[,-filter$cons_cols]
  case5_cons = as.data.frame(Xf[,-1])
  
  PWO_out_BB <- design_to_PWO(D_out_BB)
  PWO_out_BB <- PWO_out_BB[,-filter$cons_cols] 
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separaPWO regressed and predicPWOd data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
    #Y_out_BB[i] = Y[which(apply(D_out_BB, 1, function(x) identical(id,pasPWO(x,sep='',collapse=''))))]
  }
  PWO_out_BB = PWO_out_BB[,-1]
  esc_BB = as.data.frame(cbind(PWO_out_BB,Y_out_BB))
  cons_BB_PWO_model = lm(Y_out_BB ~ ., data = esc_BB)
  BB_PWO_pred = predict(cons_BB_PWO_model,newdata = case5_cons)
  rank_BB_PWO = which(rank == min(Y[which(BB_PWO_pred==min(BB_PWO_pred))]))[1]
  result2[iter,5] = rank_BB_PWO
  
  # PWO*GRASP
  PWO_out_GRASP <- design_to_PWO(D_out_GRASP)
  PWO_out_GRASP <- PWO_out_GRASP[,-filter$cons_cols] 
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separaPWO regressed and predicPWOd data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
    #Y_out_GRASP[i] = Y[which(apply(D_out_GRASP, 1, function(x) identical(id,pasPWO(x,sep='',collapse=''))))]
  }
  PWO_out_GRASP = PWO_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(PWO_out_GRASP,Y_out_GRASP))
  cons_GRASP_PWO_model = lm(Y_out_GRASP ~ ., data = esc_GRASP)
  GRASP_PWO_pred = predict(cons_GRASP_PWO_model,newdata = case5_cons)
  rank_GRASP_PWO = which(rank == min(Y[which(GRASP_PWO_pred==min(GRASP_PWO_pred))]))[1]
  result2[iter,6] = rank_GRASP_PWO
  
}

colMeans(result1)
colMeans(result2)


colMeans(resultR2)
colMeans(resultadjR2)
colMeans(resultAIC)
colMeans(resultBIC)






