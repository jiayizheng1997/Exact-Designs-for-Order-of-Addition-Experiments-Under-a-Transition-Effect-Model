# This file contains codes for experiments on the shipping data
# CP and FP algorithm results can be generated using this file
# To run this file, please run TE_PWO_CP_FP_modelsource.R first
# This file together with ship_PWO_TE.R are used to generate AIC table and optimal rank table

# To begin with, let's have a look at how the table ExpectedTravelTime table was obtained
# There are six different terminals and one start point
# The first row means, the time it takes from starting point to finish every cargo/nodes
# The second to the 21th row mean that, the time it takes from starting point of cargo i to the end point of cargo j
# For example, 10.275, the 2nd number from 2nd row to 6th row, means that it would take 10.275 from the start point of cargos 1-5 to the end point of cargo 1
# It is just the 1st number in the LoadAndCleanTime table, as no travelling is required.
# Similarly, 18.254, the 2nd number from 7th to 9th row, equal to 7.979+10.27, the time to travel from terminal 2 to terminal 1 plus the load-and-clean time for cargo 1
# From this, we can deduce that the SD table can be reduced to the travelling time standard deviation from terminal to terminal.

# Reminder: manually change n every time if needed, table B2 and table B3 e.g.

set.seed(1234)

SD <- read.table("SD.txt")
ExpectedTerminalTravelTime <- read.table("ExpectedTerminalTravelTime.txt")

m=7
D_full <- generate_full_design(m)
block = matrix(c(1,rep(0,5),seq(2,7)),nrow=2,byrow = TRUE)
cons <-cons_pair(block)

Xf = OofA_transition_model_matrix(D_full)
filter <- TE_filter1(cons=cons, D=D_full, Z=Xf)
D_full <- D_full[-filter$cons_rows,]
Xf <- Xf[,-filter$cons_cols]
B = t(Xf)%*%Xf/nrow(Xf)
id <- apply(D_full,1, function(x) paste(x,sep='',collapse=''))



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
result1MSPE = data.frame(matrix(nrow=100,ncol = 6))

colnames(result1) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')
colnames(resultAIC) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')
colnames(resultBIC) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')
colnames(resultR2) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')
colnames(resultadjR2) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')
colnames(result1MSPE) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')

N = 100 # the number of runs we want to run for each combination

for (iter in 1:N) {  # make 100 runs for each combination
  
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
  cons_BB_TE_model = lm(Y_out_BB ~., data = esc_BB)
  BB_TE_pred = predict(cons_BB_TE_model,newdata = case5_cons)
  rank_BB_TE = which(rank == min(Y[which(BB_TE_pred==min(BB_TE_pred))]))[1]
  result1[iter,1] = rank_BB_TE
  resultAIC[iter,1] = AIC(cons_BB_TE_model)
  resultBIC[iter,1] = BIC(cons_BB_TE_model)
  resultR2[iter,1] = summary(cons_BB_TE_model)$r.squared
  resultadjR2[iter,1] = summary(cons_BB_TE_model)$adj.r.squared
  result1MSPE[iter,1] = mean((BB_TE_pred-Y)^2)
  
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
  cons_GRASP_TE_model = lm(Y_out_GRASP ~., data = esc_GRASP)
  GRASP_TE_pred = predict(cons_GRASP_TE_model,newdata = case5_cons)
  rank_GRASP_TE = which(rank == min(Y[which(GRASP_TE_pred==min(GRASP_TE_pred))]))[1]
  result1[iter,2] = rank_GRASP_TE
  resultAIC[iter,2] = AIC(cons_GRASP_TE_model)
  resultBIC[iter,2] = BIC(cons_GRASP_TE_model)
  resultR2[iter,2] = summary(cons_GRASP_TE_model)$r.squared
  resultadjR2[iter,2] = summary(cons_GRASP_TE_model)$adj.r.squared
  result1MSPE[iter,2] = mean((GRASP_TE_pred-Y)^2)
  
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
  cons_BB_TE_model = lm(Y_out_BB ~., data = esc_BB)
  BB_TE_pred = predict(cons_BB_TE_model,newdata = case5_cons)
  rank_BB_TE = which(rank == min(Y[which(BB_TE_pred==min(BB_TE_pred))]))[1]
  result1[iter,3] = rank_BB_TE
  resultAIC[iter,3] = AIC(cons_BB_TE_model)
  resultBIC[iter,3] = BIC(cons_BB_TE_model)
  resultR2[iter,3] = summary(cons_BB_TE_model)$r.squared
  resultadjR2[iter,3] = summary(cons_BB_TE_model)$adj.r.squared
  result1MSPE[iter,3] = mean((BB_TE_pred-Y)^2)
  
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
  cons_GRASP_TE_model = lm(Y_out_GRASP ~., data = esc_GRASP)
  GRASP_TE_pred = predict(cons_GRASP_TE_model,newdata = case5_cons)
  rank_GRASP_TE = which(rank == min(Y[which(GRASP_TE_pred==min(GRASP_TE_pred))]))[1]
  result1[iter,4] = rank_GRASP_TE
  resultAIC[iter,4] = AIC(cons_GRASP_TE_model)
  resultBIC[iter,4] = BIC(cons_GRASP_TE_model)
  resultR2[iter,4] = summary(cons_GRASP_TE_model)$r.squared
  resultadjR2[iter,4] = summary(cons_GRASP_TE_model)$adj.r.squared
  result1MSPE[iter,4] = mean((GRASP_TE_pred-Y)^2)
  
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
  cons_BB_PWO_model = lm(Y_out_BB ~., data = esc_BB)
  BB_PWO_pred = predict(cons_BB_PWO_model,newdata = case5_cons)
  rank_BB_PWO = which(rank == min(Y[which(BB_PWO_pred==min(BB_PWO_pred))]))[1]
  result1[iter,5] = rank_BB_PWO
  resultAIC[iter,5] = AIC(cons_BB_PWO_model)
  resultBIC[iter,5] = BIC(cons_BB_PWO_model)
  resultR2[iter,5] = summary(cons_BB_PWO_model)$r.squared
  resultadjR2[iter,5] = summary(cons_BB_PWO_model)$adj.r.squared
  result1MSPE[iter,5] = mean((BB_PWO_pred-Y)^2)
  
  
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
  cons_GRASP_PWO_model = lm(Y_out_GRASP ~., data = esc_GRASP)
  GRASP_PWO_pred = predict(cons_GRASP_PWO_model,newdata = case5_cons)
  rank_GRASP_PWO = which(rank == min(Y[which(GRASP_PWO_pred==min(GRASP_PWO_pred))]))[1]
  result1[iter,6] = rank_GRASP_PWO
  resultAIC[iter,6] = AIC(cons_GRASP_PWO_model)
  resultBIC[iter,6] = BIC(cons_GRASP_PWO_model)
  resultR2[iter,6] = summary(cons_GRASP_PWO_model)$r.squared
  resultadjR2[iter,6] = summary(cons_GRASP_PWO_model)$adj.r.squared
  result1MSPE[iter,6] = mean((GRASP_PWO_pred-Y)^2)
  
}


n = 100

result2 = data.frame(matrix(nrow=100,ncol = 6))
colnames(result2) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')
result2MSPE = data.frame(matrix(nrow=100,ncol = 6))
colnames(result2MSPE) = c('TE1*BB', 'TE1*GRASP','TE2*BB', 'TE2*GRASP', 'PWO*BB', 'PWO*GRASP')

N = 100 # the number of runs we want to run for each combination

for (iter in 1:N) {  # make 100 runs for each combination
  
  # Generate Y for each iteration as Y is random
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
  cons_BB_TE_model = lm(Y_out_BB ~., data = esc_BB)
  BB_TE_pred = predict(cons_BB_TE_model,newdata = case5_cons)
  rank_BB_TE = which(rank == min(Y[which(BB_TE_pred==min(BB_TE_pred))]))[1]
  result2[iter,1] = rank_BB_TE
  result2MSPE[iter,1] = mean((BB_TE_pred-Y)^2)
  
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
  cons_GRASP_TE_model = lm(Y_out_GRASP ~., data = esc_GRASP)
  GRASP_TE_pred = predict(cons_GRASP_TE_model,newdata = case5_cons)
  rank_GRASP_TE = which(rank == min(Y[which(GRASP_TE_pred==min(GRASP_TE_pred))]))[1]
  result2[iter,2] = rank_GRASP_TE
  result2MSPE[iter,2] = mean((GRASP_TE_pred-Y)^2)
  
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
  cons_BB_TE_model = lm(Y_out_BB ~., data = esc_BB)
  BB_TE_pred = predict(cons_BB_TE_model,newdata = case5_cons)
  rank_BB_TE = which(rank == min(Y[which(BB_TE_pred==min(BB_TE_pred))]))[1]
  result2[iter,3] = rank_BB_TE
  result2MSPE[iter,3] = mean((BB_TE_pred-Y)^2)
  
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
  cons_GRASP_TE_model = lm(Y_out_GRASP ~., data = esc_GRASP)
  GRASP_TE_pred = predict(cons_GRASP_TE_model,newdata = case5_cons)
  rank_GRASP_TE = which(rank == min(Y[which(GRASP_TE_pred==min(GRASP_TE_pred))]))[1]
  result2[iter,4] = rank_GRASP_TE
  result2MSPE[iter,4] = mean((GRASP_TE_pred-Y)^2)
  
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
  cons_BB_PWO_model = lm(Y_out_BB ~., data = esc_BB)
  BB_PWO_pred = predict(cons_BB_PWO_model,newdata = case5_cons)
  rank_BB_PWO = which(rank == min(Y[which(BB_PWO_pred==min(BB_PWO_pred))]))[1]
  result2[iter,5] = rank_BB_PWO
  result2MSPE[iter,5] = mean((BB_PWO_pred-Y)^2)
  
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
  cons_GRASP_PWO_model = lm(Y_out_GRASP ~., data = esc_GRASP)
  GRASP_PWO_pred = predict(cons_GRASP_PWO_model,newdata = case5_cons)
  rank_GRASP_PWO = which(rank == min(Y[which(GRASP_PWO_pred==min(GRASP_PWO_pred))]))[1]
  result2[iter,6] = rank_GRASP_PWO
  result2MSPE[iter,6] = mean((GRASP_PWO_pred-Y)^2)
  
}

colMeans(result1)
colMeans(result2)


colMeans(resultR2)
colMeans(resultadjR2)
colMeans(resultAIC)
colMeans(resultBIC)

colMeans(result1MSPE)
colMeans(result2MSPE)


# 1. AIC, BIC comparison between all methods
# 2. Explicitly write the formulas for the CP and FP models in Section 5 (or, alternatively, in Section 2)
# 3. Try simulating 1 dataset and compare the MSPE, AIC, BIC, and show the parameter estimates for TE1
# 4. Optimal order to appendix










































