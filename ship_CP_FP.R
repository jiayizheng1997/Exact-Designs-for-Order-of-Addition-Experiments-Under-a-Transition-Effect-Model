# This file contains codes for experiments on the shipping data
# CP and FP algorithm results can be generated using this file
# To run this file, please run TE_PWO_CP_FP_modelsource.R first
# This file together with ship_PWO_TE.R are used to generate AIC table and optimal rank table

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

result1 = data.frame(matrix(nrow=100,ncol = 4))
colnames(result1) = c('CP*BB', 'CP*GRASP','FP*BB', 'FP*GRASP')

result1MSPE = data.frame(matrix(nrow=100,ncol = 4))
colnames(result1MSPE) = c('CP*BB', 'CP*GRASP','FP*BB', 'FP*GRASP')

resultAIC = data.frame(matrix(nrow=100,ncol = 4))
resultBIC = data.frame(matrix(nrow=100,ncol = 4))
resultR2 = data.frame(matrix(nrow=100,ncol = 4))
resultadjR2 = data.frame(matrix(nrow=100,ncol = 4))
colnames(resultAIC) = c('CP*BB', 'CP*GRASP','FP*BB', 'FP*GRASP')
colnames(resultBIC) = c('CP*BB', 'CP*GRASP','FP*BB', 'FP*GRASP')
colnames(resultR2) = c('CP*BB', 'CP*GRASP','FP*BB', 'FP*GRASP')
colnames(resultadjR2) = c('CP*BB', 'CP*GRASP','FP*BB', 'FP*GRASP')

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
  result1[iter,1] = rank_BB_CP
  resultAIC[iter,1] = AIC(cons_BB_CP_model)
  resultBIC[iter,1] = BIC(cons_BB_CP_model)
  resultR2[iter,1] = summary(cons_BB_CP_model)$r.squared
  resultadjR2[iter,1] = summary(cons_BB_CP_model)$adj.r.squared
  result1MSPE[iter,1] = mean((BB_CP_pred-Y)^2)
  
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
  result1[iter,2] = rank_GRASP_CP
  resultAIC[iter,2] = AIC(cons_GRASP_CP_model)
  resultBIC[iter,2] = BIC(cons_GRASP_CP_model)
  resultR2[iter,2] = summary(cons_GRASP_CP_model)$r.squared
  resultadjR2[iter,2] = summary(cons_GRASP_CP_model)$adj.r.squared
  result1MSPE[iter,2] = mean((GRASP_CP_pred-Y)^2)
  
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
  cons_BB_FP_model = lm(Y_out_BB ~., data = esc_BB)
  BB_FP_pred = predict(cons_BB_FP_model,newdata = case5_cons)
  rank_BB_FP = which(rank == min(Y[which(BB_FP_pred==min(BB_FP_pred))]))[1]
  result1[iter,3] = rank_BB_FP
  resultAIC[iter,3] = AIC(cons_BB_FP_model)
  resultBIC[iter,3] = BIC(cons_BB_FP_model)
  resultR2[iter,3] = summary(cons_BB_FP_model)$r.squared
  resultadjR2[iter,3] = summary(cons_BB_FP_model)$adj.r.squared
  result1MSPE[iter,3] = mean((BB_FP_pred-Y)^2)
  
  # FP*GRASP
  FP_out_GRASP <- design_to_FP(D_out_GRASP)
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  FP_out_GRASP = FP_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(FP_out_GRASP,Y_out_GRASP))
  cons_GRASP_FP_model = lm(Y_out_GRASP ~., data = esc_GRASP)
  GRASP_FP_pred = predict(cons_GRASP_FP_model,newdata = case5_cons)
  rank_GRASP_FP = which(rank == min(Y[which(GRASP_FP_pred==min(GRASP_FP_pred))]))[1]
  result1[iter,4] = rank_GRASP_FP
  resultAIC[iter,4] = AIC(cons_GRASP_FP_model)
  resultBIC[iter,4] = BIC(cons_GRASP_FP_model)
  resultR2[iter,4] = summary(cons_GRASP_FP_model)$r.squared
  resultadjR2[iter,4] = summary(cons_GRASP_FP_model)$adj.r.squared
  result1MSPE[iter,4] = mean((GRASP_FP_pred-Y)^2)
}


n = 100

result2 = data.frame(matrix(nrow=100,ncol = 4))
colnames(result2) = c('CP*BB', 'CP*GRASP','FP*BB', 'FP*GRASP')

result2MSPE = data.frame(matrix(nrow=100,ncol = 4))
colnames(result2MSPE) = c('CP*BB', 'CP*GRASP','FP*BB', 'FP*GRASP')

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
  result2[iter,1] = rank_BB_CP
  result2MSPE[iter,1] = mean((BB_CP_pred-Y)^2)
  
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
  result2[iter,2] = rank_GRASP_CP
  result2MSPE[iter,2] = mean((GRASP_CP_pred-Y)^2)
  
  # FP*BB
  case5_cons = as.data.frame(design_to_FP(D_full))
  Xf = OofA_transition_model_matrix_2(D_full)
  filter <- TE_filter1(cons=cons, block=block, D=D_full, Z=Xf)
  Xf <- Xf[,-filter$cons_cols]
  
  FP_out_BB <- design_to_FP(D_out_BB)
  Y_out_BB <- rep(0,nrow(D_out_BB))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_BB)) {
    Y_out_BB[i] = Y[which(id==paste(D_out_BB[i,],sep='',collapse=''))]
  }
  FP_out_BB = FP_out_BB[,-1]
  esc_BB = as.data.frame(cbind(FP_out_BB,Y_out_BB))
  cons_BB_FP_model = lm(Y_out_BB ~., data = esc_BB)
  BB_FP_pred = predict(cons_BB_FP_model,newdata = case5_cons)
  rank_BB_FP = which(rank == min(Y[which(BB_FP_pred==min(BB_FP_pred))]))[1]
  result2[iter,3] = rank_BB_FP
  result2MSPE[iter,3] = mean((BB_FP_pred-Y)^2)
  
  # FP*GRASP
  FP_out_GRASP <- design_to_FP(D_out_GRASP)
  Y_out_GRASP <- rep(0,nrow(D_out_GRASP))
  # separate regressed and predicted data
  for (i in 1:nrow(D_out_GRASP)) {
    Y_out_GRASP[i] = Y[which(id==paste(D_out_GRASP[i,],sep='',collapse=''))]
  }
  FP_out_GRASP = FP_out_GRASP[,-1]
  esc_GRASP = as.data.frame(cbind(FP_out_GRASP,Y_out_GRASP))
  cons_GRASP_FP_model = lm(Y_out_GRASP ~., data = esc_GRASP)
  GRASP_FP_pred = predict(cons_GRASP_FP_model,newdata = case5_cons)
  rank_GRASP_FP = which(rank == min(Y[which(GRASP_FP_pred==min(GRASP_FP_pred))]))[1]
  result2[iter,4] = rank_GRASP_FP
  result2MSPE[iter,4] = mean((GRASP_FP_pred-Y)^2)
  
}

colMeans(result1)
colMeans(result2)

colMeans(resultR2)
colMeans(resultadjR2)
colMeans(resultAIC)
colMeans(resultBIC)

colMeans(result1MSPE)
colMeans(result2MSPE)















































