# This file is used to create tables 4.1, 4.2, 4.3, 4.4 in dissertation
# manual change parameters is required.


# TEmodel_unconstrained_sim.R

rm(list = ls()) # clear environment


library(matrixcalc)
library(foreach)
library(doParallel)
library(doRNG)

source("TEmodelsource.R")

#######################
# Simulation Settings
#######################

m = 11 # m = 9, 10, 11, 12, .... components
n = 500 # 400, 500, 600
n_search_attempts = 20
criterion = "I" # can be "D" or "I"
max_iter = 10000 # maximum number of iterations for SA

Mfull <- get_moment_matrix(m)
B <- Mfull
Deff_full <- det(Mfull)^(1/ncol(Mfull))
Ieff_full <- 2*choose(m,2)


# test code
initial_design_found <- FALSE
while(!initial_design_found){
  D0 <- generate_random_design(m,n)
  X0 <- OofA_transition_model_matrix(D0)
  M0 <- t(X0)%*%X0/n
  if(rcond(M0) > 1e-10){
    print("Initial design found")
    initial_design_found <- TRUE
  }
}

effchanges <- get_random_effchanges(D0, B, criterion = criterion)
tseq <- seq(from = 1e-5, to = 10000, by = 0.01)
aprobs <- exp(-(max(effchanges, na.rm = T)/(tseq/2)))
temp_init = max(tseq[which(abs(aprobs-0.98) < 0.02)])

if(is.infinite(temp_init)){
    temp_init = 0.05
}

n.cores <- 20
my.cluster <- parallel::makeCluster(n.cores)
registerDoParallel(my.cluster)


simout <- foreach(i = 1:n_search_attempts, .packages = c("matrixcalc","combinat"), .combine = 'rbind') %dorng% {
  
  source("TEmodelsource.R", local = TRUE)
  print(paste("Beginning Search: ", i, sep = ""))
  
  # randomly select initial design
  initial_design_found <- FALSE
  while(!initial_design_found){
    D0 <- generate_random_design(m,n)
    X0 <- OofA_transition_model_matrix(D0)
    M0 <- t(X0)%*%X0/n
    if(rcond(M0) > 1e-10){
      print("Initial design found")
      initial_design_found <- TRUE
    }
  }
  
  t1 = Sys.time()
  D_SA <- simulated_annealing(D0,B = B, temp_init = temp_init, max_iter = max_iter, criterion = criterion)
  X_SA <- OofA_transition_model_matrix(D_SA)
  t2 = Sys.time()
  t_SA = as.numeric(t2-t1)
  
  t1 = Sys.time()
  D_bubble <- bubble_sort_design(D0, B = B, criterion = criterion, n_repeats = 1)
  X_bubble <- OofA_transition_model_matrix(D_bubble)
  t2 = Sys.time()
  t_bubble = as.numeric(t2-t1)
  
  t1 = Sys.time()
  D_grasp <- grasp_v3(D=D0, B=B, criterion  = criterion, n_rand = 30, max_iter = 100)
  X_grasp <- OofA_transition_model_matrix(D_grasp)
  t2 = Sys.time()
  t_grasp = as.numeric(t2-t1)
  
  
  


  rowout <- c(I_efficiency(X_SA, B), I_efficiency(X_bubble, B), I_efficiency(X_grasp,B), D_efficiency(X_SA), D_efficiency(X_bubble), D_efficiency(X_grasp),t_SA, t_bubble, t_grasp)
    

  
  
  
  rowout
}


if(criterion == "I"){
  Ieff_SA <- Ieff_full/simout[,1]
  Ieff_bubble <- Ieff_full/simout[,2]
  Ieff_grasp <- Ieff_full/simout[,3]
  
  simout2 <- cbind(Ieff_SA,Ieff_bubble,Ieff_grasp,simout[,4],simout[,5],simout[,6])

  
}

if(criterion == "D"){
  Deff_SA <- simout[,1]/Deff_full
  Deff_bubble <- simout[,2]/Deff_full
  Deff_grasp <- simout[,3]/Deff_full
  
  simout2 <- cbind(Deff_SA,Deff_bubble,Deff_grasp,simout[,4],simout[,5],simout[,6])

}

colnames(simout2) = c("SA","Bubble","GRASP","t_SA","t_Bubble","t_GRASP")

#namestr_simout <- paste("rel_effs",m,"n",n,criterion,".csv",sep="")
namestr_simout <- paste("robust_",m,"n",n,criterion,".csv",sep="")
write.csv(simout2, namestr_simout)

parallel::stopCluster(my.cluster)
