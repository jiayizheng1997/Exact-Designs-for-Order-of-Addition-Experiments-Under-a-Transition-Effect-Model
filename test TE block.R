library(matlib)

block = matrix(1:4,nrow = 2)
D = generate_full_design(9)
Xf = OofA_transition_model_matrix(D)
cons <-cons_pair(block)
filter <- TE_filter1(cons=cons, block=block, D=D, Z=Xf)
D <- D[-filter$cons_rows,]
#Xf <- Xf[,-filter$cons_cols]
#Xf <- Xf[-filter$cons_rows,]
Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
B = t(Xf)%*%Xf/nrow(Xf)
Dtest = D[sample(nrow(D),100),]
simulated_annealing_B1(Dtest,B,max_iter = 100, block, temp_init = 1, criterion = "I")
bubble_sort_design_B1(init_runs = Dtest,B,block=block)



echelon(B)

# partially constrained model still not work

block = matrix(1:6,nrow = 2)
D = generate_full_design(9)
Xf = OofA_transition_model_matrix(D)
cons <-cons_pair(block)
filter <- TE_filter1(cons=cons, block=block, D=D, Z=Xf)
D <- D[-filter$cons_rows,]
#Xf <- Xf[,-filter$cons_cols]
#Xf <- Xf[-filter$cons_rows,]
Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
B = t(Xf)%*%Xf/nrow(Xf)
Dtest = D[sample(nrow(D),30),]
simulated_annealing_B1(Dtest,B,max_iter = 100, block, temp_init = 1, criterion = "I")



























