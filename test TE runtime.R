mlist = c(9,10,11)
nlist = c(400,500,600)
t = 20
t_SA = c()
t_BB = c()
t_GRASP = c()
for (m in mlist) {
  Dfull = generate_full_design(m)
  Xfull = OofA_transition_model_matrix_2(Dfull)
  B = t(Xfull)%*%Xfull/nrow(Xfull)
  t1 = Sys.time()
  for (n in nlist) {
    for (i in 1:t) {
      D = Dfull[sample(nrow(Dfull),n),]
      simulated_annealing(D,B, max_iter = 1e4, temp_init = 1, criterion = "I")
    }
  }
  t2 = Sys.time()
  t_SA = c(t_SA,t2-t1)
  t1 = Sys.time()
  for (n in nlist) {
    for (i in 1:t) {
      init_runs = sample(nrow(Dfull),n)
      bubble_sort_designv2(init_runs, B, criterion = "I", n_repeats = 5)
    }
  }
  t2 = Sys.time()
  t_BB = c(t_BB,t2-t1)
  t1 = Sys.time()
  for (n in nlist) {
    for (i in 1:t) {
      D = Dfull[sample(nrow(Dfull),n),]
      grasp_v3(D, B, criterion = "I", max_iter = 10, n_rand = 250)
    }
  }
  t2 = Sys.time()
  t_GRASP = c(t_GRASP,t2-t1)
}
















