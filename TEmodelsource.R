library(combinat)

trace <- function(X){
  return(sum(diag(X)))
}

I_efficiency <- function(X, B){
  
  M <- t(X)%*%X/nrow(X)
  if(rcond(M) < 1e-10){
    return(NA)
  }
  
  return(trace( solve(M)%*%B))
  
}

# Takes a model matrix X
# returns D_efficiency
D_efficiency <- function(X){
  k <- ncol(X)
  N <- nrow(X)
  return( (1/N)*det(t(X)%*%X)^(1/k) )
}



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
    
    X2 <- c(1,X2) # add the intercept
    
    names(X2) <- c("Intercept",x2namevec)
    
    return(X2[-length(X2)])
    
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
  
  
  X <- cbind(1,X2)
  colnames(X) <- c("Intercept",x2namevec)
  X <- X[,-ncol(X)]
  return(X)
}


# OofA_transition_model_matrix_2
# Input: D, an n by m matrix of permutations
# Output: X, an n by 2*m*(m-1) matrix of model parameters
# m*(m-1) direct transition effects (length 1)
# m*(m-1) length 2 transition effects
## This design is not identifiable for m = 4, but for m > 4 it is fine.
OofA_transition_model_matrix_2 <- function(D){
  n <- nrow(D)
  m <- ncol(D)
  X2 <- matrix(0, nrow = n, ncol = m*(m-1))
  X3 <- matrix(0, nrow = n, ncol = m*(m-1))
  x2namevec <- c()
  x3namevec <- c()
  counter <- 1
  for(i in 1:m){
    for(j in 1:m){
      if(i != j){
        x2namevec <- c(x2namevec, paste("x", i, "|", j, sep = ""))
        x3namevec <- c(x3namevec, paste("z", i, "|", j, sep = ""))
      }
    }
  }
  
  
  for(i in 1:n){
    
    for(j in 2:m){
      
      prev <- D[i,j-1]
      curr <- D[i,j]
      col_ind <- which(x2namevec == paste("x", prev, "|", curr, sep = ""))
      X2[i,col_ind] <- 1
      
      if(j < m){
        
        next_component <- D[i,j+1]
        col_ind2 <- which(x3namevec == paste("z",prev,"|",next_component,sep = ""))
        X3[i,col_ind2] <- 1
        
      }
      
    }
    
    
  }
  
  
  
  colnames(X2) <- x2namevec
  X2 <- X2[,-ncol(X2)]
  colnames(X3) <- x3namevec
  X3 <- X3[,-ncol(X3)] # remove last column for identifiability, since every perm must have m-2 length 2 transitions.
  X <- cbind(1,X2,X3)
  colnames(X) <- c("Intercept",colnames(X2),colnames(X3))
  return(X)
}


### generates the full design with all m! permutations
generate_full_design <- function(m){
  Plist <- permn(1:m) ## get all permutations of 1:m
  return(matrix(unlist(Plist),byrow=TRUE,nrow=factorial(m)))
}


### simulated annealing
# B = t(Xfull)%*%Xfull/nrow(Xfull), which is required for I-optimality
simulated_annealing <- function(D,B, max_iter = 1e4, temp_init = 1, criterion = "I"){
  Dold <- D
  m <- ncol(Dold)
  n <- nrow(Dold)
  Xold <- OofA_transition_model_matrix(Dold)
  Cold <- solve(t(Xold)%*%Xold)
  if(criterion == "D"){
    det_old <- det(t(Xold)%*%Xold)
    criterion_old <- det_old^(1/ncol(Xold))/n
    criterion_best <- criterion_old
  }
  if(criterion == "I"){
    criterion_old <- trace(Cold%*%B)/n # update: works much better
    criterion_best <- criterion_old 
  }
  
  objective_function <- numeric(max_iter)
  Dbest <- Dold
  for (iter in 1:max_iter) {
    # print(paste("Iteration", iter, sep = " "))
    Dnew <- Dold
    i <- sample(1:n, 1)  # random row
    j<- sample(1:(m-1),1) # random column
    # exchange component j with component (j+1)
    Dnew[i,j] <- Dold[i,j+1]
    Dnew[i,j+1] <- Dold[i, j]
    XY <- OofA_transition_model_matrix(rbind(Dold[i,], Dnew[i,]))
    x <- XY[1,]
    y <- XY[2,]
    F1 <- cbind(y,-x)
    F2 <- cbind(y,x)
    I2 <- diag(2)
    Cnew <- Cold - Cold%*%F1%*%solve(I2 + t(F2)%*%Cold%*%F1)%*%t(F2)%*%Cold
    
    temp <- temp_init/(iter+1)
    
    if(criterion == "D"){
      f1 <- F2[,2]
      f2 <- F2[,1]
      v12 <- t(f1)%*%Cold%*%f2
      v1 <- t(f1)%*%Cold%*%f1
      v2 <- t(f2)%*%Cold%*%f2
      delta <- 1 + v2 - v1 + v12^2 - v2*v1
      det_new <- det_old*delta
      criterion_new <-  det_new^(1/ncol(Xold))/n
      if(!is.na(criterion_new) & !is.na(criterion_old)){
        aprob <- exp(-(criterion_old - criterion_new)/temp)
      }
      else{
        aprob <- -1
      }
      
    }
    if(criterion == "I"){
      # Inew <- trace(Cnew%*%B)
      criterion_new <- trace(Cnew%*%B)/n # update: works much better
      aprob <- exp(-(criterion_new-criterion_old)/temp)
    }
    u = runif(1)
    
    if(criterion == "I" & criterion_new < criterion_old){
      criterion_best <- criterion_old 
      Dbest <- Dold
    }
    if(criterion == "D" & criterion_new > criterion_old){
      criterion_best <- criterion_old
      Dbest <- Dold
    }
    
    if(u<aprob){
      criterion_old <- criterion_new
      Dold <- Dnew
      Cold <- Cnew
      if(criterion == "D"){
        det_old <- det_new
      }
    }
    
    objective_function[iter] <- criterion_best
    
    
  }
  plot(1:max_iter, objective_function, type = "l", ylab = "Objective Function", xlab = "Iteration")
  # return(Dbest)
  return(Dbest)
}

# generate_random_design
# m: number of components
# n: sample size
# output: an n times m design matrix where each row is a permutation of 1:m
generate_random_design <- function(m, n){
  
  D0 <- matrix(NA, nrow = n, ncol = m)
  D0 <- t(apply(D0,1,function(x){sample(1:m)}))
  
  
  return(D0)
  
}

# uses a closed-form construction to get the moment matrix for the
# full design with all m! permutations
get_moment_matrix <- function(m){
  q <- 2*choose(m,2) - 1
  Iq <- diag(q)
  V <- matrix(1, nrow = q, ncol = q)
  colnames(V) <- names(OofA_transition_model_matrix(generate_random_design(m,1)))[-1]
  rownames(V) <- colnames(V)
  strsbefore <- sub("\\|.*", "", rownames(V))
  rowis <- as.numeric(substr(strsbefore,2,nchar(strsbefore)))
  strsafter <- sub(".*\\|","",rownames(V))
  rowjs <- as.numeric(strsafter)
  colks <- rowis
  colells <- rowjs
  for(row in 1:nrow(V)){
    i <- rowis[row]
    j <- rowjs[row]
    for(col in 1:ncol(V)){
      k <- colks[col]
      ell <- colells[col]
      if(i == k || j == ell){
        V[row,col] = 0
      }
      if(i == ell && j == k){
        V[row,col] = 0
      }
      
    }
    
  }
  
  onesvec <- matrix(1, nrow = 1, ncol = q+1)
  M22 <- (factorial(m-1)*Iq + factorial(m-2)*V)/factorial(m)
  M <- rbind(onesvec, cbind(1,M22))
  M[1,1] <- 1
  M[1,-1] <- 1/m
  M[-1,1] <- 1/m
  return(M)
}



# Simulate the distribution of random changes in I-efficiency
get_random_effchanges <- function(D0, B, n_sample = 100, criterion = "I"){
  
  
  current_design <- D0
  n <- nrow(D0)
  m <- ncol(D0)
  eff_changes <- numeric(n_sample)
  
  Z_current <- OofA_transition_model_matrix(current_design)
  C_current <- solve(t(Z_current)%*%Z_current)
  
  
  if(criterion == "I"){
    I_current <- trace(C_current%*%B)/n
  }
  if(criterion == "D"){
    det_current <- det(t(Z_current)%*%Z_current)
    d_current <- det_current^(1/ncol(Z_current))/nrow(Z_current)
  }
  
  
  for(i in 1:n_sample){
    
    valid_swap <- FALSE
    while(!valid_swap){
      
      tmp_design <- current_design
      # only say the swap is valid if the resulting matrix is not singular
      i <- sample(1:n,1)
      j <- sample(1:(m-1), 1)
      tmp_design[i,j] <- current_design[i,j+1]
      tmp_design[i,j+1] <- current_design[i,j]
      # Z_tmp <- OofA_transition_model_matrix(tmp_design)
      z_tmpi <- OofA_transition_model_matrix(tmp_design[i,])
      Z_tmp <- Z_current
      Z_tmp[i,] <- z_tmpi
      F1_tmp <- cbind(Z_tmp[i,], -Z_current[i,])
      F2_tmp <- cbind(Z_tmp[i,], Z_current[i,])
      if(rcond(diag(2) + t(F2_tmp)%*%C_current%*%F1_tmp) > 1e-10){
        F1 <- F1_tmp
        F2 <- F2_tmp
        tmp <- solve(diag(2) + t(F2)%*%C_current%*%F1)
        C_new <- C_current - C_current%*%F1%*%tmp%*%t(F2)%*%C_current
        
        if(criterion == "I"){
          I_new <- trace(C_new%*%B)/n
          eff_changes[i] <- abs(I_current - I_new)
        }
        if(criterion == "D"){
          f1 <- F2[,2]
          f2 <- F2[,1]
          v12 <- t(f1)%*%C_current%*%f2
          v1 <- t(f1)%*%C_current%*%f1
          v2 <- t(f2)%*%C_current%*%f2
          delta <- 1 + v2 - v1 + v12^2 - v2*v1
          det_new <- det_current*delta
          d_new <-  det_new^(1/ncol(Z_current))/nrow(Z_current)
          eff_changes[i] <- abs(d_new - d_current)
        }
        
        valid_swap <- TRUE
      }
      
    }
    
    
  }
  
  
  
  
  return(eff_changes)
  
  
}


# Implements the bubble sort design form Lin & Peng 2019
# init_runs: initial design with desired number of runs
# B: matrix required to find I-efficiency
# criterion: must be one of "D" or "I"
# n_repeats: number of times the procedure is repeated
bubble_sort_design <- function(init_runs, B, criterion = "D", n_repeats = 5){
  m <- ncol(init_runs)
  current_design <- init_runs
  n <- nrow(init_runs)
  Z_current <- OofA_transition_model_matrix(current_design)
  C_current <- solve(t(Z_current)%*%Z_current)
  if(criterion == "D"){
    det_current <- det(t(Z_current)%*%Z_current)
    d_current <- D_efficiency(Z_current)
  }
  if(criterion == "I"){
    d_current <- I_efficiency(Z_current, B)/n
    # d_current <- A_efficiency(Z_current)
  }
  for(t in 1:n_repeats){
    
    
    for(i in 1:n){
      
      swapped <- TRUE
      while(swapped){
        
        swapped <- FALSE
        for(j in 1:(m-1)){
          
          current_component <- current_design[i,j]
          next_component <- current_design[i,j+1]
          # try exchanging
          exchange_possible <- FALSE
          tmp_design <- current_design
          tmp_design[i,j] <- current_design[i,j+1]
          tmp_design[i,j+1] <- current_design[i,j]
          z_tmpi <- OofA_transition_model_matrix(tmp_design[i,])
          Z_tmp <- Z_current
          Z_tmp[i,] <- z_tmpi
          # Z_tmp <- OofA_transition_model_matrix(tmp_design)
          # F1 <- cbind(Z_tmp[i,], -Z_current[i,])
          # F2 <- cbind(Z_tmp[i,], Z_current[i,])
          F1 <- cbind(z_tmpi, -Z_current[i,])
          F2 <- cbind(z_tmpi, Z_current[i,])
          if(!is.singular.matrix(diag(2) + t(F2)%*%C_current%*%F1)){
            tmp <- solve(diag(2) + t(F2)%*%C_current%*%F1)
            C_tmp <- C_current - C_current%*%F1%*%tmp%*%t(F2)%*%C_current
            if(criterion == "D"){
              f1 <- F2[,2]
              f2 <- F2[,1]
              v12 <- t(f1)%*%C_current%*%f2
              v1 <- t(f1)%*%C_current%*%f1
              v2 <- t(f2)%*%C_current%*%f2
              delta <- 1 + v2 - v1 + v12^2 - v2*v1
              det_tmp <- det_current*delta
              d_tmp <-  det_tmp^(1/ncol(Z_current))/nrow(Z_current)
              # if the exchange gives a better design, do it
              if(d_tmp > d_current){
                current_design <- tmp_design
                d_current <- d_tmp
                det_current <- det_tmp
                Z_current <- Z_tmp
                C_current <- C_tmp
                swapped <- TRUE
                # print(d_current)
              }
              
            }
            
            if(criterion == "I"){
              # d_tmp <- sum(diag(C_tmp)*nrow(Z_tmp))
              d_tmp <- trace(C_tmp%*%B)/n
              # if the exchange gives a better design, do it
              if(d_tmp < d_current){
                current_design <- tmp_design
                d_current <- d_tmp
                Z_current <- Z_tmp
                C_current <- C_tmp
                swapped <- TRUE
                # print(d_current)
              }
              
            }
            
            
            
            
          }
          
          
          
          
          
          
          
        }
        
        
        
        
      }
      
      
      
      
      
      
    }
    
    
    
  }
  
  
  return(current_design)
  
}


# grasp: Greedy Randomized Adaptive Search Procedure
# D: initial design with desired number of runs
# B: matrix required to find I-efficiency
# criterion: must be one of "D" or "I"
# tol: convergence criterion. Default is tol = 1e-8.
# max_iter: maximum numer of iterations, default is 100.
# n_rand: number of random candidates generated at each iteration. Default is n_rand = 250.
grasp <- function(D, B, criterion = "D", tol = 1e-8, max_iter = 100, n_rand = 250){
  
  delta_converge <- Inf
  t <- 1
  D0 <- D
  n <- nrow(D0)
  m <- ncol(D0)
  X0 <- OofA_transition_model_matrix(D0)
  M0 <- t(X0)%*%X0
  det0 <- det(M0)
  C0 <- solve(M0)
  q_pick <- 10
  best_design <- D0
  if(criterion == "D"){
    D_eff_prev <- det0^(1/ncol(X0))/nrow(X0)
    D_eff_best <- D_eff_prev
  }
  if(criterion == "I"){
    I_eff_prev <- trace(M0%*%B/n)
    I_eff_best <- I_eff_prev
  }
  
  
  while(delta_converge > tol && t < max_iter){
    
    print(paste("Iteration", t, sep = " "))
    
    rand_designs <- vector("list",n_rand)  
    rand_efficiencies <- numeric(n_rand)
    rand_Xs <- vector("list",n_rand)
    rand_Ctmps <- vector("list",n_rand)
    rand_dets <- numeric(n_rand)
    
    ## Step 2: Generate random candidate designs
    for(i in 1:n_rand){
      
      
      found_swap <- FALSE
      while(!found_swap){
        
        D_tmp <- D0
        rowi <- sample(1:n, 1)
        colj <- sample(1:(m-1),1)
        D_tmp[rowi,colj] <- D0[rowi,colj+1]
        D_tmp[rowi,colj+1] <- D0[rowi,colj]
        xi_tmp <- OofA_transition_model_matrix(D_tmp[rowi,])
        # X_tmp <- OofA_transition_model_matrix(D_tmp)
        # F1 <- cbind(X_tmp[rowi,], -X0[rowi,])
        # F2 <- cbind(X_tmp[rowi,], X0[rowi,])
        F1 <- cbind(xi_tmp, -X0[rowi,])
        F2 <- cbind(xi_tmp, X0[rowi,])
        if(rcond(diag(2) + t(F2)%*%C0%*%F1) > 1e-10){
          tmp <- solve(diag(2) + t(F2)%*%C0%*%F1)
          C_tmp <- C0 - C0%*%F1%*%tmp%*%t(F2)%*%C0
          rand_Ctmps[[i]] <- C_tmp
          rand_designs[[i]] <- D_tmp
          X_tmp <- X0
          X_tmp[rowi,] <- xi_tmp
          rand_Xs[[i]] <- X_tmp
          if(criterion == "D"){
            f1 <- F2[,2]
            f2 <- F2[,1]
            v12 <- t(f1)%*%C0%*%f2
            v1 <- t(f1)%*%C0%*%f1
            v2 <- t(f2)%*%C0%*%f2
            delta <- 1 + v2 - v1 + v12^2 - v2*v1
            det_tmp <- det0*delta
            d_tmp <-  det_tmp^(1/ncol(X0))/nrow(X0)
            rand_efficiencies[i] <- d_tmp
            rand_dets[i] <- det_tmp
          }
          
          if(criterion == "I"){
            d_tmp <- trace(C_tmp%*%B)
            rand_efficiencies[i] <- d_tmp
          }
          
          found_swap <- TRUE
          
          
          
        }
        
        
        
      }
      
      # Choose the best q designs
      if(criterion == "D"){
        sorted_efficiency_inds <- sort(rand_efficiencies, decreasing = TRUE, index.return = TRUE)$ix
      }
      if(criterion == "I"){
        sorted_efficiency_inds <- sort(rand_efficiencies, decreasing = FALSE, index.return = TRUE)$ix
      }
      
      
    }
    
    
    top_q_inds <- sorted_efficiency_inds[1:q_pick]
    
    # Randomly select one of the best q designs
    # and then update the design and its associated inverse moment matrix and determinant
    next_design_ind <- sample(top_q_inds, 1)
    D1 <- rand_designs[[next_design_ind]]
    D0 <- D1
    X0 <- rand_Xs[[next_design_ind]]
    C0 <- rand_Ctmps[[next_design_ind]]
    det0 <- rand_dets[next_design_ind]
    
    # If the proposed design is better, update the best design.
    if(criterion == "D"){
      D_eff_next <- rand_efficiencies[next_design_ind]
      # Update delta_converge
      delta_converge <- abs(D_eff_next - D_eff_prev)
      D_eff_prev <- D_eff_next
      if(D_eff_next > D_eff_best){
        D_eff_best <- D_eff_next
        best_design <- D0
      }
    }
    if(criterion == "I"){
      I_eff_next <- rand_efficiencies[next_design_ind]
      # Update delta
      delta_converge <- abs(I_eff_next - I_eff_prev)
      I_eff_prev <- I_eff_next
      if(I_eff_next < I_eff_best){
        I_eff_best <- I_eff_next
        best_design <- D0
      }
    }
    
    # Update the iteration counter
    t <- t+1
    
    # Update q_pick
    if(t %% 5 == 0){
      q_pick <- max(q_pick - 1,1)
    }
    
    
    
  }
  
  return(best_design)
  
}


# Implements the bubble sort design form Lin & Peng 2019
# init_runs: initial design with desired number of runs
# B: matrix required to find I-efficiency
# criterion: must be one of "D" or "I"
# n_repeats: number of times the procedure is repeated
# v2: Returns a list of [[1]] the design, [[2]] inverse of moment matrix, and (if criterion == "D") the determinant
bubble_sort_designv2 <- function(init_runs, B, criterion = "D", n_repeats = 1){
  m <- ncol(init_runs)
  current_design <- init_runs
  n <- nrow(init_runs)
  Z_current <- OofA_transition_model_matrix(current_design)
  C_current <- solve(t(Z_current)%*%Z_current)
  if(criterion == "D"){
    det_current <- det(t(Z_current)%*%Z_current)
    d_current <- D_efficiency(Z_current)
  }
  if(criterion == "I"){
    d_current <- I_efficiency(Z_current, B)/n
  }
  for(t in 1:n_repeats){
    
    
    for(i in 1:n){
      
      swapped <- TRUE
      while(swapped){
        
        swapped <- FALSE
        for(j in 1:(m-1)){
          
          # try exchanging
          exchange_possible <- FALSE
          tmp_design <- current_design
          tmp_design[i,j] <- current_design[i,j+1]
          tmp_design[i,j+1] <- current_design[i,j]
          z_tmpi <- OofA_transition_model_matrix(tmp_design[i,])
          Z_tmp <- Z_current
          Z_tmp[i,] <- z_tmpi
          F1 <- cbind(z_tmpi, -Z_current[i,])
          F2 <- cbind(z_tmpi, Z_current[i,])
          if(rcond(diag(2) + t(F2)%*%C_current%*%F1) > 1e-10){
            tmp <- solve(diag(2) + t(F2)%*%C_current%*%F1)
            C_tmp <- C_current - C_current%*%F1%*%tmp%*%t(F2)%*%C_current
            if(criterion == "D"){
              f1 <- F2[,2]
              f2 <- F2[,1]
              v12 <- t(f1)%*%C_current%*%f2
              v1 <- t(f1)%*%C_current%*%f1
              v2 <- t(f2)%*%C_current%*%f2
              delta <- 1 + v2 - v1 + v12^2 - v2*v1
              det_tmp <- det_current*delta
              d_tmp <-  det_tmp^(1/ncol(Z_current))/nrow(Z_current)
              # if the exchange gives a better design, do it
              if(d_tmp > d_current){
                current_design <- tmp_design
                d_current <- d_tmp
                det_current <- det_tmp
                Z_current <- Z_tmp
                C_current <- C_tmp
                swapped <- TRUE
              }
              
            }
            
            if(criterion == "I"){
              # d_tmp <- sum(diag(C_tmp)*nrow(Z_tmp))
              d_tmp <- trace(C_tmp%*%B)/n
              # if the exchange gives a better design, do it
              if(d_tmp < d_current){
                current_design <- tmp_design
                d_current <- d_tmp
                Z_current <- Z_tmp
                C_current <- C_tmp
                swapped <- TRUE
                # print(d_current)
              }
              
            }
            
            
            
            
          }
          
          
          
          
          
          
          
        }
        
        
        
        
      }
      
      
      
      
      
      
    }
    
    
    
  }
  
  if(criterion == "D"){
    return(list(current_design, C_current, det_current))
  }
  
  return(list(current_design, C_current))
  
}


# grasp: Greedy Randomized Adaptive Search Procedure
# D: initial design with desired number of runs
# B: matrix required to find I-efficiency
# criterion: must be one of "D" or "I"
# tol: convergence criterion. Default is tol = 1e-8.
# max_iter: maximum numer of iterations, default is 100.
# n_rand: number of random candidates generated at each iteration. Default is n_rand = 250.
grasp_v2 <- function(D, B, criterion = "D", max_iter = 10, n_rand = 250){
  
  D0 <- D
  n <- nrow(D0)
  m <- ncol(D0)
  X0 <- OofA_transition_model_matrix(D0)
  M0 <- t(X0)%*%X0
  det0 <- det(M0)
  C0 <- solve(M0)
  q_pick <- 10
  best_design <- D0
  if(criterion == "D"){
    D_eff_best <- det0^(1/ncol(X0))/nrow(X0)
  }
  if(criterion == "I"){
    I_eff_best <- trace(M0%*%B/n)
  }
  
  for(t in 1:max_iter){
    # while(delta_converge > tol && t < max_iter){
    
    print(paste("Iteration", t, sep = " "))
    if(criterion == "D"){
      print(paste("Best D-efficiency: ", D_eff_best, sep = ""))
    }
    if(criterion == "I"){
      print(paste("Best I-efficiency: ", I_eff_best, sep = ""))
    }
    
    # rand_designs <- vector("list",n_rand)  
    rand_exchanges <- matrix(NA, nrow = n_rand, ncol = 2)  # col 1: row index of exchange, col2: col index of exchange
    rand_efficiencies <- numeric(n_rand)
    
    ####################################################
    # Part 1: Construct Greedy Randomized Solution
    ####################################################
    
    for(i in 1:n_rand){
      
      
      found_swap <- FALSE
      while(!found_swap){
        
        D_tmp <- D0
        rowi <- sample(1:n, 1)
        colj <- sample(1:(m-1),1)
        D_tmp[rowi,colj] <- D0[rowi,colj+1]
        D_tmp[rowi,colj+1] <- D0[rowi,colj]
        xi_tmp <- OofA_transition_model_matrix(D_tmp[rowi,])
        F1 <- cbind(xi_tmp, -X0[rowi,])
        F2 <- cbind(xi_tmp, X0[rowi,])
        if(rcond(diag(2) + t(F2)%*%C0%*%F1) > 1e-10){
          tmp <- solve(diag(2) + t(F2)%*%C0%*%F1)
          C_tmp <- C0 - C0%*%F1%*%tmp%*%t(F2)%*%C0
          # rand_designs[[i]] <- D_tmp
          rand_exchanges[i,] <- c(rowi, colj)
          X_tmp <- X0
          X_tmp[rowi,] <- xi_tmp
          if(criterion == "D"){
            f1 <- F2[,2]
            f2 <- F2[,1]
            v12 <- t(f1)%*%C0%*%f2
            v1 <- t(f1)%*%C0%*%f1
            v2 <- t(f2)%*%C0%*%f2
            delta <- 1 + v2 - v1 + v12^2 - v2*v1
            det_tmp <- det0*delta
            d_tmp <-  det_tmp^(1/ncol(X0))/nrow(X0)
            rand_efficiencies[i] <- d_tmp
          }
          
          if(criterion == "I"){
            d_tmp <- trace(C_tmp%*%B/n)
            rand_efficiencies[i] <- d_tmp
          }
          
          found_swap <- TRUE
          
          
          
        }
        
        
        
      }
      
      
    }
    
    # Choose the best q designs
    if(criterion == "D"){
      sorted_efficiency_inds <- sort(rand_efficiencies, decreasing = TRUE, index.return = TRUE)$ix
    }
    if(criterion == "I"){
      sorted_efficiency_inds <- sort(rand_efficiencies, decreasing = FALSE, index.return = TRUE)$ix
    }
    
    
    top_q_inds <- sorted_efficiency_inds[1:q_pick]
    
    # Randomly select one of the best q designs
    next_design_ind <- sample(top_q_inds, 1)
    selected_exchange <- rand_exchanges[next_design_ind,]
    # now perform the exchange on D0
    rowi <- selected_exchange[1]
    colj <- selected_exchange[2]
    D_tmp <- D0
    D_tmp[rowi,colj] <- D0[rowi,colj+1]
    D_tmp[rowi,colj+1] <- D0[rowi,colj]
    D0 <- D_tmp
    
    
    ####################################################
    # Part 2: Use local search to select new design
    ####################################################
    
    
    out_bubble <- bubble_sort_designv2(init_runs = D0, B = B, criterion = criterion)
    D0 <- out_bubble[[1]]
    C0 <- out_bubble[[2]]
    if(criterion == "D"){
      det0 <- out_bubble[[3]]
      D_eff_next = det0^(1/ncol(X0))/n
    } 
    if(criterion == "I"){
      I_eff_next = trace(C1%*%B/n)
    }
    
    
    
    
    ####################################################
    # Part 3: Update the best solution
    ####################################################
    
    
    # If the proposed design is better, update the best design.
    if(criterion == "D"){
      if(D_eff_next > D_eff_best){
        D_eff_best <- D_eff_next
        best_design <- D0
      }
    }
    if(criterion == "I"){
      if(I_eff_next < I_eff_best){
        I_eff_best <- I_eff_next
        best_design <- D0
      }
    }
    
    # Update q_pick
    q_pick <- max(q_pick - 1,1)
    
    
  }
  
  return(best_design)
  
}


# grasp: Greedy Randomized Adaptive Search Procedure
# D: initial design with desired number of runs
# B: matrix required to find I-efficiency
# criterion: must be one of "D" or "I"
# tol: convergence criterion. Default is tol = 1e-8.
# max_iter: maximum numer of iterations, default is 100.
# n_rand: number of random candidates generated at each iteration. Default is n_rand = 250.
grasp_v3 <- function(D, B, criterion = "D", max_iter = 1000, n_rand = 250){
  
  D0 <- D
  n <- nrow(D0)
  m <- ncol(D0)
  X0 <- OofA_transition_model_matrix(D0)
  M0 <- t(X0)%*%X0
  det0 <- det(M0)
  C0 <- solve(M0)
  q_pick <- 10
  best_design <- D0
  if(criterion == "D"){
    D_eff_best <- det0^(1/ncol(X0))/nrow(X0)
  }
  if(criterion == "I"){
    I_eff_best <- trace(M0%*%B/n)
  }
  
  for(t in 1:max_iter){
    # while(delta_converge > tol && t < max_iter){
    
    print(paste("Iteration", t, sep = " "))
    if(criterion == "D"){
      print(paste("Best D-efficiency: ", D_eff_best, sep = ""))
    }
    if(criterion == "I"){
      print(paste("Best I-efficiency: ", I_eff_best, sep = ""))
    }
    
    # rand_designs <- vector("list",n_rand)  
    rand_exchanges <- matrix(NA, nrow = n_rand, ncol = 2)  # col 1: row index of exchange, col2: col index of exchange
    rand_efficiencies <- numeric(n_rand)
    
    ####################################################
    # Part 1: Construct Greedy Randomized Solution
    ####################################################
    
    for(i in 1:n_rand){
      
      
      found_swap <- FALSE
      while(!found_swap){
        
        D_tmp <- D0
        rowi <- sample(1:n, 1)
        colj <- sample(1:(m-1),1)
        D_tmp[rowi,colj] <- D0[rowi,colj+1]
        D_tmp[rowi,colj+1] <- D0[rowi,colj]
        xi_tmp <- OofA_transition_model_matrix(D_tmp[rowi,])
        F1 <- cbind(xi_tmp, -X0[rowi,])
        F2 <- cbind(xi_tmp, X0[rowi,])
        if(rcond(diag(2) + t(F2)%*%C0%*%F1) > 1e-15){
          tmp <- solve(diag(2) + t(F2)%*%C0%*%F1)
          C_tmp <- C0 - C0%*%F1%*%tmp%*%t(F2)%*%C0
          # rand_designs[[i]] <- D_tmp
          rand_exchanges[i,] <- c(rowi, colj)
          X_tmp <- X0
          X_tmp[rowi,] <- xi_tmp
          if(criterion == "D"){
            f1 <- F2[,2]
            f2 <- F2[,1]
            v12 <- t(f1)%*%C0%*%f2
            v1 <- t(f1)%*%C0%*%f1
            v2 <- t(f2)%*%C0%*%f2
            delta <- 1 + v2 - v1 + v12^2 - v2*v1
            det_tmp <- det0*delta
            d_tmp <-  det_tmp^(1/ncol(X0))/nrow(X0)
            rand_efficiencies[i] <- d_tmp
          }
          
          if(criterion == "I"){
            d_tmp <- trace(C_tmp%*%B/n)
            rand_efficiencies[i] <- d_tmp
          }
          
          found_swap <- TRUE
          
          
          
        }
        
        
        
      }
      
      
    }
    
    # Choose the best q designs
    if(criterion == "D"){
      sorted_efficiency_inds <- sort(rand_efficiencies, decreasing = TRUE, index.return = TRUE)$ix
    }
    if(criterion == "I"){
      sorted_efficiency_inds <- sort(rand_efficiencies, decreasing = FALSE, index.return = TRUE)$ix
    }
    
    
    top_q_inds <- sorted_efficiency_inds[1:q_pick]
    
    # Randomly select one of the best q designs
    next_design_ind <- sample(top_q_inds, 1)
    selected_exchange <- rand_exchanges[next_design_ind,]
    # now perform the exchange on D0
    rowi <- selected_exchange[1]
    colj <- selected_exchange[2]
    D_tmp <- D0
    D_tmp[rowi,colj] <- D0[rowi,colj+1]
    D_tmp[rowi,colj+1] <- D0[rowi,colj]
    D0 <- D_tmp
    xi_tmp <- OofA_transition_model_matrix(D_tmp[rowi,])
    F1 <- cbind(xi_tmp, -X0[rowi,])
    F2 <- cbind(xi_tmp, X0[rowi,])
    tmp <- solve(diag(2) + t(F2)%*%C0%*%F1)
    C_tmp <- C0 - C0%*%F1%*%tmp%*%t(F2)%*%C0
    X_tmp <- X0
    X_tmp[rowi,] <- xi_tmp
    X0 <- X_tmp
    if(criterion == "D"){
      f1 <- F2[,2]
      f2 <- F2[,1]
      v12 <- t(f1)%*%C0%*%f2
      v1 <- t(f1)%*%C0%*%f1
      v2 <- t(f2)%*%C0%*%f2
      delta <- 1 + v2 - v1 + v12^2 - v2*v1
      det0 <- det0*delta
      D_eff_next = det0^(1/ncol(X0))/n
    }
    C0 <- C_tmp
    if(criterion == "I"){
      I_eff_next = trace(C0%*%B/n)
    }
    
    ####################################################
    # Part 2: Use local search to select new design
    ####################################################
    
    
    for(i in 1:n){
      current_row <- D0[i,]
      colj <- sample(1:(m-1),1)
      D_tmp <- D0
      D_tmp[i,colj] <- D0[i,colj+1]
      D_tmp[i,colj+1] <- D0[i,colj]
      xi_tmp <- OofA_transition_model_matrix(D_tmp[i,])
      F1 <- cbind(xi_tmp, -X0[i,])
      F2 <- cbind(xi_tmp, X0[i,])
      if(rcond(diag(2) + t(F2)%*%C0%*%F1) > 1e-15){
        tmp <- solve(diag(2) + t(F2)%*%C0%*%F1)
        C_tmp <- C0 - C0%*%F1%*%tmp%*%t(F2)%*%C0
        X_tmp <- X0
        X_tmp[i,] <- xi_tmp
        if(criterion == "D"){
          f1 <- F2[,2]
          f2 <- F2[,1]
          v12 <- t(f1)%*%C0%*%f2
          v1 <- t(f1)%*%C0%*%f1
          v2 <- t(f2)%*%C0%*%f2
          delta <- 1 + v2 - v1 + v12^2 - v2*v1
          det_tmp <- det0*delta
          d_tmp <-  det_tmp^(1/ncol(X0))/nrow(X0)
          if(d_tmp > D_eff_next){
            D_eff_next <- d_tmp
            D0 <- D_tmp
            X0 <- X_tmp
            C0 <- C_tmp
            det0 <- det_tmp
          }
        }
        
        if(criterion == "I"){
          d_tmp <- trace(C_tmp%*%B/n)
          if(d_tmp < I_eff_next){
            I_eff_next <- d_tmp
            D0 <- D_tmp
            X0 <- X_tmp
            C0 <- C_tmp
          }
        }
        
        
        
        
      }
      
    }
    
    
    
    
    ####################################################
    # Part 3: Update the best solution
    ####################################################
    
    
    # If the proposed design is better, update the best design.
    if(criterion == "D"){
      if(D_eff_next > D_eff_best){
        D_eff_best <- D_eff_next
        best_design <- D0
      }
    }
    if(criterion == "I"){
      if(I_eff_next < I_eff_best){
        I_eff_best <- I_eff_next
        best_design <- D0
      }
    }
    
    # Update q_pick
    if(t %% 10 == 0){
      q_pick <- max(q_pick - 1,1)
    }
    
    
    
  }
  
  return(best_design)
  
}

