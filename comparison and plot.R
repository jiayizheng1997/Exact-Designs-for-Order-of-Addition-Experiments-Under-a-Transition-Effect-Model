#--------------------------------------------
### Comparison between PWO and Simulated Annealing
#--------------------------------------------
t1 = Sys.time()
set.seed(1234)
tot_iter = 100

m = 7
nlist = c(60,90,120)
block = matrix(c(2,4,3,0),byrow = TRUE,nrow = 2)

I_list_PWO = matrix(NA,nrow = tot_iter, ncol = length(nlist))
I_list_annealing = matrix(NA,nrow = tot_iter, ncol = length(nlist))

I_name = rep(0,length(nlist))
for (j in 1:length(nlist)) {
  I_name[j] <- c(paste("m=", m, ",n=", nlist[j], sep = ""))
}


colnames(I_list_PWO) = I_name
colnames(I_list_annealing) = I_name

for (n in 1:length(nlist)) {
  Df <- generate_full_design(m)
  Xf <- design_to_PWO(Df)
  cons <- cons_pair(block)
  filter <- PWO_filter(cons, Xf)
  Df <- Df[-filter$cons_rows,]
  Xf <- Xf[-filter$cons_rows,-filter$cons_cols]
  Bf <- t(Xf)%*%Xf  # this is the B matrix for I optimality
  I_full <- I_efficiency(Xf, Bf) # I optimality criterion for the full design
  for (iter in 1:tot_iter) {
    start_found <- FALSE
    # we need to find a starting design with n rows that is not singular
    while(!start_found){
      D_start <- Df[sample(1:nrow(Df),nlist[n]),]
      X_start <- design_to_PWO(D_start)
      X_start <- X_start[,-filter$cons_cols]
      M_start <- t(X_start)%*%X_start
      if(!is.singular.matrix(M_start)){
        start_found <- TRUE
        print("Found initial matrix")
      }
    }
    D_out_PWO <- random_exchange(D_start, Bf, block = block)
    X_out_PWO <- design_to_PWO(D_out_PWO)
    X_out_PWO <- X_out_PWO[,-filter$cons_cols]
    I_out_PWO <- I_efficiency(X_out_PWO, Bf)
    I_list_PWO[iter,n]=I_full/I_out_PWO
    D_out_annealing <- simulated_annealing(D_start, Bf, block = block)
    X_out_annealing <- design_to_PWO(D_out_annealing)
    X_out_annealing <- X_out_annealing[,-filter$cons_cols]
    I_out_annealing <- I_efficiency(X_out_annealing, Bf)
    I_list_annealing[iter,n]=I_full/I_out_annealing
  }
}



#calculate mean for each n
I_mean_PWO=matrix(apply(I_list_PWO, 2, mean),ncol = length(nlist),byrow = FALSE)
I_mean_annealing=matrix(apply(I_list_annealing, 2, mean),ncol = length(nlist),byrow = FALSE)


rownames(I_mean_PWO)=c(paste("m=", m, sep = ""))
colnames(I_mean_PWO)=c(paste("n=", nlist, sep = ""))

rownames(I_mean_annealing)=c(paste("m=", m, sep = ""))
colnames(I_mean_annealing)=c(paste("n=", nlist, sep = ""))

#Box plots
preplot = cbind(matrix(I_list_PWO,ncol = 1),matrix(I_list_annealing,ncol = 1))

preplot = as.data.frame(cbind(preplot,rep(nlist,each=tot_iter)))

colnames(preplot)=c("RandomExchange","AnnealingSimulation","n")
preplot$n=as.factor(preplot$n)

p1 <- ggplot(preplot, aes(x=n, y=RandomExchange)) + 
  geom_boxplot()+ggtitle("Box plot for different subset sizes") +
  ylab("I-efficiency")
p1
p2 <- ggplot(preplot, aes(x=n, y=AnnealingSimulation)) + 
  geom_boxplot()+ggtitle("Box plot for different subset sizes") +
  ylab("I-efficiency")
p2
preplot_long = melt(preplot,id.vars = c("n"))
p <- ggplot(preplot_long, aes(x=n,y=value,fill=n))+
  geom_boxplot() + labs(title=paste("Algorithm Efficiency Comparison When m = ",m)) + facet_wrap(~variable) +  ylab("Relative I-efficiency")
p
t2 = Sys.time()
t2-t1