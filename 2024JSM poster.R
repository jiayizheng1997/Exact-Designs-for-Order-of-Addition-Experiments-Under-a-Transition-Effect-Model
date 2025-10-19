library(ggplot2)
library(readr)

case2studyextremesigma0_01 <- read_csv("case2studyextremesigma0.01.csv")
data = case2studyextremesigma0_01[,7:11]
alg = c('SA','SA_Voting','IS','Fedrov','Rand')
list1 = c()
list2 = rep(alg,each=nrow(data))
for (i in 1:5) {
  list1 = rbind(list1,as.matrix(data[,i]))
}
data = cbind(list1,as.data.frame(list2))
colnames(data) <- c("Algorithms","AMSPE")
ggplot(data, aes(x=Algorithms, y=AMSPE, color=AMSPE))+geom_boxplot()+xlim(0,5)+coord_flip()+labs(y="Algorithms",x="AMSPE",title = "Extreme Constraints")

case2studynonextremesigma0_01 <- read_csv("case2studynonextremesigma0.01.csv")
data = case2studynonextremesigma0_01[,7:11]
alg = c('SA','SA_Voting','IS','Fedrov','Rand')
list1 = c()
list2 = rep(alg,each=nrow(data))
for (i in 1:5) {
  list1 = rbind(list1,as.matrix(data[,i]))
}
data = cbind(list1,as.data.frame(list2))
colnames(data) <- c("Algorithms","AMSPE")
ggplot(data, aes(x=Algorithms, y=AMSPE, color=AMSPE))+geom_boxplot()+xlim(0,5)+coord_flip()+labs(y="Algorithms",x="AMSPE",title = "Non-extreme Constraints")




df2_out10000max_iter <- read_csv("df2_out10000max_iter.csv")
data = cbind(df2_out10000max_iter$SimulatedAnnealing.8,df2_out10000max_iter$SA_voting.8,df2_out10000max_iter$InsertionSort.8,df2_out10000max_iter$optFederov.8,df2_out10000max_iter$n.32)
alg = c('SA','SA_Voting','IS','Fedrov')
list1 = c()
list2 = rep(alg,each=nrow(data))
list3 = rep(data[,5],4)
for (i in 1:4) {
  list1 = rbind(list1,as.matrix(data[,i]))
}
data = cbind(list1,as.data.frame(list2),as.data.frame(list3))
colnames(data) <- c("Relative_Iefficiency","Algorithms","n")
data$n <- as.factor(data$n)
data <- data[which(data$Relative_Iefficiency>0.1),]


ggplot(data, aes(x=n, y=Relative_Iefficiency, fill=Algorithms))+geom_boxplot()+ylim(0,1)



