Y <- c(t(read.table("ps1.phenos")))
X <- read.table("ps1.genos")
N <- length(Y)
B <- 1e+4

T_mat<-matrix(vector(), nrow=B, ncol=10)
for (i in 1:B) {
  Y_permutation <- sample(Y)
  for(j in 1:10){
  	T_mat[i,j] <- N*(cor(X[,j],Y_permutation)^2)
  }
}

p_mat<-matrix(vector(), nrow=B, ncol=10)
for (i in 1:B) {
   for(j in 1:10){
   	p_mat[i,j]<-sum(T_mat[,j]>=T_mat[i,j])/B
   }
}

fp<-0
p_threshold<-0.0087
for(i in 1:B){
	for(j in 1:10){
		if(p_mat[i,j] <= p_threshold){
			fp<-fp+1
			break
		}
	}
}
fp

jpeg('rplot.jpg')
h <- hist(T1,ylab="Value of T1")
abline(v=T1_observed, col="blue")
dev.off()
T1_observed <- N*(cor(X[,1],Y)^2)
sum(T1>=T1_observed)/B #0.00761

library(ggplot2)
ggplot(data.frame(x = c(0, 20)), aes(x = x)) +
     stat_function(fun = dchisq, args = list(df = 1))
ggsave("chi_squared_1.png")

pchisq(T1_observed, df=1, lower.tail=FALSE) # [1] 0.007608864

##0.005