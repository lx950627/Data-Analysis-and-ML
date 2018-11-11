sigmoid <- function(z){
	g <- 1/(1+exp(-z))
	return(g)
}

m <- 1000
X_train<-t(read.fwf("hw.2-1.geno", widths = rep(1,m)))
y_train<-read.table("hw.2-1.pheno")
X_train<-as.matrix(X_train)
y_train<-as.matrix(y_train)

X_train<-cbind(rep(1,m),X_train)

getNegativeLogLikelihood<-function(w,X,y){
	z <- sigmoid(X %*% w)
	return(1/m * sum(-y*log(z) - (1-y)*log(1-z)))
}

loss_df<-data.frame()
lr_list <- 10^c(-1:-5)

for(lr in lr_list){
	w <- matrix(rep(0,11),11,1)
    loss <- rep(0,100)
    
    for(i in 1:100)
    {
		y_pre <- sigmoid(X_train %*% w)
		w <- w + lr * 1/m * t(X_train) %*% (y_train - y_pre)
		loss[i] <- getNegativeLogLikelihood(w,X_train,y_train)
    }
   loss_df <- rbind(loss_df,data.frame(1:50,loss[1:50],rep(lr,50)))
}
colnames(loss_df) <-c("Iteration","NLL","StepSize")
loss_df$StepSize<-as.factor(loss_df$StepSize)

p<-ggplot(loss_df,aes(x=Iteration,y=NLL,group=StepSize))+geom_line(aes(color=StepSize))
p<-p + scale_color_brewer(palette="Dark2")
p<-p + labs(x="Number of Iterations",title="NLL-Loss against iterations under different step sizes")
ggsave("Gradient Descent Plot.png")

library(numDeriv)

w <- matrix(rep(0,11),11,1)
newton_loss <- rep(0,100)
for(i in 1:100){
	y_pre <- sigmoid(X_train %*% w)
	grad <- (1/m) * t(X_train) %*% (y_pre - y_train)
	H <- hessian(getNegativeLogLikelihood, w, method = "Richardson", X = X_train, y = y_train)
	w <- w - ginv(H) %*% grad
	newton_loss[i] <- getNegativeLogLikelihood(w,X_train,y_train)
}
newton_loss_df <- data.frame(1:100,newton_loss)
colnames(newton_loss_df)[1]<-"Iteration"
p<-ggplot(newton_loss_df,aes(x=Iteration,y=newton_loss))+geom_line()
p<-p + labs(x="Number of Iterations",y="NLL",title="NLL-Loss against iterations using Newton's method")
ggsave("Newton's Method Plot.png")

'''
> w
            [,1]
 [1,] -1.0882920
 [2,] -0.3738480  -> SNP 1
 [3,] -1.1349158
 [4,]  0.4249776
 [5,] -4.0561327
 [6,] -1.5191263
 [7,]  1.7831843
 [8,]  0.8848188  -> SNP 7
 [9,]  0.9040501
[10,] -0.2883710
[11,]  2.7223962
'''

#png(paste(lr,'-gd.png'))
#plot(loss[1:50],main="Loss(NLL) against Iterations",xlab="Iteration",ylab="Negative Log-Likelihood")
#dev.off()
