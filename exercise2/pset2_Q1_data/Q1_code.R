X_train<-t(read.fwf("Q1.training.geno", widths = rep(1,900)))
y_train<-read.table("Q1.training.pheno")
X_train<-as.matrix(X_train)
y_train<-unlist(y_train)

X_test<-t(read.fwf("Q1.test.geno", widths = rep(1,100)))
y_test<-read.table("Q1.test.pheno")
X_test<-as.matrix(X_test)
y_test<-unlist(y_test)

library(glmnet)
lambdas<-c(0,2,5,8)

for (lambda in lambdas){
	print(paste("lambda:",lambda))
	ridge_fit<-glmnet(X_train,y_train,alpha = 0,lambda = lambda)
	
	y_predicted<-predict(ridge_fit, s = lambda, newx = X_train)
	train_mse <- mean((y_predicted - y_train)^2)  # train MSE 
	print(paste("train_mse:",train_mse))
	y_predicted<-predict(ridge_fit, newx = X_test)
	test_mse <- mean((y_predicted - y_test)^2)  # test MSE 
	print(paste("test_mse:",test_mse))
}

lambdas<-c(2,5,8)
regularization_lambda<-vector()
error<-vector()
for (lambda in lambdas){
	print(paste("lambda:",lambda))
	ridge_fit<-glmnet(X_train,y_train,alpha = 0,lambda = lambda)
	
	regularization_lambda<-c(regularization_lambda,lambda)
	y_predicted<-predict(ridge_fit, s = lambda, newx = X_train)
	train_mse <- mean((y_predicted - y_train)^2)  # train MSE 
	print(paste("train_mse:",train_mse))
	error<-c(error,train_mse)
	
	regularization_lambda<-c(regularization_lambda,lambda)
	y_predicted<-predict(ridge_fit, newx = X_test)
	test_mse <- mean((y_predicted - y_test)^2)  # test MSE 
	print(paste("test_mse:",test_mse))
	error<-c(error,test_mse)
}

set<-rep(c("train","test"),3)
data<-data.frame(regularization_lambda,set,error)

library(ggplot2)
data$set <- factor(data$set,levels=c("train","test"))
p <- ggplot(data=data,aes(x=regularization_lambda,y=error,fill=set))+geom_bar(stat="identity", position=position_dodge())
p <- p +labs(y="MSE",fill="dataset",title = "The mean squared error (MSE) against different regularization parameters")
p <- p + scale_x_discrete(limits=lambdas)
ggsave("MSE.png") 

'''
[1] "lambda: 0"
[1] "train_mse: 0.153748302437673"
[1] "test_mse: 375.025527566054"
[1] "lambda: 2"
[1] "train_mse: 5.97855832096759"
[1] "test_mse: 193.29291235802"
[1] "lambda: 5"
[1] "train_mse: 17.9557796556861"
[1] "test_mse: 199.642852505887"
[1] "lambda: 8"
[1] "train_mse: 31.0965331863536"
[1] "test_mse: 212.660520330567"
'''