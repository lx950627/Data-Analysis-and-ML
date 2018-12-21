data <- read.table("ps5_q1.genos",header=T)
y <- read.table("ps5_q1.phenos",header=F)

X_train<-data[1:800,]
y_train<-y[1:800,]
X_test<-data[801:1000,]
y_test<-y[801:1000,]

train_data<-cbind(X_train,y_train)
logitModel<-glm(y_train~.,data=train_data,family=binomial())

pre_prob<-predict(logitModel, X_train,type="response")
pre<-as.numeric(pre_prob>0.5)
mean(pre==y_train)#train accuracy 67.5%

pre_prob<-predict(logitModel, X_test,type="response")
pre<-as.numeric(pre_prob>0.5)
mean(pre==y_test) #test accuracy 65.5%

#nn accuracy 0.9547738693467337

library(tidyverse)
df<-cbind(select(data, snp_4,snp_7),y)
colnames(df)[3]<-"y"

df[which(df$snp_4==1 & df$snp_7==1),]
df[which(df$snp_4==0 & df$snp_7==0),]
df[which(df$snp_4==0 & df$snp_7==1),]
df[which(df$snp_4==1 & df$snp_7==0),]