LL <- 0
Q <- 0
f <- matrix(runif(M*K),nrow=M,ncol=K)
pie <- rdirichlet(1, c(1,1))
r <- array(0,dim=c(N,K))
iteration <- 1
change <- 1
LL[1] <- LogLikelihood()
print(paste("Initial LogLikelihood:",LL[iteration]))

while(iteration <= 100 && change >= 1e-8){
	# E-step
	small <- 1e-6
	for(i in 1:N){
        max_likelihood <- -Inf
		for(k in 1:K){
			likelihood <- log(pie[k])
			for(j in 1:M){
				p <- f[j,k]
				x <- X[i,j]
				likelihood <- likelihood + as.numeric(x * log(p + small) + (1-x) * log(1-p + small))
			}
			max_likelihood <- max(max_likelihood,likelihood)
		}


		for(k in 1:K){
			numerator <- log(pie[k])
			for(j in 1:M){
				p <- f[j,k]
				x <- X[i,j]
				numerator <- numerator + as.numeric(x * log(p + small) + (1-x) * log(1-p + small))
			}

			denominator <- 0 
			for(kk in 1:K){
				likelihood <- log(pie[kk])
				for(j in 1:M){
					p <- f[j,kk]
					x <- X[i,j]
					likelihood <- likelihood + as.numeric(x * log(p + small) + (1-x) * log(1-p + small))
				}
				
				denominator <- denominator + exp(likelihood - max_likelihood)
			}
			#print(denominator)
			r[i,k] <- exp(numerator - max_likelihood) / denominator
		}
	}


    # M-step
    for(k in 1:K){
    	pie[k] <- sum(r[,k])/ N
    }

    for(j in 1:M){
    	for(k in 1:K){
    		f[j,k] <- sum(r[,k] * X[,j])/sum(r[,k])
    	}
    }

    # Log-Likelihood updates
	iteration <- iteration + 1
	LL[iteration] <- LogLikelihood()
    #Q[iteration] <- Q_Value()

    print(paste("Iteration",iteration-1,"finishes"))
    print(pie)
    print(paste("Current LogLikelihood:",LL[iteration]))
    #print(paste("Current Q value:",Q[iteration]))
   
    change <- LL[iteration] - LL[iteration-1]
    
}

library(ggplot2)
library(ramify)

plotdata<-data.frame(as.factor(0:(iteration-1)),LL)
colnames(plotdata)[1]<-"iteration"
p <- ggplot(plotdata,aes(x=iteration,y=LL,group=1)) + geom_line(linetype="dotted") + geom_point()
p + labs(x="Number of Iteration",y="Log-Likelihood",title="Log-Likelihood against number of iterations in EM") 
ggsave("a.png")


truelabel<-read.table("mixture1.ganc")
truelabel<-argmax(truelabel)
inferlabel<-argmax(r)
mean(truelabel==inferlabel)  #100%