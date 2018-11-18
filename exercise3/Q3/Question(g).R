LogLikelihood <- function(pie,f,K){
	LL <- 0
	small <- 1e-5
    for(i in 1:N)
    {
    	max_likelihood <- -Inf
    	small <- 1e-5
		for(k in 1:K){
			likelihood <- log(pie[k])
			for(j in 1:M){
				p <- f[j,k]
				x <- X[i,j]
				likelihood <- likelihood + as.numeric(x * log(p + small) + (1-x) * log(1-p + small))
			}
			max_likelihood <- max(max_likelihood,likelihood)
		}

    	k_likelihood <- 0 
		for(k in 1:K)
		{
			likelihood <- log(pie[k])
			for(j in 1:M)
			{
				p <- f[j,k]
				x <- X[i,j]
				likelihood <- likelihood + as.numeric(x * log(p + small) + (1-x) * log(1-p + small))
			}
			
			k_likelihood <- k_likelihood + exp(likelihood - max_likelihood)
			
		}
        
		LL <- LL + log(k_likelihood) + max_likelihood
    }

	return(LL)
}

EM <- function(K){ # return final log-likelihood
	LL <- 0
	f <- matrix(runif(M*K),nrow=M,ncol=K)
	pie <- rdirichlet(1, rep(1,K))
	r <- array(0,dim=c(N,K))
	iteration <- 1
	change <- 1
	LL[1] <- LogLikelihood(pie,f,K)
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
		LL[iteration] <- LogLikelihood(pie,f,K)
	    #Q[iteration] <- Q_Value()

	    print(paste("Iteration",iteration-1,"finishes"))
	    print(pie)
	    print(paste("Current LogLikelihood:",LL[iteration]))
	    #print(paste("Current Q value:",Q[iteration]))
	   
	    change <- LL[iteration] - LL[iteration-1]
    
	}

	
	print(paste("K",K,"finishes"))
	print("")
	return(LL[iteration]) 
}

k_list<-1:4
log_likelihood<-rep(0,4)
i<-1
for(k in k_list){
  log_likelihood[i] <- EM(k)
  i <- i+1
}

plotdata<-data.frame(k_list,log_likelihood)
p <- ggplot(plotdata,aes(x=k_list,y=log_likelihood,group=1)) + geom_line(linetype="dotted") + geom_point()
p + labs(x="number of mixture components (K)",y="Log-Likelihood",title="Log-Likelihood against number of mixture components (K)") 
ggsave("g.png")


