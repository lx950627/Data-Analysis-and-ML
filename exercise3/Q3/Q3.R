library(MCMCpack)

K<-2
N<-1000
M<-5000
X<-t(read.fwf("mixture1.geno", widths = rep(1,N)))

LogLikelihood <- function(){
	LL <- 0
    
    for(i in 1:N)
    {
    	max_likelihood <- -Inf
		for(k in 1:K){
			likelihood <- log(pie[k])
			for(j in 1:M){
				p <- f[j,k]
				x <- X[i,j]
				likelihood <- likelihood + as.numeric(x * log(p) + (1-x) * log(1-p))
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
				likelihood <- likelihood + as.numeric(x * log(p) + (1-x) * log(1-p))
			}
			
			k_likelihood <- k_likelihood + exp(likelihood - max_likelihood)
			
			#print(paste("likelihood:",likelihood))
			#print(paste("k_likelihood:",k_likelihood))
			
		}
        #print(paste("k_likelihood",k_likelihood))
		LL <- LL + log(k_likelihood) + max_likelihood
    }

	return(LL)
}

Q_Value <- function(){
    Qvalue <- 0

	for(i in 1:N)
	{
		for(k in 1:K)
		{
		   sum_result <- log(pie[k])
           for(j in 1:M)
           {
           	  p <- f[j,k]
           	  x <- X[i,j]
           	  sum_result <- sum_result + as.numeric(x * log(p) + (1-x) * log(1-p))
		   }

		   Qvalue <- Qvalue + r[i,k] * sum_result
      
		}
	}

	return(Qvalue)
}


LL <- 0
Q <- 0
f <- matrix(runif(M*K),nrow=M,ncol=K)
pie <- rdirichlet(1, c(1,1))
r <- array(0,dim=c(N,K))
iteration <- 0
change <- 1

while(iteration < 100){
	# E-step
	
	for(i in 1:N){
        max_likelihood <- -Inf
		for(k in 1:K){
			likelihood <- log(pie[k])
			for(j in 1:M){
				p <- f[j,k]
				x <- X[i,j]
				likelihood <- likelihood + as.numeric(x * log(p) + (1-x) * log(1-p))
			}
			max_likelihood <- max(max_likelihood,likelihood)
		}


		for(k in 1:K){
			numerator <- log(pie[k])
			for(j in 1:M){
				p <- f[j,k]
				x <- X[i,j]
				numerator <- numerator + as.numeric(x * log(p) + (1-x) * log(1-p))
			}

			denominator <- 0 
			for(kk in 1:K){
				likelihood <- log(pie[kk])
				for(j in 1:M){
					p <- f[j,kk]
					x <- X[i,j]
					likelihood <- likelihood + as.numeric(x * log(p) + (1-x) * log(1-p))
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
    Q[iteration] <- Q_Value()

    print(paste("Iteration",iteration,"finishes"))
    print(pie)
    print(paste("Current LogLikelihood:",LL[iteration]))
    print(paste("Current Q value:",Q[iteration]))
    
   # if(iteration >= 2){
    #	change <- LL[iteration] - LL[iteration-1]
    #}   
}