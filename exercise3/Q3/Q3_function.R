LogLikelihood <- function(){
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
    small <- 1e-5

	for(i in 1:N)
	{
		for(k in 1:K)
		{
		   sum_result <- log(pie[k])
           for(j in 1:M)
           {
           	  p <- f[j,k]
           	  x <- X[i,j]
           	  sum_result <- sum_result + as.numeric(x * log(p + small) + (1-x) * log(1-p + small))
		   }

		   Qvalue <- Qvalue + r[i,k] * sum_result
      
		}
	}

	return(Qvalue)
}
