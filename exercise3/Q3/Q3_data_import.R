library(MCMCpack)

K<-2
N<-1000
M<-5000
X<-t(read.fwf("mixture1.geno", widths = rep(1,N)))


X<-t(read.fwf("mixture2.geno", widths = rep(1,N)))