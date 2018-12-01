library(tidyverse)
library(EnvStats)
data <- data.frame(read_tsv("pset4.q1.data.tsv"))

M <- dim(data)[2] - 1

getPList<-function(){
    plist <- vector()
    for(i in 1:M){
    	model<-lm(data[,M+1] ~ data[,i])
		f <- summary(model)$fstatistic
	    p <- pf(f[1],f[2],f[3],lower.tail=F)
	    plist <- c(plist,p) 
    }
    return(plist)
}
plist <- getPList()

p_expected<-(-log((1:M)/(M+1)))
p_observed<-(-log(plist))
png("qq-1.png")
qqPlot(p_expected,p_observed,qq.line.type="0-1",add.line = TRUE,main="p-value QQ plot")
dev.off()

ks.test(plist,"punif",0,1)

thres <- 0.05 / M
length(plist[plist <= thres]) # 14

pca_result<-prcomp(data[,1:M],center=T,scale.=T)

p <- ggplot(data = data.frame(pca_result$x)) 
p <- p +  geom_point(aes(x = PC1, y = PC2))
p <- p +  labs(title="Population Visualization using PC1 and PC2")
p <- p +  theme_bw() 
ggsave("pca.png")

pc1 <- pca_result$x[,1]

getPList_pc<-function(){
    plist <- vector()
    for(i in 1:M){
    	model<-lm(data[,M+1] ~ data[,i] + pc1)
	    p <- summary(model)$coefficients[2,4]
	    plist <- c(plist,p) 
    }
    return(plist)
}

plist <- getPList_pc()

p_expected<-(-log((1:M)/(M+1)))
p_observed<-(-log(plist))
png("qq-2.png")
qqPlot(p_expected,p_observed,qq.line.type="0-1",add.line = TRUE,main="p-value QQ plot (PC1 included)")
dev.off()

length(plist[plist <= thres]) # 0