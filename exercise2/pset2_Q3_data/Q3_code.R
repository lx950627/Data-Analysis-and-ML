m <- 500
X_train<-t(read.fwf("hw.2-2.geno", widths = rep(1,m)))
y_train<-read.table("hw.2-2.pheno")
X_train<-as.matrix(X_train)
y_train<-as.matrix(y_train)

p_thres<-0.05/dim(X_train)[2]
pheno1.model<-lm(y_train[,1]~X_train)
pheno2.model<-lm(y_train[,2]~X_train)
pheno3.model<-lm(y_train[,3]~X_train)
pheno4.model<-lm(y_train[,4]~X_train)

pheno1.plist<-summary(pheno1.model)$coefficients[,4]
pheno2.plist<-summary(pheno2.model)$coefficients[,4]
pheno3.plist<-summary(pheno3.model)$coefficients[,4]
pheno4.plist<-summary(pheno4.model)$coefficients[,4]

pheno1.plist[pheno1.plist<p_thres]
pheno2.plist[pheno2.plist<p_thres]
pheno3.plist[pheno3.plist<p_thres]
pheno4.plist[pheno4.plist<p_thres]

getPList<-function(x){
    plist <- vector()
    for(i in 1:dim(X_train)[2]){
    	model<-lm(y_train[,x]~X_train[,i])
		f <- summary(model)$fstatistic
	    p <- pf(f[1],f[2],f[3],lower.tail=F)
	    plist <- c(plist,p) 
    }
    return(plist)
}

pheno1.plist<-getPList(1)
pheno2.plist<-getPList(2)
pheno3.plist<-getPList(3)
pheno4.plist<-getPList(4)

library(EnvStats)
p_expected<-(-log((1:m)/(m+1)))

p1_observed<-(-log(pheno1.plist))
png("qq1.png")
qqPlot(p_expected,p1_observed,add.line = TRUE,main="phenotype1")
dev.off()

p2_observed<-(-log(pheno2.plist))
png("qq2.png")
qqPlot(p_expected,p2_observed,add.line = TRUE,main="phenotype2")
dev.off()

p3_observed<-(-log(pheno3.plist))
png("qq3.png")
qqPlot(p_expected,p3_observed,add.line = TRUE,main="phenotype3")
dev.off()

p4_observed<-(-log(pheno4.plist))
png("qq4.png")
qqPlot(p_expected,p4_observed,add.line = TRUE,main="phenotype4")
dev.off()

ks.test(pheno1.plist,"punif")
ks.test(pheno2.plist,"punif")
ks.test(pheno3.plist,"punif")
ks.test(pheno4.plist,"punif")

phenotype2_258<-data.frame(X_train[,258],y_train[,2])
colnames(phenotype2_258)<-c("geno","pheno")
phenotype2_258$geno<-as.factor(phenotype2_258$geno)
p <- ggplot(phenotype2_258, aes(x=geno, y=pheno,color=geno)) + geom_boxplot()
p<-p + labs(x="SNP258_Genotype",y="Phenotype2",color="Genotype",title="The Relationship between of Phenotype2 and SNP258")
ggsave("Phenotype2_SNP258.png")

phenotype3_119<-data.frame(X_train[,119],y_train[,3])
colnames(phenotype3_119)<-c("geno","pheno")
phenotype3_119$geno<-as.factor(phenotype3_119$geno)
p <- ggplot(phenotype3_119, aes(x=geno, y=pheno,color=geno)) + geom_boxplot(outlier.size = 0)
p <- p + labs(x="SNP119_Genotype",y="Phenotype3",color="Genotype",title="The Relationship between of Phenotype3 and SNP119")
ggsave("Phenotype3_SNP119.png")

phenotype4_43<-data.frame(X_train[,43],y_train[,4])
colnames(phenotype4_43)<-c("geno","pheno")
phenotype4_43$geno<-as.factor(phenotype4_43$geno)
p <- ggplot(phenotype4_43, aes(x=geno, y=pheno,color=geno)) + geom_boxplot()
p <- p + labs(x="SNP43_Genotype",y="Phenotype4",color="Genotype",title="The Relationship between of Phenotype4 and SNP43")
ggsave("Phenotype4_SNP43.png")

