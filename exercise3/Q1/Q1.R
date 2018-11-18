library('tidyverse')
set.seed(0)
data <-read_tsv('q1.data')
snps <-data %>% select(starts_with('rs'))

pca_result<-prcomp(snps,center=T,scale.=T)

p <- ggplot(data = data.frame(pca_result$x)) 
p <- p +  geom_point(aes(x = PC1, y = PC2,color=data$population)) 
p <- p +  scale_color_manual(values=c("green", "orange", "blue","red"))
p <- p +  labs(title="Population Visualization using PC1 and PC2",color='population')
p <- p +  theme_bw() 

ggsave("pca.png")

kmeans_result <- kmeans(snps, centers = 4, nstart = 5)
cluster_assign<-kmeans_result$cluster
idx <- sort(kmeans_result$size,index.return=T,decreasing=T)$ix
for(i in 1:length(idx)){
	cluster_assign[which(cluster_assign==idx[i])]<-paste("Cluster",i)
}

p <- ggplot(data = data.frame(pca_result$x)) 
p <- p +  geom_point(aes(x = PC1, y = PC2,color=cluster_assign)) 
p <- p +  scale_color_manual(values=c("red", "blue", "green","orange"))
p <- p +  labs(title="Cluster Visualization using PC1 and PC2",color='cluster')
p <- p +  theme_bw() 

ggsave("cluster.png")

cluster2population<-function(cluster){
	number <- substr(cluster,nchar(cluster),nchar(cluster))
	if(number == '1'){
	   return("EUR")
	}
	else if(number == '2'){
	   return("ASN")
	}
	else if(number == '3'){
	   return("AFR")
	}
	else if(number == '4'){
	   return("AMR")
	}
}

new_assign<-sapply(cluster_assign,cluster2population)
same <- 0
for (i in 1:length(new_assign)){
	if(new_assign[i] == data$population[i]){
		same <- same+1
	}
}
same/length(new_assign) # 94.14%