library(tidyverse)
data <- read.table("hw.4-1.txt",col.names=c("X5","X4","X3"))

d <-mean(with(data,(X5-X4)*X3)) #0.01075889

pnorm(d, mean = 0, sd = 4e-3, lower.tail = FALSE) * 2  # 0.007151155
