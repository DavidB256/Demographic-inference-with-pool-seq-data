library(tidyverse)

ns <- c(10, 50, 100, 200, 400)
ds <- c(10, 40, 70, 100)
data <- merge(ns, ds)
colnames(data) <- c("n", "d")
data["n_eff"] <- (data$d * data$n) / (data$d + data$n)
data["loss"] <- (data$n - data$n_eff) / data$n %>% abs
data

ggplot(data, aes(x=as.factor(n), y=as.factor(d), fill=loss)) + geom_tile() +
  xlab("Number of haploids sampled") +
  ylab("Pool-seq depth") +
  labs(fill="Loss of effective sample size") +
  ggtitle("Loss of effective sample size from pool-seq depth")
