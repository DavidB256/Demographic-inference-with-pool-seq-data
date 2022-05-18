library(ggplot2)
library(scales)
library(tidyverse)

setwd("C:\\Users\\David\\Desktop\\Bergland\\Default R wd")
data <- read.table("combined_pi_output.txt", header=TRUE)
data$rr_exp <- as.factor(data$rr_exp)
data$mu <- as.factor(data$mu)
plot_line <- ggplot(data, aes(x=BIN_END, y=log(PI), color=rr_exp, shape=mu)) + geom_line() + geom_point()
plot_line

#####
plot_box <- ggplot(data, aes(x=mu, fill=rr_exp, y=log(PI))) + 
  labs(title= "Mutation rate versus windowed pi values for different recombination rates",
       fill = expression(rho)) +
  scale_fill_manual(labels = c(expression(10 ^ -6), expression(10 ^ -7),
                               expression(10 ^ -8), expression(10 ^ -9),
                               expression(10 ^ -10)),
                    values = hue_pal()(5)) +
  geom_boxplot(width = 0.7) +
  xlab(expression(log[10] (-mu))) +
  ylab(log[10] ~ pi) +
  theme_bw()
plot_box

lm(data = data,
   PI ~ as.numeric(mu) + as.numeric(rr_exp) ) %>%
  summary()

