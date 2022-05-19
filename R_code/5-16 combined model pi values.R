library(ggplot2)

setwd("C:\\Users\\David\\Desktop\\Bergland\\R")
data <- read.table("5-16 combined_models.pi", header=TRUE)
plot_line <- ggplot(data, aes(x=BIN_END, y=(PI), color=MODEL)) + geom_line() + geom_point()
plot_line
