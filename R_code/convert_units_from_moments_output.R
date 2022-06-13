#Import data
setwd("C:/Users/David/Downloads/Cyberduck Downloads")
data01 <- read.table("bounds01_optimizelog.txt", header=TRUE)
data1 <- read.table("bounds1_optimizelog.txt", header=TRUE)
data10 <- read.table("bounds10_optimizelog.txt", header=TRUE)
data100 <- read.table("bounds100_optimizelog.txt", header=TRUE)

# Create merged data frame
data_combined <- rbind(data01, data1, data10, data100)
data_combined["bound"] <- rep(c(0.1, 1, 10, 100), each=5)
data_combined

# Convert units from moments units into msprime units
data_combined_adj <- data_combined[,5:11]
data_combined_adj[,1:2] <- data_combined_adj[,1:2] * data_combined$theta / (4 * 2.9e-6 * 1e5)
data_combined_adj[,3] <- data_combined_adj[,3] * data_combined$theta / (2 * 2.9e-6 * 1e5)
data_combined_adj[,4] <- data_combined_adj[,4] * (2 * 2.9e-6 * 1e5) / data_combined$theta
data_combined_adj

# Select only the rows with the highest log-likelihoods
data_combined_adj_high_ll <- data_combined_adj[which(data_combined_adj$ll > -400),]
data_combined_adj_high_ll
write.csv(data_combined_adj_high_ll, "first_successful_moments_inference_two_pop_split_model.csv", row.names=FALSE)






