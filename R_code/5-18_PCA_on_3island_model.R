require(vcfR)
require(adegenet)
require(FactoMineR)
require(tidyverse)
require(gridExtra)

setwd("C:\\Users\\David\\Desktop\\Bergland\\data")

# Load in and perform PCA on VCF files
get_pca_coords <- function(vcf_name) {
  vcf <- read.vcfR(vcf_name)
  my_genind <- vcfR2genind(vcf)
  vcf_table <- as.data.frame(my_genind@tab)
  vcf_table <- vcf_table[,grep(".0", names(vcf_table))]
  pca_fig <- PCA(vcf_table, graph=FALSE, scale.unit=F, ncp=100)
  pca_fig$ind$coord %>%
    as.data.frame() %>% 
    mutate(pop.lab = rep(0:2, each = 10))
}

coords1 <- get_pca_coords("5-18 3islands VCFs\\3islands_1mig.vcf")
coords2 <- get_pca_coords("5-18 3islands VCFs\\3islands_2mig.vcf")
coords3 <- get_pca_coords("5-18 3islands VCFs\\3islands_3mig.vcf")
coords4 <- get_pca_coords("5-18 3islands VCFs\\3islands_4mig.vcf")

# Plot sample points in VCF files on dominant PCs
get_pca_plot <- function(coords) {
  coords %>%
    ggplot(aes(
      x=Dim.1,
      y=Dim.2,
      color=as.factor(pop.lab)
    )) +
    geom_point( size = 4) +
    geom_density2d() +
    theme(legend.position="none") -> pca_plot
  pca_plot
}

grid.arrange(get_pca_plot(coords1) + ggtitle("migration_rate=0.1"),
             get_pca_plot(coords2) + ggtitle("migration_rate=0.01"),
             get_pca_plot(coords3) + ggtitle("migration_rate=0.001"),
             get_pca_plot(coords4) + ggtitle("migration_rate=0.0001"),
             nrow=2)

# Determine how clustered the subpopulations are along the dominant PCs
get_mean_dist_from_centroid <- function(coords) {
  coords_subpops <- sapply(0:2, function(i, coords) 
    {coords[, c("Dim.1", "Dim.2", "pop.lab")] %>% filter(pop.lab==i)}, coords)
  mean_dists <- apply(coords_subpops, 2, function(coords) {
                      mean(sqrt((coords$Dim.1 - mean(coords$Dim.1))^2 + 
                                  (coords$Dim.2 - mean(coords$Dim.2))^2))})
  mean(mean_dists)
}

coords_vector <- list(coords1, coords2, coords3, coords4)
mean_dists_from_centroids <- sapply(coords_vector, get_mean_dist_from_centroid)
ggplot(data.frame(mig_rate=1:4, dists=as.vector(mean_dists_from_centroids)),
       aes(x=mig_rate, y=dists)) +
  geom_line() +
  geom_point() +
  labs(y="Mean distance from subpopulation centroid along dominant PCs",
       x=expression(-log[10] ~ "(migration rate)"),
       title="Mean distance of individuals from respective subpopulation centroids vs.
migration rate in 3-island population models")










