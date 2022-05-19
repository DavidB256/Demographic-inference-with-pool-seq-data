require(vcfR)
require(adegenet)
require(FactoMineR)
require(tidyverse)

# Import VCF files
setwd("C:\\Users\\David\\Desktop\\Bergland\\R")
vcf <- read.vcfR("d3_model.vcf")

# Load VCF tables and perform PCA
my_genind <- vcfR2genind(vcf)
vcf_table <- as.data.frame(my_genind@tab)
vcf_table <- vcf_table[,grep(".0", names(vcf_table))]
pca_fig <- PCA(vcf_table, graph=FALSE, scale.unit=F, ncp=100)

plot(pca_fig$eig[, 3])

# PC plot
coords <- as.data.frame(pca_fig$ind$coord)
coords["pop"] = rep(0:8, each=10)
plot(coords[, 1], coords[, 2], col=coords$pop)


vcf_low_mig <- read.vcfR("d3_model_low_mig.vcf")

# Load VCF tables and perform PCA
my_genind2 <- vcfR2genind(vcf_low_mig)
vcf_table2 <- as.data.frame(my_genind2@tab)
vcf_table2 <- vcf_table2[,grep(".0", names(vcf_table2))]
pca_fig2 <- PCA(vcf_table2, graph=FALSE, scale.unit=F, ncp=100)

# PC plot for low mig
coords2 <- as.data.frame(pca_fig2$ind$coord)
coords2["pop"] = rep(0:8, each=10)
plot(coords2[, 1], coords2[, 2], col=coords2$pop)



pca_fig2$ind$coord %>%
  as.data.frame() %>% 
  mutate(pop.lab = rep(0:8, each = 10)) %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    color=as.factor(pop.lab)
  )) +
  geom_point( size = 4) +
  geom_density2d()

# junk
dim(vcf_table)
