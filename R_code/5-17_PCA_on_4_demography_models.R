require(vcfR)
require(adegenet)
require(FactoMineR)

# Import VCF files
setwd("C:\\Users\\David\\Desktop\\Bergland\\R")
vcf1 <- read.vcfR("5-16 model VCFs\\control.vcf")
vcf2 <- read.vcfR("5-16 model VCFs\\continent_islands.vcf")
vcf3 <- read.vcfR("5-16 model VCFs\\equal_islands.vcf")
vcf4 <- read.vcfR("5-16 model VCFs\\stepping_stone.vcf")

# Load VCF tables and perform PCA
my_genind1 <- vcfR2genind(vcf1)
vcf_table1 <- as.data.frame(my_genind1@tab)
vcf_table1 <- vcf_table1[,grep(".0", names(vcf_table1))]
pca_fig1 <- PCA(vcf_table1, graph=FALSE, scale.unit=F, ncp=100)

my_genind2 <- vcfR2genind(vcf2)
vcf_table2 <- as.data.frame(my_genind2@tab)
vcf_table2 <- vcf_table2[,grep(".0", names(vcf_table2))]
pca_fig2 <- PCA(vcf_table2, graph=FALSE, scale.unit=F, ncp=100)

my_genind3 <- vcfR2genind(vcf3)
vcf_table3 <- as.data.frame(my_genind3@tab)
vcf_table3 <- vcf_table3[,grep(".0", names(vcf_table3))]
pca_fig3 <- PCA(vcf_table3, graph=FALSE, scale.unit=F, ncp=100)

my_genind4 <- vcfR2genind(vcf4)
vcf_table4 <- as.data.frame(my_genind4@tab)
vcf_table4 <- vcf_table4[,grep(".0", names(vcf_table4))]
pca_fig4 <- PCA(vcf_table4, graph=FALSE, scale.unit=F, ncp=100)

# Plot cumulative variance explanation by PC factors
plot(pca_fig1$eig[, 3])

# PC plot 2
coords2 <- as.data.frame(pca_fig2$ind$coord)
coords2["pop"] = rep(c(0, 1, 2, 3), each=100)
plot(coords2[, 1], coords2[, 2], col=coords2$pop)

# PC plot 3
coords3 <- as.data.frame(pca_fig3$ind$coord)
coords3["pop"] = rep(c(0, 1, 2, 3), each=100)
plot(coords3[, 1], coords3[, 2], col=coords3$pop)

# PC plot 4
coords4 <- as.data.frame(pca_fig4$ind$coord)
coords4["pop"] = rep(c(0, 1, 2, 3), each=100)
plot(coords4[, 1], coords4[, 2], col=coords4$pop)

# Junk
pca_fig1
vcf_table1[1:10, 1:10]
dim(vcf_table1)





