require(vcfR)
require(adegenet)
require(FactoMineR)
require(tidyverse)
require(gridExtra)

setwd("C:\\Users\\David\\Desktop\\Bergland\\R data")

# Load in and perform PCA on VCF files
get_pca_fig <- function(vcf_name) {
  vcf <- read.vcfR(vcf_name)
  my_genind <- vcfR2genind(vcf)
  vcf_table <- as.data.frame(my_genind@tab)
  vcf_table <- vcf_table[,grep(".0", names(vcf_table))]
  pca_fig <- PCA(vcf_table, graph=FALSE, scale.unit=F, ncp=100)
  pca_fig
}

pca_fig1 <- get_pca_fig("5-18 3islands VCFs\\3islands_1mig.vcf")
pca_fig2 <- get_pca_fig("5-18 3islands VCFs\\3islands_2mig.vcf")
pca_fig3 <- get_pca_fig("5-18 3islands VCFs\\3islands_3mig.vcf")
pca_fig4 <- get_pca_fig("5-18 3islands VCFs\\3islands_4mig.vcf")

# Plot sample points in VCF files on dominant PCs
get_pca_plot <- function(pca_fig) {
  pca_fig$ind$coord %>%
    as.data.frame() %>% 
    mutate(pop.lab = rep(0:2, each = 10)) %>%
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

grid.arrange(get_pca_plot(pca_fig1) + ggtitle("migration_rate=0.1"),
             get_pca_plot(pca_fig2) + ggtitle("migration_rate=0.01"),
             get_pca_plot(pca_fig3) + ggtitle("migration_rate=0.001"),
             get_pca_plot(pca_fig4) + ggtitle("migration_rate=0.0001"),
             nrow=2)

