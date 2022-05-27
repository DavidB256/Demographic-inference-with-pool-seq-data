require(vcfR)
require(adegenet)
require(FactoMineR)
require(tidyverse)
require(gridExtra)

setwd("C:\\Users\\David\\Desktop\\Bergland\\data")

# /scratch/djb3ve/data/2island_0mig_model.vcf:    0.092145
# /scratch/djb3ve/data/2island_1mig_model.vcf:    0.030151
# /scratch/djb3ve/data/2island_2mig_model.vcf:    0.045646
# /scratch/djb3ve/data/2island_3mig_model.vcf:    0.577786
# /scratch/djb3ve/data/2island_4mig_model.vcf:    0.656413
# /scratch/djb3ve/data/2island_5mig_model.vcf:    0.941419

x <- data.frame(mig_rate=0:5, fst=c(0.092145, 0.030151, 0.045646, 0.577786,
                               0.656413, 0.941419))
plot(x)
ggplot(x, aes(x=mig_rate, y=fst)) + geom_line() + geom_point()
