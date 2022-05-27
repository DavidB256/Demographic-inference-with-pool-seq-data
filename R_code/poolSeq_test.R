install.packages("C:/Users/David/Documents/poolSeq-0.3.5.tar.gz", repos=NULL, type="source")

library(poolSeq)
library(genomalicious)

# generate random allele frequencies
af <- runif(10000)

# introduce sampling variance to mimic Pool-seq of the entire population at 100X coverage
afSeq <- sample.alleles(af, size=100, mode="coverage")

# plot distribution of differences in allele frequency before and after sampling
hist(af-afSeq$p.smpld, main="Sequencing at 100X", xlab="Error in allele frequency (%)", ylab="Occurrences")
head(afSeq)
