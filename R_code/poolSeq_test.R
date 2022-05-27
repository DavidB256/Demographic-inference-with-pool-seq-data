library(poolSeq)
library(vcfR)
library(stringr)

setwd("C:/Users/David/Desktop/Bergland/data")

get_pooled_allele_counts <- function(vcf_name, poolSeq_coverage, poolSeq_mode) {
  vcf <- read.vcfR(vcf_name)
  vcf_genind <- vcfR2genind(vcf)
  vcf_table <- as.data.frame(vcf_genind@tab)
  polarized_vcf_table <- vcf_table[ ,which(sapply(names(vcf_table),
                                                  str_sub, -1, -1) == "1")]
  allele_freqs <- apply(polarized_vcf_table, 2, sum) / 10
  pooled_allele_freqs <- sample.alleles(allele_freqs, 
                                        size=poolSeq_coverage, 
                                        mode=poolSeq_mode)
  pooled_allele_counts <- pooled_allele_freqs$p.smpld * pooled_allele_freqs$size
  pooled_allele_counts
}

get_pooled_allele_counts("small_vcf.vcf", 100, "coverage")


