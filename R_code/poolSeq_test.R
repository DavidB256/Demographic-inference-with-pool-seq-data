library(poolSeq)
library(vcfR)
library(stringr)

setwd("C:/Users/David/Desktop/Bergland/data")

get_folded_fs_1d <- function(vcf_name, haploid_count) {
  vcf <- read.vcfR(vcf_name, verbose=FALSE)
  vcf_genind <- vcfR2genind(vcf)
  vcf_table <- as.data.frame(vcf_genind@tab)
  polarized_vcf_table <- vcf_table[ ,which(sapply(names(vcf_table),
                                                  str_sub, -1, -1) == "1")]
  allele_counts <- apply(polarized_vcf_table, 2, sum)
  fs <- sapply(0:haploid_count, 
               function(i) { count(allele_counts, value=i) })
  folded_fs <- sapply(0:haploid_count,
                      function(i) { ifelse(i <= (haploid_count / 2) + 0.5,
                                           fs[i] + fs[-i],
                                           0 })
  folded_fs
}

get_pooled_fs_1d <- 
  function(vcf_name, poolSeq_coverage, poolSeq_mode, haploid_count) {
  vcf <- read.vcfR(vcf_name, verbose=FALSE)
  vcf_genind <- vcfR2genind(vcf)
  vcf_table <- as.data.frame(vcf_genind@tab)
  polarized_vcf_table <- vcf_table[ ,which(sapply(names(vcf_table),
                                                  str_sub, -1, -1) == "1")]
  allele_freqs <- apply(polarized_vcf_table, 2, sum) / haploid_count
  pooled_allele_freqs <- sample.alleles(allele_freqs, 
                                        size=poolSeq_coverage, 
                                        mode=poolSeq_mode)
  pooled_allele_counts <- round(pooled_allele_freqs$p.smpld * 
                                  pooled_allele_freqs$size / haploid_count)
  fs <- sapply(0:haploid_count, 
               function(i) { count(pooled_allele_counts, value=i) })
}

pooled_fs <- get_pooled_fs_1d("small_vcf.vcf", 100, "coverage", 10)
pooled_fs
sum(fs)

fs <- get_fs_1d("small_vcf.vcf", 10)
fs
sum(fs)




