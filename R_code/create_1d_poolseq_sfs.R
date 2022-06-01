library(poolSeq)
library(vcfR)
library(stringr)

setwd("C:/Users/David/Desktop/Bergland/data")

get_pooled_folded_fs_1d <- function(vcf_name, poolSeq_coverage, haploid_count) {
  # Import VCF file "vcf_name" as "vcf_table"
  vcf <- read.vcfR(vcf_name, verbose=FALSE)
  vcf_genind <- vcfR2genind(vcf)
  vcf_table <- as.data.frame(vcf_genind@tab)
  # Get counts of minor allele at each SNP locus
  polarized_vcf_table <- vcf_table[ ,which(sapply(names(vcf_table),
                                                  str_sub, -1, -1) == "1")]
  allele_counts <- apply(polarized_vcf_table, 2, sum)
  allele_counts <- sapply(allele_counts, 
                          function(i) { min(i, haploid_count - i) })
  # Apply noise in order to emulate the effects of Pool-seq
  allele_freqs <- allele_counts / haploid_count
  pooled_allele_freqs <- sample.alleles(allele_freqs, 
                                        size=poolSeq_coverage, 
                                        mode="coverage")
  pooled_allele_counts <- round(pooled_allele_freqs$p.smpld * haploid_count)
  # Assemble site frequence spectrum
  fs <- sapply(0:haploid_count, 
               function(i) { matrixStats::count(pooled_allele_counts, value=i) })
  fs
}

pooled_fs <- get_pooled_folded_fs_1d("small_vcf.vcf", 100, 10)
pooled_fs
sum(pooled_fs)

# Obtained with moments
actual_fs <- c(0, 43, 6, 64, 119, 0, 0, 0, 0, 0, 0)

coverages <- c(10, 20, 50, 80, 100, 200, 500, 800, 1000)
pooled_fs_at_diff_coverage <- 
  sapply(coverages, get_pooled_folded_fs_1d,
         vcf_name="small_vcf.vcf", haploid_count=10)
pooled_fs_errors <- apply(pooled_fs_at_diff_coverage, 2,
                         function(x) { dist(rbind(x, actual_fs)) })
data <- data.frame(coverage=coverages, error=pooled_fs_errors)
ggplot(data, aes(x=coverage, y=error)) + geom_point() + geom_line() +
  xlab("Average Pool-seq coverage") + 
  ylab("Euclidean distance between actual SFS and SFS created from Pool-seq")

