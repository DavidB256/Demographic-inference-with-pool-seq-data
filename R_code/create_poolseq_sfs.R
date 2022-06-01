library(poolSeq)
library(vcfR)
library(stringr)

setwd("C:/Users/David/Desktop/Bergland/data")

get_pooled_folded_fs <- function(vcf_name, poolSeq_coverage, popinfo, 
                                 haploid_counts) {
  num_of_pops <- length(haploid_counts)
  # Import VCF file "vcf_name" as "vcf_table"
  vcf <- read.vcfR(vcf_name, verbose=FALSE)
  vcf_genind <- vcfR2genind(vcf)
  vcf_table <- as.data.frame(vcf_genind@tab)
  # Polarize "vcf_table" to remove repeat columns
  polarized_vcf_table <- vcf_table[ , which(sapply(names(vcf_table),
                                                  str_sub, -1, -1) == "1")]
  # Subdivide VCF file by population
  populations <- lapply(0:max(popinfo), 
                        function(i) 
                          { polarized_vcf_table[which(popinfo==i),] })
  allele_counts <- lapply(1:num_of_pops,
                          function(i) {apply(populations[[i]], 2, sum) })
  allele_counts <- lapply(1:num_of_pops,
                          function(i) {sapply(allele_counts[[i]], function(j) 
                                         { min(j, haploid_counts[i] - j) } ) })
  # Apply noise in order to emulate the effects of Pool-seq
  allele_freqs <- lapply(1:num_of_pops, 
                         function(i) { allele_counts[[i]] / haploid_counts[i] } )
  pooled_allele_freqs <- lapply(allele_freqs, 
                                function(x) 
                                  {sample.alleles(x, 
                                                  size=poolSeq_coverage, 
                                                  mode="coverage") } )
  pooled_allele_counts <- lapply(1:num_of_pops,
                                 function(i) 
                                   { round(pooled_allele_freqs[[i]]$p.smpld * haploid_counts[i]) } )
  # Assemble site frequence spectrum
  fs <- array(0, sapply(haploid_counts, function(i) {i + 1}))
  for (i in 1:length(pooled_allele_counts[[1]])) {
    coord <- sapply(pooled_allele_counts, function(x) { x[i] + 1 })
    fs[rbind(coord)] <- fs[rbind(coord)] + 1
  }
  fs
}

vcf_name <- "3island_small_model.vcf"
coverage <- 10
# "popinfo" codes which sample comes from which population. If the ith element
# of "popinfo" is n, then the ith sample is included in the nth population.
# "popinfo" must be 0-indexed and should not skip any numbers
popinfo <- rep(0:2, each=2)
haploid_counts <- c(4, 4, 4)
debug(get_pooled_folded_fs)
undebug(get_pooled_folded_fs)
fs <- get_pooled_folded_fs(vcf_name, coverage, popinfo, haploid_counts)
fs

# Export fs to a text file with one number per line.
# First the length of each dimension in the SFS, then a divider,
# then the values of the SFS
# This gets read into a numpy array for use in moments
write(c(dim(fs), "-----", fs), "test_poolseq_sfs.txt", ncolumns=1)

# Moments output for 3island_small_model.vcf
# 
# [[[-- 51.0 56.0 21.0 44.0]
#   [32.0 0.0 0.0 0.0 0.0]
#   [2.0 0.0 0.0 0.0 0.0]
#   [11.0 0.0 0.0 0.0 --]
#   [115.0 0.0 5.5 -- --]]
#   
#   [[144.0 14.0 1.0 1.0 12.0]
#     [0.0 0.0 0.0 0.0 0.0]
#     [0.0 0.0 0.0 0.0 --]
#     [0.0 0.0 0.0 -- --]
#     [74.0 0.0 -- -- --]]
#   
#   [[268.0 0.0 0.0 153.0 13.0]
#     [0.0 0.0 0.0 0.0 --]
#     [0.0 0.0 0.0 -- --]
#     [0.0 0.0 -- -- --]
#     [13.0 -- -- -- --]]
#   
#   [[5.0 0.0 1.0 0.0 --]
#     [0.0 0.0 0.0 -- --]
#     [0.0 0.0 -- -- --]
#     [0.0 -- -- -- --]
#     [-- -- -- -- --]]
#   
#   [[13.0 0.0 5.5 -- --]
#     [0.0 0.0 -- -- --]
#     [0.0 -- -- -- --]
#     [-- -- -- -- --]
#     [-- -- -- -- --]]]

# Moments output SFS for 2island_1mig_model.vcf
# 
# [[-- 79.0 41.0 14.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#   0.0 0.0 0.0 0.0]
#  [41.0 56.0 29.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#   0.0 0.0 0.0 --]
#  [30.0 46.0 7.0 4.0 0.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#   0.0 0.0 -- --]
#  [1.0 52.0 11.0 2.0 335.0 1.0 33.0 38.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
#   0.0 0.0 -- -- --]
#  [1.0 0.0 7.0 9.0 1.0 1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 --
#   -- -- --]
#  [2.0 5.0 0.0 0.0 0.0 0.0 46.0 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 -- --
#   -- -- --]
#  [0.0 0.0 6.0 0.0 2.0 0.0 14.0 0.0 7.0 1.0 1.0 2.0 0.0 0.0 0.0 -- -- --
#   -- -- --]
#  [0.0 0.0 0.0 3.0 0.0 0.0 3.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 -- -- -- -- --
#   -- --]
#  [0.0 0.0 0.0 0.0 0.0 0.0 0.0 45.0 0.0 0.0 2.0 0.0 0.5 -- -- -- -- -- --
#   -- --]
#  [0.0 0.0 0.0 1.0 3.0 0.0 1.0 0.0 36.0 4.0 0.0 1.5 -- -- -- -- -- -- --
#   -- --]
#  [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 5.0 0.0 -- -- -- -- -- -- -- -- --
#   --]
#  [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.5 -- -- -- -- -- -- -- -- -- --
#   --]
#  [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5 -- -- -- -- -- -- -- -- -- -- -- --]
#  [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -- -- -- -- -- -- -- -- -- -- -- -- --]
#  [0.0 0.0 0.0 0.0 0.0 0.0 0.0 -- -- -- -- -- -- -- -- -- -- -- -- -- --]
#  [0.0 0.0 0.0 0.0 0.0 0.0 -- -- -- -- -- -- -- -- -- -- -- -- -- -- --]
#  [0.0 0.0 0.0 0.0 0.0 -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --]
#  [0.0 0.0 0.0 0.0 -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --]
#  [0.0 0.0 0.0 -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --]
#  [0.0 0.0 -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --]
#  [0.0 -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --]]







