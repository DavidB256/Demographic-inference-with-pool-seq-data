library(stringr)
library(vcfR)

# This function was copied from Thomas Taus' poolSeq R package source code
sample.alleles <- function(p, size, mode=c("coverage", "individuals"), Ncensus=NA, ploidy=2) {
  # determine number of return values
  maxlen <- max(length(p), length(size))
  # check length of a, n and size parameter
  if(maxlen %% length(p) != 0 || maxlen %% length(size) != 0)
    warning("Parameters differ in length and are not multiple of one another: length(p)=", length(p), ", length(size)=", length(size))
  
  # check mode parameter
  mode <- match.arg(mode)
  
  # sample allele counts according to the specified method
  return(switch(mode,
                coverage={
                  # if length of 'size' equals '1' then generate target coverage values using the Poisson distribution, otherwise use values of 'size' directly
                  cov <- if(length(size) == 1) rpois(n=maxlen, lambda=size) else size
                  # sample allele frequencies from Binomial distribution based on 'cov' and 'p'
                  p.smpld <- rbinom(n=maxlen, size=cov, prob=p) / cov
                  # return results, including coverage values if they were drawn from a Poisson distribution
                  if(length(size) == 1) data.table(p.smpld=p.smpld, size=cov) else p.smpld
                },
                individuals={
                  # if length of 'size' is larger than 1, then send warning message
                  if(length(size) > 1)
                    warning("Only the first element in 'size' will be used, because sampling mode is 'individuals'")
                  
                  # sample random allele frequencies from Hypergeometric distribution
                  rhyper(nn=maxlen, m=p*Ncensus*ploidy, n=(1-p)*Ncensus*ploidy, k=size[1]*ploidy)/(size[1]*ploidy)
                }))
}

# This function mimics moments' "Spectrum.from_data_dict" function, but with the
# application of pool-seq noise to allele frequencies via the "sample.alleles" function
get_pooled_folded_fs <- function(vcf_name, popinfo, haploid_counts, poolseq_coverage) {
  num_of_pops <- length(haploid_counts)
  # Import VCF file "vcf_name" as "vcf_table"
  vcf <- read.vcfR(vcf_name, verbose=FALSE)
  vcf_genind <- vcfR2genind(vcf)
  vcf_table <- as.data.frame(vcf_genind@tab)
  # Polarize "vcf_table" to remove repeat columns
  polarized_vcf_table <- vcf_table[ , which(sapply(names(vcf_table),
                                                   str_sub, -1, -1) == "0")]
  # Subdivide VCF file by population
  populations <- lapply(0:max(popinfo), 
                        function(i) 
                        { polarized_vcf_table[which(popinfo==i),] })
  allele_counts <- lapply(1:num_of_pops,
                          function(i) {apply(populations[[i]], 2, sum) })
  # Apply noise in order to emulate the effects of pool-seq with the sample.alleles
  # function from Thomas Taus' poolSeq package
  allele_freqs <- lapply(1:num_of_pops, 
                         function(i) { allele_counts[[i]] / haploid_counts[i] } )
  pooled_allele_freqs <- lapply(allele_freqs, 
                                function(x) 
                                {sample.alleles(x, 
                                                size=poolseq_coverage, 
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

setwd("/scratch/djb3ve/data/first_models/")

# Hands command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) { stop("Error: Five command line arguments must be supplied.", call.=FALSE) }
vcf_name <- as.character(args[1])
popinfo <- eval(parse(text=args[2])) 
haploid_counts <- eval(parse(text=args[3]))
poolseq_coverage <- as.numeric(args[4])

output_file_name <- paste(str_sub(vcf_name, end=-5), "_pooled_sfs_serialized.txt", sep="")

print(getwd())

fs <- get_pooled_folded_fs(vcf_name, popinfo, haploid_counts, poolseq_coverage)
# Serialize SFS for use in Python with moments
setwd("/scratch/djb3ve/data/first_models/serialized_pooled_sfss/")
write(c(dim(fs), "-----", rev(fs)), output_file_name, ncolumns=1)






