library(stringr)
library(vcfR, quietly=TRUE)
library(data.table)

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
                  # Change all zeros in 'cov' to ones. This prevents division by zero in the next line, but makes the assumption that all sites have
                  # coverage of at least one. Note that sites with zero coverage are extremely rare, especially for larger values of 'size'.
                  cov[which(cov==0)] <- 1
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
# application of pool-seq noise to allele frequencies via the "sample.alleles" function.
# The rounding methods are taken from genomalicious' "dadi_inputs_pools" function
# (Thia, J.A. & Riginos, C. (2019) genomalicious: serving up a smorgasbord of R functions
# for population genomic analyses. bioRxiv. doi: https://doi.org/10.1101/667337).
get_pooled_folded_fs <- function(vcf_name, popinfo, haploid_counts, poolseq_depth, method="counts") {
  num_of_pops <- length(haploid_counts)
  # Import VCF file "vcf_name" as "vcf_table"
  vcf <- read.vcfR(vcf_name, verbose=FALSE)
  vcf_genind <- vcfR2genind(vcf)
  vcf_table <- as.data.frame(vcf_genind@tab)
  # Polarize "vcf_table" to remove repeat columns
  polarized_vcf_table <- vcf_table[, which(lapply(names(vcf_table),
                                                  str_sub, -1, -1) == "0")]
  # Subdivide VCF file by population
  populations <- lapply(unique(popinfo),
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
                                                size=poolseq_depth,
                                                mode="coverage") } )
  # Convert frequencies to counts with one of two different methods.
  if (method == "counts") {
    pooled_allele_counts <- lapply(1:num_of_pops,
                                   function(i)
                                   { lapply(pooled_allele_freqs[[i]]$p.smpld,
                                            function(freq)
                                            { count <- freq * haploid_counts[i]
                                            if (count > 0 && count < 0) 1 else round(count)
                                            }) } )
  } else if (method == "probs") {
    pooled_allele_counts <- lapply(1:num_of_pops,
                                   function(i)
                                   { return(rbinom(n=length(pooled_allele_freqs[[i]]$p.smpld),
                                                   size=haploid_counts[i],
                                                   prob=pooled_allele_freqs[[i]]$p.smpld)) } )
  } else {
    stop("Error: 'method' argument must be either 'counts' or 'probs'.", call.=FALSE)
  }
  # Assemble site frequence spectrum
  fs <- array(0, lapply(haploid_counts, function(i) {i + 1}))
  for (i in 1:length(pooled_allele_counts[[1]])) {
    coord <- sapply(pooled_allele_counts, function(x) { x[[i]] + 1 })
    fs[rbind(coord)] <- fs[rbind(coord)] + 1
    rm(coord)
    gc()
  }

  return(fs)
}

generate_pooled_sfs_from_pipeline_instruction <-
  function(vcf_file, popinfo_eval_string, haploid_counts_eval_string, poolseq_depth,
           iterations, wd_string,
           seed) {
  message(seed)
  set.seed(seed)
  fs <- get_pooled_folded_fs(paste(wd_string, "vcfs/", vcf_file, sep=""),
                             eval(parse(text=popinfo_eval_string)),
                             eval(parse(text=haploid_counts_eval_string)),
                             poolseq_depth)
  output_file_name <- paste(wd_string, "serialized_pooled_sfss/",
                            str_sub(vcf_file, end=-5),
                            "_depth", poolseq_depth,
                            "_seed", seed,
                            "_pooled_sfs_serialized.txt", sep="")
  write(c(dim(fs), "-----", rev(fs)), output_file_name, ncolumns=1)
}

# REQUIRES SUBDIRECTORY "serialized_pooled_sfss/" IN WORKING DIRECTORY

# Operational variables
wd_string <- "/scratch/djb3ve/data/second_pipeline/"
setwd(wd_string)
iterations <- 3

# Handle command line arguments
args <- commandArgs(trailingOnly=TRUE)
vcf_file <- args[1]
popinfo_eval_string <- args[2]
haploid_counts_eval_string <- args[3]
poolseq_depth <- as.numeric(args[4])

lapply(1:iterations, function(seed)
  {generate_pooled_sfs_from_pipeline_instruction(vcf_file, popinfo_eval_string,
                                                 haploid_counts_eval_string, poolseq_depth,
                                                 iterations, wd_string, seed)})
