
args = commandArgs(trailingOnly=TRUE)

# This function was copied from Thomas Taus' poolSeq R package source code
sample.alleles <- function(p, size, mode=c("coverage", "individuals"), Ncensus=NA, ploidy=2) {
  # determin number of return values
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


