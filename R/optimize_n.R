#' Optimize the n parameter during denovo stacks assembly
#'
#' This function requires the path to stacks vcf file(s) as input.
#' There are slots for varying the n parameter across M-1, M, and M-1 (as recommended by Paris et al. 2017).
#' After running stacks with each of the n options, plug the output vcf files into this
#' function to visualize the effect of varying m on number of SNPs and loci built to
#' recognize which value optimizes the n parameter for your dataset at the 'R80' cutoff (Paris et al. 2017).
#'
#' @param nequalsMminus1 Path to the input vcf file for a run when n=M-1
#' @param nequalsM Path to the input vcf file for a run when n=M
#' @param nequalsMplus1 Path to the input vcf file for a run when n=M+1
#' @return A dataframe showing the number of SNPs and loci retained across filtering levels for each n value
#' @examples
#' optimize_n(nequalsM =
#' system.file("extdata","nequalsm.vcf.gz",package="RADstackshelpR",mustWork=TRUE))
#' @export
optimize_n <- function(nequalsMminus1=NULL,nequalsM=NULL,nequalsMplus1=NULL){
  #initialize empty snp.df
  snp.df<- data.frame(var=character(), snps=numeric())
  #initialize empty loci.df
  loci.df<- data.frame(var=character(), loci=numeric())
  #initialize empty snp.80.df
  snp.80.df<- data.frame(var=character(), snps.80=numeric())
  #initialize empty loci.80.df
  loci.80.df<- data.frame(var=character(), loci.80=numeric())
  #set vector of m identifiers
  ms<-c("n=M-1","n=M","n=M+1")
  #start on first position in vector of m identifiers
  j=1

  #open for loop for each m identifier
  for(x in list(nequalsMminus1,nequalsM,nequalsMplus1)){
    #open if else statement, if no m of given value, move j up to next m identifier, else calculate snps/loci retained
    if(is.null(x)){j=j+1} else{
      ##read in vcfR
      invisible(utils::capture.output(vcf.r<- vcfR::read.vcfR(x))) #read in all data
      #initialize vectors to hold filt level, snps retained, poly loci retained
      snps<- vector("numeric", length = ncol(vcf.r@gt)-1)
      poly.loci<- vector("numeric", length = ncol(vcf.r@gt)-1)
      ###rep m identifier, times = number of samples in the vcf
      m<- rep(ms[j], times = ncol(vcf.r@gt)-1)
      ##run loop to fill up vectors with a value for each filter level
      k=1
      for (i in 2:ncol(vcf.r@gt)){
        #calculate the number of non-missing SNPs present in the given sample
        snps[k]<-sum(is.na(vcf.r@gt[,i]) == FALSE)
        #calculate number of polymorphic loci present in the given sample
        poly.loci[k]<-length(unique(vcf.r@fix[,1][is.na(vcf.r@gt[,i]) == FALSE]))
        k=k+1
        #close for loop
      }
      #append each to existing df
      snp.df<- rbind(snp.df, as.data.frame(cbind(m, snps)))
      loci.df<- rbind(loci.df, as.data.frame(cbind(m, poly.loci)))

      #calculate the number of loci and SNPs retained in the 80% complete dataset for the given m value
      snps.80<-nrow(vcf.r@gt[(rowSums(is.na(vcf.r@gt))/ncol(vcf.r@gt) <= .2),])
      #calculate number of polymorphic loci retained at this cutoff
      poly.loci.80<-length(unique(vcf.r@fix[,1][(rowSums(is.na(vcf.r@gt))/ncol(vcf.r@gt) <= .2)]))
      #append each to existing df
      snp.80.df<- rbind(snp.80.df, as.data.frame(cbind(ms[j], snps.80)))
      loci.80.df<- rbind(loci.80.df, as.data.frame(cbind(ms[j], poly.loci.80)))

      #set j for the next m identifier for next time we go through this loop
      j=j+1
      #close if else statement
    }
    #close for loop
  }

  #fix colnames
  colnames(snp.80.df)<-c("var","snps.80")
  colnames(loci.80.df)<-c("var","poly.loci.80")
  colnames(snp.df)<-c("var","snps")
  colnames(loci.df)<-c("var","poly.loci")

  #return the depth and snp/loci dataframes in case you want to do your own visualizations
  out <- list()
  out$snp<-snp.df
  out$loci<-loci.df
  out$snp.R80<-snp.80.df
  out$loci.R80<-loci.80.df
  return(out)
  #close function
}
