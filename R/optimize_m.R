#' Optimize the m parameter during denovo stacks assembly
#'
#' This function requires the path to stacks vcf file(s) as input.
#' There are slots for varying the m parameter from 3-7 (as recommended by Paris et al. 2017).
#' After running stacks with each of the m options, plug the output vcf files into this
#' function to calculate the effect of varying m on depth and number of SNPs/loci built. Plug the output of this function into
#' vis_loci() to visualize the optimal the m parameter for your dataset at the 'R80' cutoff (Paris et al. 2017).
#'
#' @param m3 Path to the input vcf file for a run when m=3
#' @param m4 Path to the input vcf file for a run when m=4
#' @param m5 Path to the input vcf file for a run when m=5
#' @param m6 Path to the input vcf file for a run when m=6
#' @param m7 Path to the input vcf file for a run when m=7
#' @return A list containing five summary dataframes, 'depth' showing depth per sample for each m value,
#' 'snp' showing the number of non-missing SNPs retained in each sample at each m value, 'loci' showing the
#' number of non-missing loci retained in each sample at each m value, 'snp.R80' showing the total number of SNPs
#' retained at an 80% completeness cutoff, and 'loci.R80' showing the total number of polymorphic loci
#' retained at an 80% completeness cutoff.
#' @examples
#' optimize_m(m3=system.file("extdata","m3.vcf.gz",package="RADstackshelpR",mustWork=TRUE))
#' @export
optimize_m <- function(m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL){
  #initialize empty depth.df
  depth.df<- data.frame(m=character(), avg.depth=numeric())
  #initialize empty snp.df
  snp.df<- data.frame(var=character(), snps=numeric())
  #initialize empty loci.df
  loci.df<- data.frame(var=character(), loci=numeric())
  #initialize empty snp.80.df
  snp.80.df<- data.frame(var=character(), snps.80=numeric())
  #initialize empty loci.80.df
  loci.80.df<- data.frame(var=character(), loci.80=numeric())
  #set vector of m identifiers
  ms<-c("m3","m4","m5","m6","m7")
  #start on first position in vector of m identifiers
  j=1

  #open for loop for each m identifier
  for(x in list(m3,m4,m5,m6,m7)){
    #open if else statement, if no m of given value, move j up to next m identifier, else calculate snps/loci retained
    if(is.null(x)){j=j+1} else{
      #calculate depth first
      ##read in vcfR
      invisible(utils::capture.output(vcf.r<- vcfR::read.vcfR(x))) #read in all data
      ###calc avg depth of each individual
      dep<- (colSums(vcfR::extract.gt(vcf.r, element='DP', as.numeric=TRUE), na.rm = T)) / (colSums(is.na(vcfR::extract.gt(vcf.r, element='DP', as.numeric=TRUE)) == "FALSE"))
      ###rep m identifier, times = number of samples in the vcf
      m<- rep(ms[j], times = length(dep))
      ###cbind depth and m identifier into depth df
      depth.df<- rbind(depth.df, as.data.frame(cbind(m,dep)))

      #initialize vectors to hold filt level, snps retained, poly loci retained
      snps<- vector("numeric", length = length(dep))
      poly.loci<- vector("numeric", length = length(dep))
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
  out$depth<-depth.df
  out$snp<-snp.df
  out$loci<-loci.df
  out$snp.R80<-snp.80.df
  out$loci.R80<-loci.80.df
  return(out)
  #close function
}
