#' Optimize the M parameter during denovo stacks assembly
#'
#' This function requires the path to stacks vcf file(s) as input.
#' There are slots for varying the M parameter from 1-8 (as recommended by Paris et al. 2017).
#' After running stacks with each of the M options, plug the output vcf files into this
#' function to visualize the effect of varying m on the number of SNPs and loci built to
#' recognize which value optimizes the M parameter for your dataset at the 'R80' cutoff (Paris et al. 2017).
#'
#' @param M1 Path to the input vcf file for a run when M=1
#' @param M2 Path to the input vcf file for a run when M=2
#' @param M3 Path to the input vcf file for a run when M=3
#' @param M4 Path to the input vcf file for a run when M=4
#' @param M5 Path to the input vcf file for a run when M=5
#' @param M6 Path to the input vcf file for a run when M=6
#' @param M7 Path to the input vcf file for a run when M=7
#' @param M8 Path to the input vcf file for a run when M=8
#' @return A dataframe showing the number of SNPs and loci retained across filtering levels for each M value
#' @export
optimize_M <- function(M1=NULL,M2=NULL,M3=NULL,M4=NULL,M5=NULL,M6=NULL,M7=NULL,M8=NULL){
  #initialize empty M.df
  M.df<- data.frame(M=character(), filt=as.numeric(), snps=numeric(), V4=character())
  #set vector of m identifiers
  Ms<-c("M1","M2","M3","M4","M5","M6","M7","M8")
  #start on first position in vector of m identifiers
  j=1

  #open for loop for each m identifier
  for(x in list(M1,M2,M3,M4,M5,M6,M7,M8)){
    #open if else statement, if no m of given value, move j up to next m identifier, else calculate snps/loci retained
    if(is.null(x)){j=j+1} else{
      ##read in vcfR
      vcf.r<- vcfR::read.vcfR(x) #read in all data

      #initialize vectors to hold filt level, snps retained, poly loci retained
      filt<- vector("numeric", length = 11)
      snps<- vector("numeric", length = 11)
      poly.loci<- vector("numeric", length = 11)
      ##run loop to fill up vectors with a value for each filter level
      k=1
      for (i in c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)){
        #calculate the completeness cutoff for each snp to be retained
        filt[k]<-i
        #calculate the number of snps retained at this cutoff
        snps[k]<-nrow(vcf.r@gt[(rowSums(is.na(vcf.r@gt))/ncol(vcf.r@gt) <= 1-i),])
        #calculate number of polymorphic loci retained at this cutoff
        poly.loci[k]<-length(unique(vcf.r@fix[,1][(rowSums(is.na(vcf.r@gt))/ncol(vcf.r@gt) <= 1-i)]))
        k=k+1
        #close for loop
      }
      ##cbind these three vectors with m identifier and append it all to the m df
      snpsubset<-as.data.frame(cbind(rep(Ms[j], times=11), filt, snps, rep("snp", times = 11)))
      locisubset<-as.data.frame(cbind(rep(Ms[j], times=11), filt, poly.loci, rep("loci", times = 11)))
      #match colnames so you can rbind these together in tidy format
      colnames(locisubset)[3]<-"snps"
      #append to existing df
      M.df<- rbind(M.df, as.data.frame(rbind(snpsubset,locisubset)))
      #set j for the next m identifier for next time we go through this loop
      j=j+1
      #close if else statement
    }
    #close for loop
  }

  print("Optimal M value returns the most polymorphic loci in the 80% complete matrix (Paris et al. 2017)")
  #take m df output from all these possibilities
  #plot number of SNPs retained colored by m at each filt level, as open circles
  #plot the number of polymorphic loci retained colored by m at each filt level, as closed circles
  #plot a vertical line at x=.8
  #return a message as output telling the user to pick the m value with the most polyloci retained at r80 (.8)
  #rename columns
  colnames(M.df)<-c("M","filt","retained","snp.locus")
  M.df$M<-as.character(M.df$M)
  M.df$filt<-as.numeric(as.character(M.df$filt))
  M.df$retained<-as.numeric(as.character(M.df$retained))
  M.df$snp.locus<-as.character(M.df$snp.locus)
  print(
    ggplot2::ggplot(M.df, ggplot2::aes(x=filt, y=retained, col = M, shape=snp.locus))+
      ggplot2::geom_point(alpha = .75, size=3)+
      ggplot2::ggtitle("total SNPs and polymorphic loci retained by filtering scheme")+
      ggplot2::xlab("fraction of non-missing genotypes required to retain each SNP (0-1)")+
      ggplot2::ylab("# SNPs/loci")+
      ggplot2::theme_light()+
      ggplot2::geom_vline(xintercept=.8)+
      ggplot2::labs(col = c("mismatches allowed\nbetween stacks\nwithin a locus"), shape="")
  )

  print("Correctly setting M requires a balance – set it too low and alleles from the same locus will not collapse, set it too high and paralogous or repetitive loci will incorrectly merge together. When alleles from the same locus are undermerged, the software will incorrectly consider them as independent loci. When loci are overmerged because they happen to be close in sequence space, an errant locus with false polymorphism will result. Therefore, M is particularly dataset‐specific because it depends on the natural levels of polymorphism in the species, as well as the amount of error generated during the preparation and sequencing of the RAD‐seq libraries.")

  #return the depth and snp/loci dataframes in case you want to do your own visualizations
  return(M.df)
    #close function
}
