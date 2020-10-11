#' Optimize the m parameter during denovo stacks assembly
#'
#' This function requires the path to stacks vcf file(s) as input.
#' There are slots for varying the m parameter from 3-7 (as recommended by Paris et al. 2017).
#' After running stacks with each of the m options, plug the output vcf files into this
#' function to visualize the effect of varying m on depth and number of SNPs/loci built to
#' recognize which value optimizes the m parameter for your dataset at the 'R80' cutoff (Paris et al. 2017).
#'
#' @param m3 Path to the input vcf file for a run when m=3
#' @param m4 Path to the input vcf file for a run when m=4
#' @param m5 Path to the input vcf file for a run when m=5
#' @param m6 Path to the input vcf file for a run when m=6
#' @param m7 Path to the input vcf file for a run when m=7
#' @return A list containing two dataframes, 'depth.df' showing depth per sample for each m value,
#' and 'm.df' showing the number of SNPs and loci retained at various filtering levels for each m value
#' @export
optimize_m <- function(m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL){
  #initialize empty depth.df
  depth.df<- data.frame(m=character(), avg.depth=numeric())
  #initialize empty m.df
  m.df<- data.frame(m=character(), filt=as.numeric(), snps=numeric(), V4=character())
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
      vcf.r<- vcfR::read.vcfR(x) #read in all data
      ###calc avg depth of each individual
      dep<- (colSums(vcfR::extract.gt(vcf.r, element='DP', as.numeric=TRUE), na.rm = T)) / (colSums(is.na(vcfR::extract.gt(vcf.r, element='DP', as.numeric=TRUE)) == "FALSE"))
      ###rep m identifier, times = number of samples in the vcf
      m<- rep(ms[j], times = length(dep))
      ###cbind depth and m identifier into depth df
      depth.df<- rbind(depth.df, as.data.frame(cbind(m,dep)))

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
      snpsubset<-as.data.frame(cbind(rep(ms[j], times=11), filt, snps, rep("snp", times = 11)))
      locisubset<-as.data.frame(cbind(rep(ms[j], times=11), filt, poly.loci, rep("loci", times = 11)))
      #match colnames so you can rbind these together in tidy format
      colnames(locisubset)[3]<-"snps"
      #append to existing df
      m.df<- rbind(m.df, as.data.frame(rbind(snpsubset,locisubset)))
      #set j for the next m identifier for next time we go through this loop
      j=j+1
      #close if else statement
    }
    #close for loop
  }

  #take depth df output from all of these possibilities
  #plot hist of depth at each m value on same plot
  depth.df$m<-as.factor(depth.df$m)
  depth.df$dep<-as.numeric(depth.df$dep)
  print("Visualize how different values of m affect average depth in each sample")
  print(
    ggplot2::ggplot(depth.df, ggplot2::aes(x = dep, y = m, fill = m, color = m)) +
      ggridges::geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .35, cex=.5) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "average depth in each sample", y = "m value (minimum stack depth)") +
      ggplot2::theme(legend.position = "none")
  )
  #take m df output from all these possibilities
  #plot number of SNPs retained colored by m at each filt level, as open circles
  #plot the number of polymorphic loci retained colored by m at each filt level, as closed circles
  #plot a vertical line at x=.8
  #return a message as output telling the user to pick the m value with the most polyloci retained at r80 (.8)
  #rename columns
  colnames(m.df)<-c("m","filt","retained","snp.locus")
  m.df$m<-as.character(m.df$m)
  m.df$filt<-as.numeric(as.character(m.df$filt))
  m.df$retained<-as.numeric(as.character(m.df$retained))
  m.df$snp.locus<-as.character(m.df$snp.locus)
  print("The optimal m value returns the most polymorphic loci in the 80% complete matrix (Paris et al. 2017)")
  print(
    ggplot2::ggplot(m.df, ggplot2::aes(x=filt, y=retained, col = m, shape=snp.locus))+
      ggplot2::geom_point(alpha = .75, size=3)+
      ggplot2::ggtitle("total SNPs and polymorphic loci retained by filtering scheme")+
      ggplot2::xlab("fraction of non-missing genotypes required to retain each SNP (0-1)")+
      ggplot2::ylab("# SNPs/loci")+
      ggplot2::theme_light()+
      ggplot2::geom_vline(xintercept=.8)+
      ggplot2::labs(col = "min. stack depth", shape="")
  )

  #return the depth and snp/loci dataframes in case you want to do your own visualizations
  out <- list()
  out$depth<-depth.df
  out$m.comparisons<-m.df
  return(out)
  #close function
}
