
#' @export
missing.by.snp <- function(vcfR, cutoff=NULL){
  
  if (!is.null(cutoff)) {
    
    #do basic vis to show cutoff
    #extract genotype from the vcf
    dp.matrix<- vcfR::extract.gt(vcfR, as.numeric=TRUE)
    
    #calculate the proportion of individuals successfully genotyped at each SNP
    miss<-rowSums(is.na(dp.matrix))/ncol(dp.matrix)
    
    #loop that stores a vector of # non-missing loci retained for each individual
    #looped over all possible completeness filter values
    #initialize df.x
    df.x<- data.frame(filt=numeric(), missingness=numeric(), snps.retained=numeric())
    #loop
    for (i in c(.3,.5,.6,.65,.7,.75,.8,.85,.9,.95,1)){
      #subset the matrix to only snps passing a given filter level
      gen.z.mat<-dp.matrix[1-miss >= i,]
      #calc the total number of snps retained at the given filter level
      snps.retained<-nrow(gen.z.mat)
      #calculate the total missing data % for the given cutoff
      missingness<-sum(is.na(gen.z.mat))/(ncol(gen.z.mat)*nrow(gen.z.mat))
      #calculate the completeness cutoff for each snp to be retained
      filt<-i
      df.x<-rbind(df.x, as.data.frame(cbind(filt, missingness, snps.retained)))
      }
    
    #make columns numeric for plotting
    df.x$filt<-as.numeric(as.character(df.x$filt))
    df.x$missingness<-as.numeric(as.character(df.x$missingness))
    df.x$snps.retained<-as.numeric(as.character(df.x$snps.retained))
    
    #visualize dotplot for total loci retained at each filter level
    plot1<-ggplot2::ggplot(df.x, ggplot2::aes(x=filt)) +
      ggplot2::scale_x_continuous(breaks=c(.3,.5,.6,.65,.7,.75,.8,.85,.9,.95,1))+
      ggplot2::geom_point(ggplot2::aes(y=snps.retained)) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "SNP completeness cutoff", y = "total SNPs retained")+
      ggplot2::geom_vline(xintercept = cutoff, color = "red")
    
    #visualize dotplot for missing data % at each filter level
    plot2<-ggplot2::ggplot(df.x, ggplot2::aes(x=filt)) +
      ggplot2::scale_x_continuous(breaks=c(.3,.5,.6,.65,.7,.75,.8,.85,.9,.95,1))+
      ggplot2::geom_point(ggplot2::aes(y=missingness)) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "SNP completeness cutoff", y = "total % missing data")+
      ggplot2::geom_vline(xintercept = cutoff, color = "red")
    
    gridExtra::grid.arrange(plot1,plot2)
    
    #calc # of SNPs filtered
    p<-round(sum(miss > 1-cutoff)/length(miss)*100, 2)
    
    #report # of SNPs filtered
    print(paste0(p,"% of SNPs fell below a completeness cutoff of ", cutoff, " and were removed from the VCF"))
    
    #do filtering
    vcfR <- vcfR[miss <= 1-cutoff, ]
    
    return(vcfR)
    
  } 
  else {
    
    ###Part 1
    #Vis to understand the interplay between retaining samples and setting a missing data cutoff
    #extract genotype from the vcf
    dp.matrix<- vcfR::extract.gt(vcfR, as.numeric=TRUE)
    
    #calculate vector containing the proportion of individuals successfully genotyped at each SNP
    miss<-rowSums(is.na(dp.matrix))/ncol(dp.matrix)
    #initialize df.x
    df.y<- data.frame(filt=numeric(), snps=numeric())
    #loop
    for (i in c(.3,.5,.6,.65,.7,.75,.8,.85,.9,.95,1)){
      #subset matrix to only SNPs retained at given filter level
      gen.z.mat<-dp.matrix[1-miss >= i,]
      #calc the total number of snps retained at the given filter level in each sample
      snps<-colSums(is.na(gen.z.mat) == TRUE)/nrow(gen.z.mat)
      #calculate the completeness cutoff for each snp to be retained
      filt<-rep(i, times= length(snps))
      df.y<-rbind(df.y, as.data.frame(cbind(filt, snps)))
      }
    
    #make columns numeric for plotting
    df.y$filt<-as.numeric(as.character(df.y$filt))
    df.y$snps<-as.numeric(as.character(df.y$snps))
    
    #visualize filtering levels as stacked histograms
    print(
      ggplot2::ggplot(df.y, ggplot2::aes(x = snps, y = as.character(filt), fill = filt, color = filt)) +
        ggridges::geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .25, cex=.5) +
        ggplot2::theme_classic() +
        ggplot2::labs(x = "missing data proportion in each sample", y = "SNP completeness cutoff") +
        ggplot2::theme(legend.position = "none")
      )
    
    ###Part 2
    #Vis to make decision on cutoff
    #calculate the proportion of individuals successfully genotyped at each SNP
    miss<-rowSums(is.na(dp.matrix))/ncol(dp.matrix)
    
    #loop that stores a vector of # non-missing loci retained for each individual
    #looped over all possible completeness filter values
    #initialize df.x
    df.x<- data.frame(filt=numeric(), missingness=numeric(), snps.retained=numeric())
    #loop
    for (i in c(.3,.5,.6,.65,.7,.75,.8,.85,.9,.95,1)){
      #subset the matrix to only snps passing a given filter level
      gen.z.mat<-dp.matrix[1-miss >= i,]
      #calc the total number of snps retained at the given filter level
      snps.retained<-nrow(gen.z.mat)
      #calculate the total missing data % for the given cutoff
      missingness<-sum(is.na(gen.z.mat))/(ncol(gen.z.mat)*nrow(gen.z.mat))
      #calculate the completeness cutoff for each snp to be retained
      filt<-i
      df.x<-rbind(df.x, as.data.frame(cbind(filt, missingness, snps.retained)))
      }
    
    #make columns numeric for plotting
    df.x$filt<-as.numeric(as.character(df.x$filt))
    df.x$missingness<-as.numeric(as.character(df.x$missingness))
    df.x$snps.retained<-as.numeric(as.character(df.x$snps.retained))
    
    #visualize dotplot for total loci retained at each filter level
    plot1<-ggplot2::ggplot(df.x, ggplot2::aes(x=filt)) +
      ggplot2::scale_x_continuous(breaks=c(.3,.5,.6,.65,.7,.75,.8,.85,.9,.95,1))+
      ggplot2::geom_point(ggplot2::aes(y=snps.retained)) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "SNP completeness cutoff", y = "total SNPs retained")

    #visualize dotplot for missing data % at each filter level
    plot2<-ggplot2::ggplot(df.x, ggplot2::aes(x=filt)) +
      ggplot2::scale_x_continuous(breaks=c(.3,.5,.6,.65,.7,.75,.8,.85,.9,.95,1))+
      ggplot2::geom_point(ggplot2::aes(y=missingness)) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "SNP completeness cutoff", y = "total % missing data")

    gridExtra::grid.arrange(plot1,plot2)
  }
  
}


