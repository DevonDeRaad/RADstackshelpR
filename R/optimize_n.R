
#open function:
optimize_n <- function(nequalsMminus1=NULL,nequalsM=NULL,nequalsMplus1=NULL){
  #initialize empty n.df
  n.df<- data.frame(m=character(), filt=as.numeric(), snps=numeric(), V4=character())
  #set vector of n identifiers
  ns<-c("nequalsMminus1","nequalsM","nequalsMplus1")
  #start on first position in vector of m identifiers
  j=1
  
  #open for loop for each m identifier
  for(x in list(nequalsMminus1,nequalsM,nequalsMplus1)){
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
        poly.loci[k]<-length(unique(stringr::str_extract(vcf.r@fix[,3][(rowSums(is.na(vcf.r@gt))/ncol(vcf.r@gt) <= 1-i)], pattern = "[0-9]+")))
        k=k+1
        #close for loop
      }
      ##cbind these three vectors with m identifier and append it all to the m df
      snpsubset<-as.data.frame(cbind(rep(ns[j], times=11), filt, snps, rep("snp", times = 11)))
      locisubset<-as.data.frame(cbind(rep(ns[j], times=11), filt, poly.loci, rep("loci", times = 11)))
      #match colnames so you can rbind these together in tidy format
      colnames(locisubset)[3]<-"snps"
      #append to existing df
      n.df<- rbind(n.df, as.data.frame(rbind(snpsubset,locisubset)))
      #close if else statement
    }
    #set j for the next m identifier for next time we go through this loop
    j=j+1
    #close for loop
  }
  
  print("Optimal n value returns the most polymorphic loci in the 80% complete matrix (Paris et al. 2017)")
  #take m df output from all these possibilities
  #plot number of SNPs retained colored by m at each filt level, as open circles
  #plot the number of polymorphic loci retained colored by m at each filt level, as closed circles
  #plot a vertical line at x=.8
  #return a message as output telling the user to pick the m value with the most polyloci retained at r80 (.8)
  #rename columns
  colnames(n.df)<-c("n","filt","retained","snp.locus")
  n.df$M<-as.character(n.df$n)
  n.df$filt<-as.numeric(as.character(n.df$filt))
  n.df$retained<-as.numeric(as.character(n.df$retained))
  n.df$snp.locus<-as.character(n.df$snp.locus)
  print(
    ggplot2::ggplot(n.df, ggplot2::aes(x=filt, y=retained, col = n, shape=snp.locus))+
      ggplot2::geom_point(alpha = .75, size=3)+
      ggplot2::ggtitle("total SNPs and polymorphic loci retained by filtering scheme")+
      ggplot2::xlab("fraction of non-missing genotypes required to retain each SNP (0-1)")+
      ggplot2::ylab("# SNPs/loci")+
      ggplot2::theme_light()+
      ggplot2::geom_vline(xintercept=.8)+
      ggplot2::labs(col = c("mismatches allowed\nbetween stacks\nduring catalogue building"), shape="")
  )
  
  #return the depth and snp/loci dataframes in case you want to do your own visualizations
  return(n.df)
  #close function
}

