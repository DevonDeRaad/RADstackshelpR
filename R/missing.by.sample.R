#leans heavily on the tutorial from Knaus:
#https://knausb.github.io/vcfR_documentation/sequence_coverage.html

#' @export
missing.by.sample <- function(vcfR, popmap=NULL, cutoff=NULL){
        
  #extract depth from the vcf
  dp<- vcfR::extract.gt(vcfR, element='DP', as.numeric=TRUE)
  
        if (is.null(cutoff)){

            if (is.null(popmap)) {
                print("No popmap provided")
                } 
            else {
                #calculate missingness by pop here and make dotplot
                  #popmap must be a two column dataframe with 'id' and 'pop' columns
                  #id's must match the ids in the vcf file
                  #calculate missingness by individual
                  miss<-colSums(is.na(dp))/nrow(dp)
                  #calculate avg depth by individual
                  avg.depth<-colMeans(dp, na.rm = TRUE)
                  #store ordered column names
                  samples<-colnames(dp)
                  #create df
                  df.x<-data.frame(id=samples,missingness=miss,avg.depth=avg.depth, row.names = NULL)
                  #bring in popmap
                  df.f<-merge(df.x, popmap, by="id")
                  #plot missingness and depth by pop
                  plot1<-ggplot2::ggplot(df.f, ggplot2::aes(x= reorder(pop, -missingness), y=missingness)) + 
                    ggplot2::geom_violin()+
                    ggplot2::theme_classic()+
                    ggplot2::theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))+
                    ggplot2::ylab("proportion missing")+
                    ggplot2::geom_dotplot(binaxis='y', stackdir='center', dotsize=.8)
                  #dotplot avg depth of sequencing       
                  plot2<-ggplot2::ggplot(df.f, ggplot2::aes(x= reorder(pop, -missingness), y=avg.depth)) + 
                    ggplot2::geom_violin()+
                    ggplot2::theme_classic()+
                    ggplot2::theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))+
                    ggplot2::ylab("average depth")+
                    ggplot2::geom_dotplot(binaxis='y', stackdir='center', dotsize=.8)
                  
                  #
                  gridExtra::grid.arrange(plot1,plot2)
                  
                }
      

    #initialize dataframe
    df.y<- data.frame(indiv=character(), snps.retained=numeric(), filt=numeric())
    #loop
    for (i in c(.5,.6,.7,.8,.9,1)){
      #get vector of individuals we are looping over
      indiv<-colnames(dp)
      #calculate the completeness cutoff for each snp to be retained in this iteration
      filt<-rep(i, times =ncol(dp))
      #calculate the number of loci successfully genotyped in each sample in this iteration
      snps.retained<-colSums(is.na(dp[(rowSums(is.na(dp))/ncol(dp) <= 1-i),]) == "FALSE")
      #append the data from this iteration to existing df
      df.y<-rbind(df.y, as.data.frame(cbind(indiv, filt, snps.retained)))
      #close for loop
    }
    
    #make columns numeric for plotting
    df.y$filt<-as.numeric(as.character(df.y$filt))
    df.y$snps.retained<-as.numeric(as.character(df.y$snps.retained))
    rownames(df.y)<-NULL
    #visualize color coded by individual
    print( 
      ggplot2::ggplot(df.y, ggplot2::aes(x=filt, y=snps.retained, color=indiv))+
        ggplot2::geom_point()+
        ggplot2::ggtitle("SNPs retained by filtering scheme") +
        ggplot2::xlab("fraction of non-missing genotypes required to retain each SNP (0-1)")+
        ggplot2::ylab("non-missing SNPs retained per sample")+
        ggplot2::theme_light()+
        ggplot2::theme(legend.position="none")
    )
    
    #calculate missingness by individual
    miss<-colSums(is.na(dp))/nrow(dp)
    #show plot with missingness by sample
    dotchart(sort(miss), cex=.5, xlab = "proportion missing data")
    abline(v=c(.5,.6,.7,.8,.9,1), lty="dashed")
    
    #return the dataframe showing the missingness and avg depth per individual
    return(df.y)
    
      }
  else {
    
    #calculate missingness by individual
    miss<-colSums(is.na(dp))/nrow(dp)
    #vis plot to show where cutoff was set
    dotchart(sort(miss), cex=.5)
    abline(v=cutoff, col="red")
    
    #print
    print(paste0(length(labels(miss)[miss > cutoff])," samples fall below ",cutoff," missing data cutoff, and were removed from VCF"))
    #drop those individuals from vcfR
    vcfR@gt <- vcfR@gt[,!(colnames(vcfR@gt) %in% labels(miss)[miss > cutoff])]

    return(vcfR)
  }
    
}

