
#' @export
hard.filter.vcf <- function(vcfR, depth=NULL, gq=NULL){
  
  #extract depth from the vcf
  dp.matrix<- vcfR::extract.gt(vcfR, element='DP', as.numeric=TRUE)
  
  #calculate the SNPs that fall below the depth filter
  i<-round((sum(dp.matrix < depth, na.rm = TRUE)/sum(!is.na(dp.matrix)))*100, 2)
  #report filter
  print(paste0(i,"% of genotypes fall below a read depth of ",depth," and were converted to NA"))
  
  #convert to NAs
  dp.matrix[dp.matrix < depth] <- NA
  vcfR@gt[,-1][ is.na(dp.matrix) == TRUE ] <- NA
  
  #extract gq from the vcf
  gq.matrix<- vcfR::extract.gt(vcfR, element='GQ', as.numeric=TRUE)
  
  #calculate the SNPs that fall below the gq filter
  j<-round((sum(gq.matrix < gq, na.rm = TRUE)/sum(!is.na(gq.matrix)))*100, 2)
  #report filter
  print(paste0(j,"% of genotypes fall below a genotype quality of ",gq," and were converted to NA"))
  
  #convert to NAs
  gq.matrix[gq.matrix < gq] <- NA
  vcfR@gt[,-1][ is.na(gq.matrix) == TRUE ] <- NA
  
  return(vcfR)
}


vcfR<-hard.filter.vcf(vcfR=x, depth = 5, gq = 30)
