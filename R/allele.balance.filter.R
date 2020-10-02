
#' @export
filter.allele.balance <- function(vcfR){
  
  #extract AD from the vcf
  ad.matrix<- vcfR::extract.gt(vcfR, element='AD')
  #extract GT from the vcf
  gt.matrix<- vcfR::extract.gt(vcfR, element='GT')
  
  #mask dp matrix to include only called hets from gt matrix
  ad.matrix[gt.matrix != "0/1"]<-NA
  
  #split allele 1 depth from allele 2 depth
  al1<-structure(as.numeric(gsub(",.*", "", ad.matrix)), dim=dim(ad.matrix))
  al2<-structure(as.numeric(gsub(".*,", "", ad.matrix)), dim=dim(ad.matrix))
  
  #calculate percentage of hets failing AB filter
  AB<-al1/(al1 + al2) > .75 | al1/(al1 + al2) <.25
  p<-round(sum(AB, na.rm = TRUE) / sum(is.na(AB) == FALSE)*100, 2)
  j<-round(sum(AB, na.rm = TRUE) / sum(is.na(gt.matrix) == FALSE)*100, 2)
  
  print(paste0(p,"% of het genotypes (",j,"% of all genotypes) fall outside of .25 - .75 allele balance and were converted to NA"))
  
  #convert failing genotypes to NA
  vcfR@gt[,-1][AB]<-NA
  
  return(vcfR)
}
  
