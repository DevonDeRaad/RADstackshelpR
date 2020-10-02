
#' @export
min.mac <- function(vcfR, popmap=NULL, min.mac=NULL){
  
  if (is.null(popmap) & is.null(min.mac)){
    stop("popmap must be provided in order to compare dapc clustering to a set of a priori defined groups")
  }
  else{}
  
  if (colnames(popmap)[1] != "id" | colnames(popmap)[2] != "pop"){
    stop("popmap must be a dataframe with two columns, 'id' and 'pop'")
  }
  else{}
  
  if (is.null(min.mac)) {
    
    #convert vcfR to matrix and make numeric  
    gt.matrix<-vcfR::extract.gt(vcfR)
    gt.matrix[gt.matrix == "0/0"]<-0
    gt.matrix[gt.matrix == "0/1"]<-1
    gt.matrix[gt.matrix == "1/1"]<-2
    class(gt.matrix) <- "numeric"
    
    #calc sfs
    sfs<-rowSums(gt.matrix, na.rm = TRUE)
    #fold sfs
    for (i in 1:length(sfs)) {
      if (sfs[i] <= sum(!is.na(gt.matrix[i,]))){}
      else {
        sfs[i]<-(sum(!is.na(gt.matrix[i,]))*2 - sfs[i])
      }
    }
    
    #hist folded mac with cutoff shown
    hist(sfs, main="folded SFS", xlab = "MAC")
    
    #convert vcfR into genlight
    genlight<-vcfR::vcfR2genlight(vcfR)

    #run dapc for mac 1,2,3,4,5,10
    for (i in c(1,2,3,4,5,10)){
      
      #filter genlight by given mac
      genlight<-genlight[,sfs >= i]
      #subset sfs vector to only samples left in the vcf
      sfs<-sfs[sfs >= i]
      
      #assign samples to the number of groups present in popmap, retain all PCAs
      grp<-adegenet::find.clusters(genlight, n.pca = ncol(gt.matrix)-1, n.clust = length(levels(popmap$pop)))
      
      #check how well that assignment matched up to the provided popmap
      samps<-merge(popmap, data.frame(group=grp$grp, id=labels(grp$grp)), by='id')
      print(paste0("for ", i, " minimum MAC cutoff, compare k means clustering to popmap assignment"))
      print(table(samps$pop, samps$group))
      
      #run dapc, retain all discriminant axes, and enough PC axes to explain 75% of variance
      dapc1<-adegenet::dapc(genlight, grp$grp, n.da = length(levels(popmap$pop))-1, pca.select = "percVar", perc.pca = 75)
      
      #plot compoplot
      adegenet::compoplot(dapc1, legend=FALSE, show.lab =TRUE, cex.names=.4, main=paste0("min. MAC ",i,", total SNPs ",length(sfs)))
      
      #print
      print(paste0("DAPC with min. MAC ", i, " and ", length(sfs), " total SNPs, complete"))
    }
    
  } 
    else {
      
    #convert vcfR to matrix and make numeric  
    gt.matrix<-vcfR::extract.gt(vcfR)
    gt.matrix[gt.matrix == "0/0"]<-0
    gt.matrix[gt.matrix == "0/1"]<-1
    gt.matrix[gt.matrix == "1/1"]<-2
    class(gt.matrix) <- "numeric"
    
    #calc sfs
    sfs<-rowSums(gt.matrix, na.rm = TRUE)
    #fold sfs
    for (i in 1:length(sfs)) {
      if (sfs[i] <= sum(!is.na(gt.matrix[i,]))){}
      else {
        sfs[i]<-(sum(!is.na(gt.matrix[i,]))*2 - sfs[i])
      }
    }
    
    #hist folded mac with cutoff shown
    hist(sfs, main="folded SFS", xlab = "MAC")
    abline(v=min.mac-1, col="red")
    
    #calculate % of SNPs to be removed, and print it
    p<-round((sum(sfs < min.mac)/length(sfs))*100, 2)
    print(paste0(p, "% of SNPs fell below a minor allele count of ", min.mac, " and were removed from the VCF"))
    
    #filter vcfR
    vcfR <- vcfR[sfs >= min.mac,]
    
    return(vcfR)

    }
  
}
  
