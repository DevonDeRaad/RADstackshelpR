## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.show='hold', dpi=150------------------------------------------------
#load RADstackshelpR package
library(RADstackshelpR)
#optimize_m function will generate summary stats on your 5 iterative runs
#input can be full path to each file, or just the file name if the vcf files are in your working directory
m.out<-optimize_m(m3 = system.file("extdata", "m3.vcf.gz", package = "RADstackshelpR"),
                     m4 = system.file("extdata", "m4.vcf.gz", package = "RADstackshelpR"),
                     m5 = system.file("extdata", "m5.vcf.gz", package = "RADstackshelpR"),
                     m6 = system.file("extdata", "m6.vcf.gz", package = "RADstackshelpR"),
                     m7 = system.file("extdata", "m7.vcf.gz", package = "RADstackshelpR"))
#Assigning the output of this function to the variable 'm.out' should generate a list containing five objects of class 'data.frame' with the following characteristics: 'depth' showing depth per sample for each m value, 'snp' showing the number of non-missing SNPs retained in each sample at each m value, 'loci' showing the number of non-missing loci retained in each sample at each m value, 'snp.R80' showing the total number of SNPs retained at an 80% completeness cutoff, and 'loci.R80' showing the total number of polymorphic loci retained at an 80% completeness cutoff.
#Use this output list as input for this function, to visualize the effect of varying m on the depth of each sample
vis_depth(output = m.out)
#visualize the effect of varying m on the number of SNPs retained
vis_snps(output = m.out, stacks_param = "m")
#visualize the effect of varying m on the number of loci retained
vis_loci(output = m.out, stacks_param = "m")
#3 is the optimal m value, and will be used next to optimize M

## ---- dpi =150----------------------------------------------------------------
#optimize_bigM function will generate summary stats on your 8 iterative runs
M.out<-optimize_bigM(M1 = system.file("extdata", "bigM1.vcf.gz", package = "RADstackshelpR"),
                     M2 = system.file("extdata", "bigM2.vcf.gz", package = "RADstackshelpR"),
                     M3 = system.file("extdata", "bigM3.vcf.gz", package = "RADstackshelpR"),
                     M4 = system.file("extdata", "bigM4.vcf.gz", package = "RADstackshelpR"),
                     M5 = system.file("extdata", "bigM5.vcf.gz", package = "RADstackshelpR"),
                     M6 = system.file("extdata", "bigM6.vcf.gz", package = "RADstackshelpR"),
                     M7 = system.file("extdata", "bigM7.vcf.gz", package = "RADstackshelpR"),
                     M8 = system.file("extdata", "bigM8.vcf.gz", package = "RADstackshelpR"))
#Assigning the output of this function to the variable 'M.out' should generate a list containing four objects of class 'data.frame' with the following characteristics: 'snp' showing the number of non-missing SNPs retained in each sample at each m value, 'loci' showing the number of non-missing loci retained in each sample at each m value, 'snp.R80' showing the total number of SNPs retained at an 80% completeness cutoff, and 'loci.R80' showing the total number of polymorphic loci retained at an 80% completeness cutoff.
#use this function to visualize the effect of varying 'M' on the number of SNPs retained
vis_snps(output = M.out, stacks_param = "M")
#visualize the effect of varying 'M' on the number of polymorphic loci retained
vis_loci(output = M.out, stacks_param = "M")
#optimal value for this dataset is M = 2

## -----------------------------------------------------------------------------
#optimize n
n.out<-optimize_n(nequalsMminus1 = system.file("extdata", "nequalsmminus1.vcf.gz", package = "RADstackshelpR"),
                        nequalsM = system.file("extdata", "nequalsm.vcf.gz", package = "RADstackshelpR"),
                        nequalsMplus1 = system.file("extdata", "nequalsmplus1.vcf.gz", package = "RADstackshelpR"))
##Assigning the output of this function to the variable 'n.out' should generate a single object of class 'data.frame' showing the number of SNPs and loci retained across filtering levels for each value of n.
#visualize the effect of varying n on the number of SNPs retained
vis_snps(output = n.out, stacks_param = "n")
#visualize the effect of varying n on the number of polymorphic loci retained
vis_loci(output = n.out, stacks_param = "n")

