
devtools::install_github("DevonDeRaad/RADstackshelpR")

#clean package
#library(pkgdown)
#clean_site(pkg = "RADstackshelpR")
#library(devtools)
#clean_vignettes(pkg = "RADstackshelpR")

#read in vcf and make subset example datasets
library(vcfR)

vcffog<-read.vcfR("~/Desktop/hipposideros/m_3.vcf")

vc.300<-vcffog[sample.int(110576, 90), c(1:20)]
vc.275<-vcffog[sample.int(110576, 80), c(1:20)]
vc.250<-vcffog[sample.int(110576, 70), c(1:20)]
vc.225<-vcffog[sample.int(110576, 60), c(1:20)]
vc.200<-vcffog[sample.int(110576, 50), c(1:20)]

#write.vcf(vc.300, file = "~/Desktop/RADstackshelpR/inst/extdata/m3.vcf.gz")
#write.vcf(vc.275, file = "~/Desktop/RADstackshelpR/inst/extdata/m4.vcf.gz")
#write.vcf(vc.250, file = "~/Desktop/RADstackshelpR/inst/extdata/m5.vcf.gz")
#write.vcf(vc.225, file = "~/Desktop/RADstackshelpR/inst/extdata/m6.vcf.gz")
#write.vcf(vc.200, file = "~/Desktop/RADstackshelpR/inst/extdata/m7.vcf.gz")

opt.m<- optimize_m(m3 = system.file("extdata", "m3.vcf.gz", package = "RADstackshelpR"),
                   m4 = system.file("extdata", "m4.vcf.gz", package = "RADstackshelpR"),
                   m5 = system.file("extdata", "m5.vcf.gz", package = "RADstackshelpR"),
                   m6 = system.file("extdata", "m6.vcf.gz", package = "RADstackshelpR"),
                   m7 = system.file("extdata", "m7.vcf.gz", package = "RADstackshelpR"))

vis_loci(opt.m)

#Now do big M
vc.1<-vcffog[sample.int(110576, 100), c(1:20)]
vc.2<-vcffog[sample.int(110576, 90), c(1:20)]
vc.3<-vcffog[sample.int(110576, 80), c(1:20)]
vc.4<-vcffog[sample.int(110576, 70), c(1:20)]
vc.5<-vcffog[sample.int(110576, 60), c(1:20)]
vc.6<-vcffog[sample.int(110576, 50), c(1:20)]
vc.7<-vcffog[sample.int(110576, 40), c(1:20)]
vc.8<-vcffog[sample.int(110576, 30), c(1:20)]

#write.vcf(vc.1, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM1.vcf.gz")
#write.vcf(vc.2, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM2.vcf.gz")
#write.vcf(vc.3, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM3.vcf.gz")
#write.vcf(vc.4, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM4.vcf.gz")
#write.vcf(vc.5, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM5.vcf.gz")
#write.vcf(vc.6, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM6.vcf.gz")
#write.vcf(vc.7, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM7.vcf.gz")
#write.vcf(vc.8, file = "~/Desktop/RADstackshelpR/inst/extdata/bigM8.vcf.gz")

opt.bigm<- optimize_bigM(M1 = system.file("extdata", "bigM1.vcf.gz", package = "RADstackshelpR"),
                   M2 = system.file("extdata", "bigM2.vcf.gz", package = "RADstackshelpR"),
                   M3 = system.file("extdata", "bigM3.vcf.gz", package = "RADstackshelpR"),
                   M4 = system.file("extdata", "bigM4.vcf.gz", package = "RADstackshelpR"),
                   M5 = system.file("extdata", "bigM5.vcf.gz", package = "RADstackshelpR"),
                   M6 = system.file("extdata", "bigM6.vcf.gz", package = "RADstackshelpR"),
                   M7 = system.file("extdata", "bigM7.vcf.gz", package = "RADstackshelpR"),
                   M8 = system.file("extdata", "bigM8.vcf.gz", package = "RADstackshelpR"))

vis_loci(opt.bigm)

#now do n
#write.vcf(vc.8, file = "~/Desktop/RADstackshelpR/inst/extdata/nequalsmminus1.vcf.gz")
#write.vcf(vc.6, file = "~/Desktop/RADstackshelpR/inst/extdata/nequalsm.vcf.gz")
#write.vcf(vc.7, file = "~/Desktop/RADstackshelpR/inst/extdata/nequalsmplus1.vcf.gz")

opt.n<- optimize_n(nequalsMminus1 = system.file("extdata", "nequalsmminus1.vcf.gz", package = "RADstackshelpR"),
                         nequalsM = system.file("extdata", "nequalsm.vcf.gz", package = "RADstackshelpR"),
                         nequalsMplus1 = system.file("extdata", "nequalsmplus1.vcf.gz", package = "RADstackshelpR"))

vis_loci(opt.n)


