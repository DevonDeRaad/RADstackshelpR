context("vis_snps")
library(RADstackshelpR)

test_that("vis_snps generates output of the appropriate class (ggplot)", {
  #find data in local directory
  #opt.m<-optimize_bigM(M1 = "~/Desktop/RADstackshelpR/inst/extdata/m3.vcf.gz")
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_m(m3 = system.file("extdata", "m3.vcf.gz", package = "RADstackshelpR"),
                     m4 = system.file("extdata", "m4.vcf.gz", package = "RADstackshelpR"),
                     m5 = system.file("extdata", "m5.vcf.gz", package = "RADstackshelpR"))
  #test that vis_snps returns an object of class "ggplot"
  expect_is(vis_snps(output = opt.m), "ggplot")
})


test_that("vis_snps generates a ggplot object with length of 9", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_m(m3 = system.file("extdata", "m3.vcf.gz", package = "RADstackshelpR"),
                     m4 = system.file("extdata", "m4.vcf.gz", package = "RADstackshelpR"),
                     m5 = system.file("extdata", "m5.vcf.gz", package = "RADstackshelpR"))
  #test that vis_snps returns an object with length = 9
  expect_equal(length(vis_snps(output = opt.m)), 9)
})


test_that("vis_snps generates an error if run with a non-vcf file", {
  #generate random vector
  x<-rnorm(100)
  #expect error trying to read this vector when you need the output list from optimize_m()
  expect_error(vis_snps(output = x))
})


test_that("vis_snps generates ggplot objects with appropriate dimensions when all slots are filled", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_m(m3 = system.file("extdata", "m3.vcf.gz", package = "RADstackshelpR"),
                     m4 = system.file("extdata", "m4.vcf.gz", package = "RADstackshelpR"),
                     m5 = system.file("extdata", "m5.vcf.gz", package = "RADstackshelpR"),
                     m6 = system.file("extdata", "m6.vcf.gz", package = "RADstackshelpR"),
                     m7 = system.file("extdata", "m7.vcf.gz", package = "RADstackshelpR"))
  #test that vis_snps returns an object of class "ggplot", with a 25 row data.frame as the first object when all slots are filled
  expect_equal(nrow(vis_snps(output = opt.m)[[1]]), 95)
})



