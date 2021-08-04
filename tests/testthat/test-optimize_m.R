context("optimize_m")
library(RADstackshelpR)

test_that("optimize_m generates output of the appropriate class", {
  #find data in local directory
  opt.m<-optimize_m(m3 = "~/Desktop/RADstackshelpR/inst/extdata/populations.snps.vcf.gz")
  #find data in package using CRAN friendly syntax
  #opt.m<- optimize_m(system.file("extdata", "populations.snps.vcf", package = "RADstackshelpR"))
  #test that optimize_m returns an object of class "list"
  expect_equal( class(opt.m), "list" ) }
)


test_that("optimize_m generates output with the appropriate dimensions", {
  #find data in local directory
  opt.m<-optimize_m(m3 = "~/Desktop/RADstackshelpR/inst/extdata/populations.snps.vcf.gz")
  #find data in package using CRAN friendly syntax
  #opt.m<- optimize_m(system.file("extdata", "populations.snps.vcf", package = "RADstackshelpR"))
  #test that optimize_m returns an object of class "list"
  expect_equal( length(opt.m), 5 ) }
)
