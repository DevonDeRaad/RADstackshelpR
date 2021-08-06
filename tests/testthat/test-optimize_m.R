context("optimize_m")
library(RADstackshelpR)

test_that("optimize_m generates output of the appropriate class (list)", {
  #find data in local directory
  #opt.m<-optimize_m(m3 = "~/Desktop/RADstackshelpR/inst/extdata/m3.vcf.gz")
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_m(m3 = system.file("extdata", "m3.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_m returns an object of class "list"
  expect_is(opt.m, "list" )
})


test_that("optimize_m generates a list with length of 5", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_m(m3 = system.file("extdata", "m3.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_m returns an object of class "list"
  expect_equal(length(opt.m), 5 )
})


test_that("optimize_m generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(optimize_m(m3 = system.file("extdata", "denovo.stacks.pipeline.sh", package = "RADstackshelpR"))
)
})


test_that("optimize_m generates a list with the appropriate names", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_m(m3 = system.file("extdata", "m3.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_m returns an object of class "list" with appropriately named components
  expect_equal(names(opt.m)[1], "depth")
  expect_equal(names(opt.m)[2], "snp")
  expect_equal(names(opt.m)[3], "loci")
  expect_equal(names(opt.m)[4], "snp.R80")
  expect_equal(names(opt.m)[5], "loci.R80")
})


test_that("optimize_m generates a list with each object inside being a dataframe", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_m(m3 = system.file("extdata", "m3.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_m returns an object of class "list", with each object inside being a "data.frame" object
  for (i in length(opt.m)){
    expect_is(opt.m[[i]], "data.frame")
  }
})


test_that("optimize_m generates dataframes with appropriate dimensions when all slots are filled", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_m(m3 = system.file("extdata", "m3.vcf.gz", package = "RADstackshelpR"),
                     m4 = system.file("extdata", "m4.vcf.gz", package = "RADstackshelpR"),
                     m5 = system.file("extdata", "m5.vcf.gz", package = "RADstackshelpR"),
                     m6 = system.file("extdata", "m6.vcf.gz", package = "RADstackshelpR"),
                     m7 = system.file("extdata", "m7.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_m returns an object of class "list", with a 25 row data.frame as the first object when all slots are filled
  expect_equal(nrow(opt.m[[1]]), 95)
})





