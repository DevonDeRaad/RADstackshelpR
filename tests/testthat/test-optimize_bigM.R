context("optimize_bigM")
library(RADstackshelpR)

test_that("optimize_bigM generates output of the appropriate class (list)", {
  #find data in local directory
  #opt.m<-optimize_bigM(M1 = "~/Desktop/RADstackshelpR/inst/extdata/bigM1.vcf.gz")
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_bigM(M1 = system.file("extdata", "bigM1.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_bigM returns an object of class "list"
  expect_is(opt.m, "list" )
})


test_that("optimize_bigM generates a list with length of 5", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_bigM(M1 = system.file("extdata", "bigM1.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_bigM returns an object of class "list"
  expect_equal(length(opt.m), 4)
})


test_that("optimize_bigM generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(optimize_bigM(M1 = system.file("extdata", "denovo.stacks.pipeline.sh", package = "RADstackshelpR"))
  )
})


test_that("optimize_bigM generates a list with the appropriate names", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_bigM(M1 = system.file("extdata", "bigM1.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_bigM returns an object of class "list" with appropriately named components
  expect_equal(names(opt.m)[1], "snp")
  expect_equal(names(opt.m)[2], "loci")
  expect_equal(names(opt.m)[3], "snp.R80")
  expect_equal(names(opt.m)[4], "loci.R80")
})


test_that("optimize_bigM generates a list with each object inside being a dataframe", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_bigM(M1 = system.file("extdata", "bigM1.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_bigM returns an object of class "list", with each object inside being a "data.frame" object
  for (i in length(opt.m)){
    expect_is(opt.m[[i]], "data.frame")
  }
})


test_that("optimize_bigM generates dataframes with appropriate dimensions when all slots are filled", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_bigM(M1 = system.file("extdata", "bigM1.vcf.gz", package = "RADstackshelpR"),
                     M2 = system.file("extdata", "bigM2.vcf.gz", package = "RADstackshelpR"),
                     M3 = system.file("extdata", "bigM3.vcf.gz", package = "RADstackshelpR"),
                     M4 = system.file("extdata", "bigM4.vcf.gz", package = "RADstackshelpR"),
                     M5 = system.file("extdata", "bigM5.vcf.gz", package = "RADstackshelpR"),
                     M6 = system.file("extdata", "bigM6.vcf.gz", package = "RADstackshelpR"),
                     M7 = system.file("extdata", "bigM7.vcf.gz", package = "RADstackshelpR"),
                     M8 = system.file("extdata", "bigM8.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_bigM returns an object of class "list", with a 25 row data.frame as the first object when all slots are filled
  expect_equal(nrow(opt.m[[1]]), 152)
})





