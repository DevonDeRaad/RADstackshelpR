context("optimize_n")
library(RADstackshelpR)

test_that("optimize_n generates output of the appropriate class (list)", {
  #find data in local directory
  #opt.m<-optimize_n(nequalsMminus1 = "~/Desktop/RADstackshelpR/inst/extdata/nequalsmminus1.vcf.gz")
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_n(nequalsMminus1 = system.file("extdata", "nequalsmminus1.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_n returns an object of class "list"
  expect_is(opt.m, "list" )
})


test_that("optimize_n generates a list with length of 5", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_n(nequalsMminus1 = system.file("extdata", "nequalsmminus1.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_n returns an object of class "list"
  expect_equal(length(opt.m), 4)
})


test_that("optimize_n generates an error if run with a non-vcf file", {
  #expect error trying to read this vector as a vcf file
  expect_error(optimize_n(nequalsMminus1 = system.file("extdata", "denovo.stacks.pipeline.sh", package = "RADstackshelpR"))
  )
})


test_that("optimize_n generates a list with the appropriate names", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_n(nequalsMminus1 = system.file("extdata", "nequalsmminus1.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_n returns an object of class "list" with appropriately named components
  expect_equal(names(opt.m)[1], "snp")
  expect_equal(names(opt.m)[2], "loci")
  expect_equal(names(opt.m)[3], "snp.R80")
  expect_equal(names(opt.m)[4], "loci.R80")
})


test_that("optimize_n generates a list with each object inside being a dataframe", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_n(nequalsMminus1 = system.file("extdata", "nequalsmminus1.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_n returns an object of class "list", with each object inside being a "data.frame" object
  for (i in length(opt.m)){
    expect_is(opt.m[[i]], "data.frame")
  }
})


test_that("optimize_n generates dataframes with appropriate dimensions when all slots are filled", {
  #find data in package using CRAN friendly syntax
  opt.m<- optimize_n(nequalsMminus1 = system.file("extdata", "nequalsmminus1.vcf.gz", package = "RADstackshelpR"),
                        nequalsM = system.file("extdata", "nequalsm.vcf.gz", package = "RADstackshelpR"),
                        nequalsMplus1 = system.file("extdata", "nequalsmplus1.vcf.gz", package = "RADstackshelpR"))
  #test that optimize_n returns an object of class "list", with a 25 row data.frame as the first object when all slots are filled
  expect_equal(nrow(opt.m[[1]]), 57)
})





