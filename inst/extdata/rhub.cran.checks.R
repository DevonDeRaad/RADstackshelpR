
#verify that package is ready for CRAN
#install rhub
#install.packages("rhub")

#library rhub
library(rhub)

#make sure to build this tarball using the following code in a terminal window, while one level above the package repository
###/Users/devder/opt/anaconda3/bin/R CMD build RADstackshelpR
#this points to the current release of R on my local machine and builds the package tarball

#now test package build on the same operating systems that CRAN uses, by pointing this function to the built tarball for the package
cran_prep <- check_for_cran(path = "~/Desktop/RADstackshelpR_0.1.0.tar.gz")

#now run cran_summary to ensure that the package is ready for CRAN submission
cran_prep$cran_summary()
