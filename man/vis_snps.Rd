% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_snps.R
\name{vis_snps}
\alias{vis_snps}
\title{Visualize the effect of varying a stacks parameter on the number of SNPs retained}
\usage{
vis_snps(output = NULL, stacks_param = NULL)
}
\arguments{
\item{output}{A list containing 5 dataframes generated by optimize_m()}

\item{stacks_param}{A character string indicating the stacks parameter iterated over}
}
\value{
A plot showing the number of SNPs retained at each given parameter value
}
\description{
This function takes the list of dataframes output by optimize_m(), optimize_M(), or optimize_n() as input.
The function then uses ggplot2 to visualize the effect of the given stacks on the number of SNPs retained.
}
\examples{
vis_snps(output=
readRDS(system.file("extdata","optimize.m.output.RDS",package="RADstackshelpR",mustWork=TRUE)),
         stacks_param = "m")
}
