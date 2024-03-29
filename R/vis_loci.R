#' Visualize the effect of varying a stacks parameter on the number of polymorphic loci retained
#'
#' This function takes the list of dataframes output by optimize_m(), optimize_M(), or optimize_n() as input.
#' The function then uses ggplot2 to visualize the effect of the given stacks on the number of polymorphic loci retained, reporting which value is optimal.
#'
#' @param output A list containing 5 dataframes generated by optimize_m()
#' @param stacks_param A character string indicating the stacks parameter iterated over
#' @return A plot showing the number of polymorphic loci retained at each given parameter value
#' @examples
#' vis_loci(output =
#' readRDS(system.file("extdata","optimize.m.output.RDS",package="RADstackshelpR",mustWork=TRUE)),
#'          stacks_param = "m")
#' @export
vis_loci <- function(output=NULL, stacks_param=NULL){
#bind these variable names to the function
poly.loci <- var <- poly.loci.80 <- NULL

  #isolate SNP df
  loci.df<-output$loci
  loci.df$var<-as.factor(loci.df$var)
  loci.df$poly.loci<-as.numeric(as.character(loci.df$poly.loci))

  #isolate SNP R80 df
  loci.R80.df<-output$loci.R80
  loci.R80.df$var<-as.factor(loci.R80.df$var)
  loci.R80.df$poly.loci.80<-as.numeric(as.character(loci.R80.df$poly.loci.80))

  #isolate optimized value
  opt<-loci.R80.df[loci.R80.df$poly.loci.80 == max(loci.R80.df$poly.loci.80), ]

  cat("Visualize how different values of", stacks_param, "affect number of polymorphic loci retained.
Density plot shows the distribution of the number of loci retained in each sample,
while the asterisk denotes the total number of loci retained at an 80% completeness cutoff. The optimal value is denoted by red color.", "\n")
  return(
    ggplot2::ggplot(loci.df, ggplot2::aes(x = poly.loci, y = var)) +
      ggridges::geom_density_ridges(jittered_points = FALSE, alpha = .5) +
      ggplot2::geom_point(loci.R80.df, mapping=ggplot2::aes(x=poly.loci.80, y=var), pch=8, cex=3)+
      ggplot2::geom_point(opt, mapping=ggplot2::aes(x=poly.loci.80, y=var, col="red"), pch=8, cex=3)+
      ggplot2::theme_classic() +
      ggplot2::labs(x = "polymorphic loci retained", y = paste(stacks_param, "value")) +
      ggplot2::theme(legend.position = "none")
  )
  #close function
}
