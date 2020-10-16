
<!-- README.md is generated from README.Rmd. Please edit that file -->
RADstackshelpR <img src="man/figures/logo.png" align="right" alt="" width="120" />
==================================================================================

RADstackshelpR is designed to integrate the typically disparate steps of parameter optimization, exploration, and filtering of a RAD dataset, into a reproducible and documented workflow in R. You can see vignettes that detail the logic behind each function and direct implementation of the functions as a cohesive pipeline at <https://devonderaad.github.io/RADstackshelpR/index.html>

Installation
------------

``` r
# Install development version from GitHub
devtools::install_github("DevonDeRaad/RADstackshelpR")
```

Usage
-----

RADstackshelpR packages a handful of useful wrapper functions together to streamline your RAD workflow.

### denovo

The denovo workflow is based heavily on the 2017 paper 'Lost in Parameter Space' <https://doi.org/10.1111/2041-210X.12775> which establishes clear recommendations for optimizing the parameters 'm', 'M', and 'n'. An example bash script for running stacks denovo iteratively over m (3-7), M (1-8), and n(m-1,m,m+1) is available from <https://github.com/DevonDeRaad/RADstackshelpR/tree/master/vignettes/>. A tutorial for visualizing each of these parameter variations and optimizing them based on the optimal 'R80' value is available at <https://devonderaad.github.io/RADstackshelpR/articles/xxx>.

### reference aligned

The reference aligned workflow simply consists of swapping the parameter iteration steps for mapping your raw reads to your reference genome using bwa, and then running gstacks and outputting an unfiltered vcf file from the STACKS populations module. An example bash script for implementing these steps is available from: <https://github.com/DevonDeRaad/RADstackshelpR/tree/master/vignettes/>

### filtering

Finally, both workflows converge on a set of filtering steps which is heavily based on the Ddocent filtering tutorial (available at <https://www.ddocent.com/filtering/>). The data visualization and vcf filtering steps are implemented via RADstackshelpR functions, and a full vignette detailing the filtering workflow is available at <https://devonderaad.github.io/RADstackshelpR/articles/filtering_vignette.html>. Much of the code in this package used to manipulate vcfR files in R and implement the various filters was adapted from the tutorial on working with reduced-representation genomic data in R by JF Tabima, BJ Knaus, and NJ Grünwald, available at <https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html>. I am extremely grateful to them for making this tutorial publicly available, as working through their exercises was the genesis of this package.
