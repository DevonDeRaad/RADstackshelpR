---
title: 'RADstackshelpR: an R package for streamlining the de novo optimization of Stacks parameters'
tags:
  - R
  - RAD sequencing
  - DNA sequencing
  - bioinformatics
  - reproducibility
authors:
  - name: Devon A. DeRaad^[Corresponding author]
    orcid: 0000-0003-3105-985X
    affiliation: 1
affiliations:
 - name: Department of Ecology & Evolutionary Biology and Biodiversity Institute, University of Kansas, Lawrence, Kansas, USA
date: 16 August 2021
bibliography: paper.bib
---

# Summary

Restriction-site Associated DNA sequencing (RADseq) is a commonly used method for rapidly and cost-effectively generating sequence data from thousands of genome-wide loci. One great advantage of this approach is that these genome-wide loci can be assembled de novo (i.e., without aligning the sequences to a reference genome) making RAD sequencing amenable for non-model organisms for which a reference genome is not available. Additionally, RADseq offers reasonable file-sizes and cost of sequencing, leading to the rapid adoption of this convenient approach for answering population genetics and phylogenetics questions over the last decade.

Unfortunately, RADseq data also comes with the downside of a dizzying array of computational choices which must be made during bioinformatic processing, especially during de novo assembly. Rampant missing data, uneven sequencing amongst samples, idiosyncratic features of the genome, and stochastic factors associated with a given sequencing run, all contribute to making the optimal parameters for assembling loci and calling SNPs from a given RAD sequencing run a moving target. The software pipeline Stacks [@Catchen:2011] is designed specifically for assembling loci and calling SNPs from RADseq data, and is designed modularly, in order to give end-users control over a variety of assembly parameters throughout the pipeline. The paper 'Lost in Parameter Space' [@Paris:2018] performed detailed explorations of parameter space for the de novo assembly of RAD loci from three empirical datasets, resulting in a clearly defined set of best practices for de novo assembly of RAD loci using the Stacks pipeline.

# Statement of need

`RADstackshelpR` is an R [@R Core Team:2019]package designed to streamline and automate the following these best practices outlined by @Paris:2018, for de novo assembly of RAD loci using the Stacks pipeline. These best practices suggest executing 16 different iterations of the Stacks pipeline with varying parameters, and this process can easily become both computationally and logistically overwhelming. `RADstackshelpR` offers an explicit, extensively documented pipeline for Stacks users to follow in order to thoroughly explore parameter space and optimize the de novo assembly and SNP calling process for their empirical dataset of choice.

`RADstackshelpR` contains a series of wrapper functions which make use internally of the R package vcfR [@Knaus:2017] to read in variant call format (vcf) files which are output directly by Stacks, and return a data frame detailing the number of polymorphic loci and SNPs contained in each vcf file at an 80% completeness cutoff ('R80') [@Paris:2018]. These output data frames can then be input directly into convenient visualization functions which internally utilize the R packages ggplot2 [@Wickham:2020], ggridges [@Claus:2018], and gridExtra [@Auguie:2017] to quickly and reproducibly generate publication quality figures detailing the optimization process for your dataset.

Now, rather than relying on series' of spreadsheets and disparate script files to haphazardly explore different parameter combinations, biologists looking to tackle de novo assembly of RADseq data can follow an explicit pipeline for thoroughly and repeatably exploring parameter space according to best practices and visualizing their results for publication.

![Fig 1. Multi-panel figure visualizing the entire optimization process facilitated by `RADstackshelpR`.\label{fig:overview}](fig1.png)

# Acknowledgements
I would like to acknowledge the authors of the R package vcfR, Brian J. Knaus and Nikolas J. Grünwald, for developing and documenting extensive resources for reading and manipulating SNP data using R, specifically the exceptional GitHub site: https://knausb.github.io/vcfR_documentation/index.html which provided invaluable guidance in the development of this package.

# References
