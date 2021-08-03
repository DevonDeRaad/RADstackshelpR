
<!-- README.md is generated from README.Rmd. Please edit that file -->
RADstackshelpR <img src="man/figures/logo.png" align="right" alt="" width="120" />
==================================================================================

RADstackshelpR is designed to streamline the process of parameter optimization for denovo assembly of rad loci using the powerful computing and data visualization offered in R. You can check out the pkgdown website associated with the package, which details each function and hosts a comprehensive vignette, here: <https://devonderaad.github.io/RADstackshelpR/index.html>

Installation
------------

``` R
# Install development version from GitHub
devtools::install_github("DevonDeRaad/RADstackshelpR")
```

Usage
-----

RADstackshelpR packages a handful of useful wrapper functions together to streamline your RAD workflow.

### The pipeline

This denovo workflow relies on the program STACKS <https://catchenlab.life.illinois.edu/stacks/> which is designed to assemble loci and call variants from restriction-enzyme associated DNA sequence data (RADseq). The pipeline implemented here is based on the 2017 paper 'Lost in Parameter Space' <https://doi.org/10.1111/2041-210X.12775> which establishes clear recommendations for optimizing the parameters 'm', 'M', and 'n', during the process of assembling loci. This pipeline involves iterating over a range of values for each paramter (m: 3-7, M: 1-8, and n: M-1,M,M+1) and choosing the optimum value for each parameter. The first step would be demultiplexing your sequence data using the 'process_radtags' function from STACKS, which might look something like this, if your raw sequence file is in your working directory, and you used the enzyme 'ndeI' as your cutter:

``` Shell
/home/d669d153/work/stacks-2.3b/process_radtags -p .  -o . -b plate.1.barcodes.txt -e ndeI -r -c -q
```
More details on demultiplexing using process_radtags can be found at: <https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php>

Optimize 'ustacks' module
-----

Once you have an individual zipped fastq file for each sample, you would then run something like the following code in a terminal window, in order to iterate over the relevant values for 'm' within the 'ustacks' module: (here using 15 threads at each step to speed up computation)

``` Shell
#designate all sample ID's to a single variable called 'files', each sample should be in the directory, and the filename should match this designation except for the extension, e.g., '908108_H_diadema_Gatokae' = '908108_H_diadema_Gatokae.fq.gz'
files="908108_H_diadema_Gatokae
908150_H_dinops_Guadalcanal
908151_H_diadema_Guadalcanal
908152_H_diadema_Guadalcanal
908153a_H_dinops_Guadalcanal"

# Build loci de novo in each sample for the single-end reads only.
# -M — Maximum distance (in nucleotides) allowed between stacks (default 2).
# -m — Minimum depth of coverage required to create a stack (default 3).
#here, vary m from 3-7, and leave all other paramaters default

for i in {3..7}
do
#create a directory to hold this unique iteration:
mkdir stacks_m$i
#run ustacks with m equal to the current iteration (3-7) for each sample
id=1
for sample in $files
do
    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_m3 -i $id -m $i -p 15
    let "id+=1"
done
## Run cstacks to compile stacks between samples.
/home/d669d153/work/stacks-2.41/cstacks -P stacks_m$i -M pipeline_popmap.txt -p 15
## Run sstacks. Match all samples supplied in the population map against the catalog.
/home/d669d153/work/stacks-2.41/sstacks -P stacks_m$i -M pipeline_popmap.txt -p 15
## Run tsv2bam to transpose the data so it is stored by locus, instead of by sample.
/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_m$i -M pipeline_popmap.txt -t 15
## Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided),
## align reads per sample, call variant sites in the population, genotypes in each individual.
/home/d669d153/work/stacks-2.41/gstacks -P stacks_m$i -M pipeline_popmap.txt -t 15
## Run populations completely unfiltered and output unfiltered vcf, for input to the RADstackshelpR package
/home/d669d153/work/stacks-2.41/populations -P stacks_m$i -M pipeline_popmap.txt --vcf -t 15
done
```

Visualize the output of these 5 runs and determine the optimal value for 'm'.
-----




Then, you would repeat the process of iterating over different parameter values, this time varying the 'M' parameter from 1-8 within the 'ustacks' module, again using 15 threads, by running the following code in a terminal window:

``` Shell
# -M — Maximum distance (in nucleotides) allowed between stacks (default 2).
# -m — Minimum depth of coverage required to create a stack (default 3).
#here, vary M from 1-8, and set m to the optimized value based on prior visualizations (here 3)

for i in {1..8}
do
#create a directory to hold this unique iteration:
mkdir stacks_bigM$i
#run ustacks with M equal to the current iteration (1-8) for each sample, and m set to the optimized value (here, m=3)
id=1
for sample in $files
do
    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_m3 -i $id -m 3 -M $i -p 15
    let "id+=1"
done
/home/d669d153/work/stacks-2.41/cstacks -P stacks_bigM$i -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/sstacks -P stacks_bigM$i -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_bigM$i -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/gstacks -P stacks_bigM$i -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/populations -P stacks_bigM$i -M pipeline_popmap.txt --vcf -t 15
done
```

Visualize the output of these 8 runs and determine the optimal value for 'M'.
-----




Optimize 'cstacks' module
-----

For the last optimization iteration, you would vary the 'n' parameter within the 'cstacks' module, using the previously optimized values for 'm' and 'M', by running the following code in a terminal window:

``` Shell
# -n — Number of mismatches allowed between sample loci when build the catalog (default 1).
#here, vary 'n' across M-1, M, and M+1 (here equal to 1, 2, and 3), with 'm' set to the optimized value based on prior visualizations (here 3), and 'M' set to the optimized value (here 2).

for i in {1..3}
do
#create a directory to hold this unique iteration:
mkdir stacks_n$i
#run ustacks with n equal to the current iteration (1-3) for each sample, m = 3, and M=2
id=1
for sample in $files
do
    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_m3 -i $id -m 3 -M 2 $i -p 15
    let "id+=1"
done
/home/d669d153/work/stacks-2.41/cstacks -n $i -P stacks_n$i -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/sstacks -P stacks_n$i -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_n$i -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/gstacks -P stacks_n$i -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/populations -P stacks_n$i -M pipeline_popmap.txt --vcf -t 15
done
```

Visualize the output of these 3 runs and determine the optimal value for 'n'.
-----



You now have a denovo-assembled, parameter-optimized, unfiltered SNP dataset for your organism. Filter this dataset using the SNPfiltR package, which is designed to pick up where this workflow leaves off, or any other program of choice.

An example, comprehensive bash script for running all of these iterations in STACKS is available at <https://github.com/DevonDeRaad/RADstackshelpR/blob/master/denovo.stacks.pipeline.sh>. You can follow these steps, commenting and uncommenting the relevant parameter iterations as you go. A vignette detailing the implementation of RADstackshelpR with STACKS to produce a cohesive and reproducible pipeline for assembling RAD loci without a reference genome, lives here: <https://devonderaad.github.io/RADstackshelpR/articles/optimize_denovo.html>.

