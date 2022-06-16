#!/bin/sh
#
#SBATCH --job-name=denovo.rad.hippo              # Job Name
#SBATCH --nodes=1             # nodes
#SBATCH --cpus-per-task=15               # CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/scratch/diadema.dinops	# Set working d$
#SBATCH --mem-per-cpu=1gb            # memory requested
#SBATCH --time=10000

src=/home/d669d153/scratch/diadema.dinops #directory that contains a folder with sample data in fasta format

files="908108_H_diadema_Gatokae
908150_H_dinops_Guadalcanal
908151_H_diadema_Guadalcanal
908152_H_diadema_Guadalcanal
908153a_H_dinops_Guadalcanal
908154_H_dinops_Guadalcanal
908155_H_dinops_Guadalcanal
908156_H_diadema_Guadalcanal
908208_H_diadema_Guadalcanal
KVO150_H_diadema_Isabel
KVO168_H_diadema_Isabel
KVO169_H_diadema_Isabel
KVO170_H_diadema_Isabel
KVO171_H_diadema_Isabel
KVO172_H_diadema_Isabel
KVO242_H_dinops_Isabel
KVO243_H_dinops_Isabel
KVO244_H_dinops_Isabel
KVO245_H_dinops_Isabel
KVO246_H_dinops_Isabel
KVO248_H_dinops_Isabel
THL1156_H_demissus_Makira
THL1172_H_dinops_Guadalcanal
THL1173_H_dinops_Guadalcanal
THL17193_H_diadema_Ngella
THL17194_H_diadema_Ngella
THL17195_H_diadema_Ngella
THL17197_H_diadema_Ngella
THL17198_H_diadema_Ngella
THL17199_H_diadema_Ngella
WD1705_H_diadema_E_New_Britain
WD2047_H_diadema_Simbu_Prov
WD2074_H_diadema_Gulf_Prov"

# Build loci de novo in each sample for the single-end reads only.
# -M — Maximum distance (in nucleotides) allowed between stacks (default 2).
# -m — Minimum depth of coverage required to create a stack (default 3).
#here, we will vary m from 3-7, and leave all other paramaters default

for i in {3..7}
do
#create a directory to hold this unique iteration:
mkdir stacks_m$i
#run ustacks with m equal to the current iteration (3-7) for each sample
id=1
for sample in $files
do
    /home/path/to/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_m$i -i $id -m $i -p 15
    let "id+=1"
done
## Run cstacks to compile stacks between samples. Popmap is a file in working directory called 'pipeline_popmap.txt'
/home/path/to/stacks-2.41/cstacks -P stacks_m$i -M pipeline_popmap.txt -p 15
## Run sstacks. Match all samples supplied in the population map against the catalog.
/home/path/to/stacks-2.41/sstacks -P stacks_m$i -M pipeline_popmap.txt -p 15
## Run tsv2bam to transpose the data so it is stored by locus, instead of by sample.
/home/path/to/stacks-2.41/tsv2bam -P stacks_m$i -M pipeline_popmap.txt -t 15
## Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided),
## align reads per sample, call variant sites in the population, genotypes in each individual.
/home/path/to/stacks-2.41/gstacks -P stacks_m$i -M pipeline_popmap.txt -t 15
## Run populations completely unfiltered and output unfiltered vcf, for input to the RADstackshelpR package
/home/path/to/stacks-2.41/populations -P stacks_m$i -M pipeline_popmap.txt --vcf -t 15
done

####set m to optimized value determined using RADstackshelpR

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
    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_bigM$i -i $id -m 3 -M $i -p 15
    let "id+=1"
done
/home/path/to/stacks-2.41/cstacks -P stacks_bigM$i -M pipeline_popmap.txt -p 15
/home/path/to/stacks-2.41/sstacks -P stacks_bigM$i -M pipeline_popmap.txt -p 15
/home/path/to/stacks-2.41/tsv2bam -P stacks_bigM$i -M pipeline_popmap.txt -t 15
/home/path/to/stacks-2.41/gstacks -P stacks_bigM$i -M pipeline_popmap.txt -t 15
/home/path/to/stacks-2.41/populations -P stacks_bigM$i -M pipeline_popmap.txt --vcf -t 15
done

##set M to optimized value
##optimize n by trying n= M-1,M,M+1

# -n — Number of mismatches allowed between sample loci when build the catalog (default 1).
#here, vary 'n' across M-1, M, and M+1 (because my optimized 'M' value = 2, I will iterate over 1, 2, and 3 here), with 'm' and 'M' set to the optimized value based on prior visualizations (here 'm' = 3, and 'M'=2).

for i in {1..3}
do
#create a directory to hold this unique iteration:
mkdir stacks_n$i
#run ustacks with n equal to the current iteration (1-3) for each sample, m = 3, and M=2
id=1
for sample in $files
do
    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_n$i -i $id -m 3 -M 2 $i -p 15
    let "id+=1"
done
/home/path/to/stacks-2.41/cstacks -n $i -P stacks_n$i -M pipeline_popmap.txt -p 15
/home/path/to/stacks-2.41/sstacks -P stacks_n$i -M pipeline_popmap.txt -p 15
/home/path/to/stacks-2.41/tsv2bam -P stacks_n$i -M pipeline_popmap.txt -t 15
/home/path/to/stacks-2.41/gstacks -P stacks_n$i -M pipeline_popmap.txt -t 15
/home/path/to/stacks-2.41/populations -P stacks_n$i -M pipeline_popmap.txt --vcf -t 15
done
