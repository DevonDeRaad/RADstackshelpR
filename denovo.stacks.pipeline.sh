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

files="908108_H_diadema_Gatokae        diadema
908150_H_dinops_Guadalcanal     dinops
908151_H_diadema_Guadalcanal    diadema
908152_H_diadema_Guadalcanal    diadema
908153a_H_dinops_Guadalcanal    dinops
908154_H_dinops_Guadalcanal     dinops
908155_H_dinops_Guadalcanal     dinops
908156_H_diadema_Guadalcanal    diadema
908208_H_diadema_Guadalcanal    diadema
KVO150_H_diadema_Isabel diadema
KVO168_H_diadema_Isabel diadema
KVO169_H_diadema_Isabel diadema
KVO170_H_diadema_Isabel diadema
KVO171_H_diadema_Isabel diadema
KVO172_H_diadema_Isabel diadema
KVO242_H_dinops_Isabel  dinops
KVO243_H_dinops_Isabel  dinops
KVO244_H_dinops_Isabel  dinops
KVO245_H_dinops_Isabel  dinops
KVO246_H_dinops_Isabel  dinops
KVO248_H_dinops_Isabel  dinops
THL1156_H_demissus_Makira       demissus
THL1172_H_dinops_Guadalcanal    dinops
THL1173_H_dinops_Guadalcanal    dinops
THL17193_H_diadema_Ngella       diadema
THL17194_H_diadema_Ngella       diadema
THL17195_H_diadema_Ngella       diadema
THL17197_H_diadema_Ngella       diadema
THL17198_H_diadema_Ngella       diadema
THL17199_H_diadema_Ngella       diadema
WD1705_H_diadema_E_New_Britain  diadema
WD2047_H_diadema_Simbu_Prov     diadema
WD2074_H_diadema_Gulf_Prov      diadema"

# Build loci de novo in each sample for the single-end reads only. If paired-end reads are available, 
# they will be integrated in a later stage (tsv2bam stage).
# This loop will run ustacks on each sample, e.g.
#   ustacks -f ./samples/sample_01.1.fq.gz -o ./stacks -i 1 --name sample_01 -M 4 --gapped -p 8
#
# -M — Maximum distance (in nucleotides) allowed between stacks (default 2).
# -m — Minimum depth of coverage required to create a stack (default 3).

#optimize ustacks here: vary m3-m7, and M1-M8
#set m3
mkdir stacks_m3
id=1
for sample in $files
do
    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_m3 -i $id -m 3 -p 15
    let "id+=1"
done
## Run cstacks to compile stacks between samples.
/home/d669d153/work/stacks-2.41/cstacks -P stacks_m3 -M pipeline_popmap.txt -p 15
## Run sstacks. Match all samples supplied in the population map against the catalog.
/home/d669d153/work/stacks-2.41/sstacks -P stacks_m3 -M pipeline_popmap.txt -p 15
## Run tsv2bam to transpose the data so it is stored by locus, instead of by sample. We will include
## paired-end reads using tsv2bam. tsv2bam expects the paired read files to be in the samples
## directory and they should be named consistently with the single-end reads,
## e.g. sample_01.1.fq.gz and sample_01.2.fq.gz, which is how process_radtags will output them.
/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_m3 -M pipeline_popmap.txt -t 15
## Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided),
## align reads per sample, call variant sites in the population, genotypes in each individual.
/home/d669d153/work/stacks-2.41/gstacks -P stacks_m3 -M pipeline_popmap.txt -t 15
## Run populations completely unfiltered and filter using RADstackshelpR package
/home/d669d153/work/stacks-2.41/populations -P stacks_m3 -M pipeline_popmap.txt --vcf -t 15


#set m4
mkdir stacks_m4
id=1
for sample in $files
do
    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_m4 -i $id -m 4 -p 15
    let "id+=1"
done
/home/d669d153/work/stacks-2.41/cstacks -P stacks_m4 -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/sstacks -P stacks_m4 -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_m4 -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/gstacks -P stacks_m4 -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/populations -P stacks_m4 -M pipeline_popmap.txt --vcf -t 15


#set m5
mkdir stacks_m5
id=1
for sample in $files
do
    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_m5 -i $id -m 5 -p 15
    let "id+=1"
done
/home/d669d153/work/stacks-2.41/cstacks -P stacks_m5 -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/sstacks -P stacks_m5 -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_m5 -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/gstacks -P stacks_m5 -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/populations -P stacks_m5 -M pipeline_popmap.txt --vcf -t 15

#set m6
mkdir stacks_m6
id=1
for sample in $files
do
    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_m6 -i $id -m 6 -p 15
    let "id+=1"
done
/home/d669d153/work/stacks-2.41/cstacks -P stacks_m6 -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/sstacks -P stacks_m6 -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_m6 -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/gstacks -P stacks_m6 -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/populations -P stacks_m6 -M pipeline_popmap.txt --vcf -t 15

#set m7
mkdir stacks_m7
id=1
for sample in $files
do
    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_m7 -i $id -m 7 -p 15
    let "id+=1"
done
/home/d669d153/work/stacks-2.41/cstacks -P stacks_m7 -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/sstacks -P stacks_m7 -M pipeline_popmap.txt -p 15
/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_m7 -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/gstacks -P stacks_m7 -M pipeline_popmap.txt -t 15
/home/d669d153/work/stacks-2.41/populations -P stacks_m7 -M pipeline_popmap.txt --vcf -t 15

####set m to optimized value
#
#set M1
#mkdir stacks_M1
#id=1
#for sample in $files
#do
#    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_M1 -i $id -m x -M 1 -p 15
#    let "id+=1"
#done
#/home/d669d153/work/stacks-2.41/cstacks -P stacks_M1 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/sstacks -P stacks_M1 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_M1 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/gstacks -P stacks_M1 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/populations -P stacks_M1 -M pipeline_popmap.txt --vcf -t 15
#
##set M2
#mkdir stacks_M2
#id=1
#for sample in $files
#do
#    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_M2 -i $id -m x -M 2 -p 15
#    let "id+=1"
#done
#/home/d669d153/work/stacks-2.41/cstacks -P stacks_M2 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/sstacks -P stacks_M2 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_M2 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/gstacks -P stacks_M2 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/populations -P stacks_M2 -M pipeline_popmap.txt --vcf -t 15
#
##set M3
#mkdir stacks_M3
#id=1
#for sample in $files
#do
#    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_M3 -i $id -m x -M 3 -p 15
#    let "id+=1"
#done
#/home/d669d153/work/stacks-2.41/cstacks -P stacks_M3 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/sstacks -P stacks_M3 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_M3 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/gstacks -P stacks_M3 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/populations -P stacks_M3 -M pipeline_popmap.txt --vcf -t 15
#
##set M4
#mkdir stacks_M4
#id=1
#for sample in $files
#do
#    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_M4 -i $id -m x -M 4 -p 15
#    let "id+=1"
#done
#/home/d669d153/work/stacks-2.41/cstacks -P stacks_M4 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/sstacks -P stacks_M4 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_M4 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/gstacks -P stacks_M4 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/populations -P stacks_M4 -M pipeline_popmap.txt --vcf -t 15
#
##set M5
#mkdir stacks_M5
#id=1
#for sample in $files
#do
#    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_M5 -i $id -m x -M 5 -p 15
#    let "id+=1"
#done
#/home/d669d153/work/stacks-2.41/cstacks -P stacks_M5 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/sstacks -P stacks_M5 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_M5 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/gstacks -P stacks_M5 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/populations -P stacks_M5 -M pipeline_popmap.txt --vcf -t 15
#
##set M6
#mkdir stacks_M6
#id=1
#for sample in $files
#do
#    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_M6 -i $id -m x -M 6 -p 15
#    let "id+=1"
#done
#/home/d669d153/work/stacks-2.41/cstacks -P stacks_M6 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/sstacks -P stacks_M6 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_M6 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/gstacks -P stacks_M6 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/populations -P stacks_M6 -M pipeline_popmap.txt --vcf -t 15
#
##set M7
#mkdir stacks_M7
#id=1
#for sample in $files
#do
#    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_M7 -i $id -m x -M 7 -p 15
#    let "id+=1"
#done
#/home/d669d153/work/stacks-2.41/cstacks -P stacks_M7 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/sstacks -P stacks_M7 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_M7 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/gstacks -P stacks_M7 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/populations -P stacks_M7 -M pipeline_popmap.txt --vcf -t 15
#
##set M8
#mkdir stacks_M8
#id=1
#for sample in $files
#do
#    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_M8 -i $id -m x -M 8 -p 15
#    let "id+=1"
#done
#/home/d669d153/work/stacks-2.41/cstacks -P stacks_M8 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/sstacks -P stacks_M8 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/tsv2bam -P stacks_M8 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/gstacks -P stacks_M8 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/populations -P stacks_M8 -M pipeline_popmap.txt --vcf -t 15
#

##optimize n by trying n= M-1,M,M+1
##M-1
#mkdir stacks_nM-1
#id=1
#for sample in $files
#do
#    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_nM-1 -i $id -m x -M x -p 15
#    let "id+=1"
#done
#/home/d669d153/work/stacks-2.41/cstacks -n M-1 -P stacksnM-1/ -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/sstacks -P stacksnM-1 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/tsv2bam -P stacksnM-1 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/gstacks -P stacksnM-1 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/populations -P stacksnM-1 -M pipeline_popmap.txt --vcf -t 15
#
##M
#mkdir stacks_nM
#id=1
#for sample in $files
#do
#    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_nM -i $id -m x -M x -p 15
#    let "id+=1"
#done
#/home/d669d153/work/stacks-2.41/cstacks -n M -P stacksnM/ -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/sstacks -P stacksnM -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/tsv2bam -P stacksnM -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/gstacks -P stacksnM -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/populations -P stacksnM -M pipeline_popmap.txt --vcf -t 15
#
##M+1
#mkdir stacks_nM+1
#id=1
#for sample in $files
#do
#    /home/d669d153/work/stacks-2.41/ustacks -f ${sample}.fq.gz -o stacks_nM+1 -i $id -m x -M x -p 15
#    let "id+=1"
#done
#/home/d669d153/work/stacks-2.41/cstacks -n M+1 -P stacks_nM+1/ -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/sstacks -P stacksnM+1 -M pipeline_popmap.txt -p 15
#/home/d669d153/work/stacks-2.41/tsv2bam -P stacksnM+1 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/gstacks -P stacksnM+1 -M pipeline_popmap.txt -t 15
#/home/d669d153/work/stacks-2.41/populations -P stacksnM+1 -M pipeline_popmap.txt --vcf -t 15
#
