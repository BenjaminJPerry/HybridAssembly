#!/bin/bash
#Nanopore + Illumina Hybrid Assembly Pipeline

#Historic: Loop for trimming used
#for i in $(ls); do porechop -i "$i"/"$i".fastq --format fastq -v 2 -t 12 -b "$i/$i"_porechop; done
#Historic: Loop for read filtering used
#for i in $(ls); do filtlong --min_length 1000 --min_window_q"$i"/"$i".adpttrim.fastq 
#Historic: Command for long read evaluation
#NanoPlot -t 4 --fastq SU343.adpttrim.fastq --plots hex

#Module 1
#SPAdes K-mer Correction of Illumina Data

#Module 2
#LoRDEC K-mer Correction of the Nanopore Long Reads

#Module 3
#Flye Genome Assembly of the LoRDEC Corrected Long Reads

#Module 4
#Unicycler Hybrid Assembly with Corrected Illumina and Nanopore Reads
