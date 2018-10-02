# Author: Benjamin J Perry
# Copyright: Copyright 2018 CC BY-NC-SA 4.0
# Credits: Benjamin J Perry
# Version: 0.1
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Production
#
#!/bin/bash
# Overview: Nanopore + Illumina Hybrid Assembly Pipeline
# Historic: Loop for trimming used
# for i in $(ls); do porechop -i "$i"/"$i".fastq --format fastq -v 2 -t 12 -b "$i/$i"_porechop; done
# Historic: Loop for read filtering used
# for i in $(ls); do filtlong --min_length 1000 --min_window_q"$i"/"$i".adpttrim.fastq
# Historic: Command for long read evaluation
# NanoPlot -t 4 --fastq SU343.adpttrim.fastq --plots hex

# Declarations
ILLUMINAR1=""
ILLUMINAR2=""
SPADESCOROUT=""
SPADESCOROUTR1=""
SPADESCOROUTR2=""

ONTFILT=""
ILLUMINASPADESCOR=""
ONTLORDEC=""

FLYEOUT=""
FLYEGFA=""

UNICYCLEROUT=""

# Module 1
# SPAdes K-mer Correction of Illumina Data py36
spades.py -1 $ILLUMINAR1 -2 $ILLUMINAR2 -o $SPADESCOROUT --only-error-correction -t 14 -m 24

# Module 2
# LoRDEC K-mer Correction of the Nanopore Long Reads py36
lordec-correct -i $ONTFILT -2 $ILLUMINASPADESCOR -k 19 -s 3 -T 14 -p -o $ONTLORDEC

# Module 3
# Flye Genome Assembly of the LoRDEC Corrected Long Reads py27
flye --nano-corr $ONTLORDEC -g 7m -o $FLYEOUT -t 14

# Module 4
# Unicycler Hybrid Assembly with Corrected Illumina and Nanopore Reads py36
unicycler -1 $SPADESCOROUTR1 -2 $SPADESCOROUTR2 -l $ONTLORDEC --verbosity 2 --vcf -t 14 --existing_long_read_assembly $FLYEGFA -o $UNICYCLEROUT

exit
