# 2018 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 0.1.0
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

clear
printf "\n"
printf "#####################################################################\n"
printf "###             ONT Illumina Hybris Assembly Pipeline             ###\n"
printf "#####################################################################\n"
printf "\n\n"
printf "2018 Benjamin J Perry - (CC BY-NC-SA 4.0)\n"
printf "Author: Benjamin .J Perry\n"
printf "Email: benjamin.perry@postgrad.otago.ac.nz\n"
printf "Revised: 3/10/2018\n\n"
sleep 5

printf "Begin Execution at: $(date)\n\n"
START=`date +%s`
sleep 1

###  Check for read files ###
if [ -f $1 ]; then
	ONTFILT=$1
	printf "ONT Read File: $ONTFILT\n\n"
    STRAIN=$(echo $ONTFILT | cut -d "." -f 1)
else
	printf "ERRORL: No file designated for ONT reads.\n"
	printf "Usage: ./HybridAssembly.sh '$STRAIN'.chop.filt.*.fastq Illumina.*.R1.fastq.gz Illumina.*.R2.fastq.gz\n\n"
	exit 1
fi

if [ -f $2 ] && [ -f $3 ]; then
	ILLUMINAR1=$2
	ILLUMINAR2=$3
	printf "Illumina Read Files: $ILLUMINAR1 $ILLUMINAR2\n\n"
else
	printf "ERROR: Illumina read files are not properly indicated.\n"
	printf "Usage: ./HybridAssembly.sh '$STRAIN'.chop.filt.*.fastq Illumina.*.R1.fastq.gz Illumina.*.R2.fastq.gz\n\n"
	exit 1
fi

# Sub-directories
SPADESCOROUT="SPAdes_Cor/" # SPAdes makes it's own dir when specified
LORDECCOROUT="LoRDEC_Cor/" # Must explicity make dir
FLYEOUT="Flye_Assembly/" # Flye makes it's own dir when specified
UNICYCLEROUT="Unicycler_Assembly/"
ORIGINFASTA=""
# Set Environment to python v3.6
source activate py36

# Module 1
# SPAdes K-mer Correction of Illumina Data py36
# spades.py -1 "$ILLUMINAR1" -2 "$ILLUMINAR2" -o "$SPADESCOROUT" --only-error-correction -t 14 -m 24
# SPADESCOROUTR1="$SPADESCOROUT"corrected/`ls "$SPADESCOROUT"corrected | grep -e "R1"`
# SPADESCOROUTR2="$SPADESCOROUT"corrected/`ls "$SPADESCOROUT"corrected | grep -e "R2"`

# Module 2
# LoRDEC K-mer Correction of the Nanopore Long Reads py36
mkdir "$LORDECCOROUT"
LORDECCOROUTFILE="$LORDECCOROUT""$STRAIN".lordec.fasta
ILLUMINASPADESCOR="$SPADESCOROUT"corrected/"$STRAIN".cat.cor.fastq.gz
# cat "$SPADESCOROUTR1" "$SPADESCOROUTR2" > "$ILLUMINASPADESCOR"
lordec-correct -i "$ONTFILT" -2 "$ILLUMINASPADESCOR" -k 19 -s 3 -T 14 -p -o "$LORDECCOROUTFILE"


# Exit py36 Environment
source deactivate
# Set Environment to python v2.7
source activate py27

# Module 3
# Flye Genome Assembly of the LoRDEC Corrected Long Reads py27
flye --nano-corr "$LORDECCOROUTFILE" -g 7m --out-dir "$FLYEOUT" -t 14

# Exit py27 Environment
source deactivate
# Set Environment to python v3.6
source activate py36

# Module 4
# Unicycler Hybrid Assembly with Corrected Illumina and Nanopore Reads py36
FLYEGFA="$FLYEOUT"2-repeat/graph_final.gfa
mkdir "$UNICYCLEROUT"

unicycler -1 "$SPADESCOROUTR1" -2 "$SPADESCOROUTR2" -l "$LORDECCOROUTFILE" --verbosity 3 --vcf -t 14 --existing_long_read_assembly "$FLYEGFA" -o "$UNICYCLEROUT"

# Exit py36 Environment
source deactivate

# Exit Pipeline
exit
