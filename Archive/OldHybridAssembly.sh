#!/bin/bash
#
# 2018 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
#
# Author Benjamin .J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Revised: 06/05/2018

# This shell script takes one fastq format ONT long reads file and two Illumina
# paired read files as arguments.
#
# ./ont_assembly_polish.sh ont.reads.fastq illumina.R1.fastq.gz illumina.R2.fastq.gz
#
# The code will return an assembled and polished assembly using multiple steps:
# 1. Reads are assembled de novo using Flye deBruijn graph mediated assembler.
#   - Flye is also used to poish the genome 3 times post assembly.
# 2. The preliminary assembly has the contig start position modified to dnaA
#   using the circularization algorithm cirlator.
# 3. The updated assembly is then polished 3 times with Racon
# 4. The polished genome is the again ploished with Pilon, this time encorperating
#   paired end information in the polishing.
# 5. The final Assemly is then copied to the working directory for the assembly.
#
# Directory Sturcuture:
# /GENOME
#   ./GENOME.ont.reads.fastq
#   ./illumina.R1.fastq.gz
#   ./illumina.R2.fastq.gz
# Note: the "GENOME" substring at the root of the ont reads will be used for
# subsiquent nameing. Illumina naming is not used.
#
# Copmuting Environment:
# All the dependencies needed for this pipline to funtction properly can be installed
# using bioconda. The following packages, and their dependencies, are required:
# 1. flye (requires python 2.7 conda env built with name "py27")
# 2. circlator
# 3. racon
# 4. pilon
# 5. zcat
# 6. gzip
# 7. bowtie2

clear
printf "\n"
printf "#####################################################################\n"
printf "###                 ONP Assembly and Polish Pipeline              ###\n"
printf "#####################################################################\n"
printf "\n\n"
printf "2018 Benjamin J Perry - (CC BY-NC-SA 4.0)\n"
printf "Author: Benjamin .J Perry\n"
printf "Email: benjamin.perry@postgrad.otago.ac.nz\n"
printf "Revised: 23/04/2018\n\n"
sleep 5

printf "Begin Execution at: $(date)\n\n"
START=`date +%s`
sleep 1

###  Check for read files ###
if [ -f $1 ]; then
	ONPREADS=$1
	printf "Oxford Nanopore Read File: $ONPREADS\n\n"
  GENOMEBASE=$(echo $ONPREADS | cut -d "." -f1)
else
	printf "ERRORL: No file designated for ONP reads.\n"
	printf "Usage: ./onp_assembly_polish.sh onpreads.fastq illuminareads.R1.fastq.gz illuminareads.R2.fastq.gz\n\n"
	exit 1
fi

if [ -f $2 ] && [ -f $3 ]; then
	READ1=$2
	READ2=$3
	printf "Illumina Read Files: $READ1 $READ2\n\n"
else
	printf "ERROR: Illumina read files are not properly indicated.\n"
	printf "Usage: ./onp_assembly_polish.sh onpreads.fastq illuminareads.R1.fastq.gz illuminareads.R2.fastq.gz\n"
	exit 1
fi

### Updated read files ###
ILLREAD1="$(echo $READ1 | cut -d "." -f1)".R1.updated.fastq.gz
ILLREAD2="$(echo $READ2 | cut -d "." -f1)".R2.updated.fastq.gz
# Update read 1 and read 2 headerlined to make them unique,
# substitute {9}0 substring with R1 or R2 respectively.
printf "Updating read 1 header line...\n"
zcat $READ1 | sed 's/\s/R1 /' | gzip > "$ILLREAD1"
printf "Created file: $ILLREAD1\n\n"

printf "Updating read 2 header line...\n"
zcat $READ2 | sed 's/\s/R2 /' | gzip > "$ILLREAD2"
printf "Created file: $ILLREAD2\n\n"

Make reads sub directory
Move reads into
mkdir reads
mv $ONPREADS reads
mv $ILLREAD1 reads
mv $ILLREAD2 reads
mv $READ1 reads
mv $READ2 reads
# Update variables
ILLREAD1=$(readlink -f reads/$ILLREAD1)
printf "Illumina read file 1: $ILLREAD1\n"
ILLREAD2=$(readlink -f reads/$ILLREAD2)
printf "Illumina read file 2: $ILLREAD2\n\n"
ONPREADS=reads/$ONPREADS
printf "ONP read file: $ONPREADS\n\n"
sleep 2

# Preparing read files for use in Racon
printf "Concatenating Illumina read files for Racon...\n"
CATREADS=reads/temp.fastq.gz
zcat $ILLREAD1 $ILLREAD2 | gzip -c > $CATREADS
CATREADS=$(readlink -f $CATREADS)
printf "Completed concatenation.\n\n"
sleep 2

### Flye Genome Assembly ###
printf "\n\n\n"
printf "#####################################################################\n"
printf "###                   1.  Flye Genome Assembly                    ###\n"
printf "#####################################################################\n"
printf "\n"

# Takes the root of the file. ie. the substring prior to the first '.' character
GENOMEBASEDIR="$GENOMEBASE"_flye
ITERGENOME="$GENOMEBASEDIR"/contigs.fasta

printf "genome base: $GENOMEBASE\n"
printf "genome base dir: $GENOMEBASEDIR\n\n"

source activate py27

flye --nano-raw $ONPREADS -g 6.5m -o $GENOMEBASEDIR -t 16 -i 3 --debug

source deactivate

printf "Completed ONP read Assembly using Flye.\n"
printf "Initial assembly: $ITERGENOME\n\n"

# printf "#####################################################################\n"
# printf "###              2.  Correct Assembly Linearization               ###\n"
# printf "#####################################################################\n"
# printf "\n"
# printf "Experimental: Adjusting assembly start position with circlator.\n\n"
#
# CIRCLATORDIR="$GENOMEBASE"_Circlator
# mkdir $CIRCLATORDIR
# cd $CIRCLATORDIR
# circlator fixstart --verbose ../"$ITERGENOME" $GENOMEBASE
# cd ..
#
# ITERGENOME="$CIRCLATORDIR"/"$GENOMEBASE".fasta
#
printf "\n\nPolishing: $ITERGENOME\n"


printf "#####################################################################\n"
printf "###                   2.  Racon Genome Polish                     ###\n"
printf "#####################################################################\n"
printf "\n"

# racon [options ...] <sequences> <overlaps> <target sequences>
for I in {1..3};
do
	printf "Racon Polish $I of 3\n"
	printf "#####################################################################\n"

	#	Build bowtie2 index for current iteration
	RACONBASEDIR=racoon0"$I"
	mkdir $RACONBASEDIR
	INDX="$RACONBASEDIR"/"$GENOMEBASE""$I"
	printf "Building bowtie2 index: $INDX\n"
	printf "Buidling...\n"
	sleep 1
	bowtie2-build $ITERGENOME $INDX
	printf "Completed buidling bowtie2 index.\n\n"

	ALIGNDREADS="$RACONBASEDIR"/"$GENOMEBASE"."$I".sam
	printf "Aligning $CATREADS to $INDX with bowtie2...\n\n"

	bowtie2 -p 16 --non-deterministic --very-sensitive -x $INDX -U $CATREADS -S $ALIGNDREADS
	printf	"\nCompleted bowtie2 alignment.\n\n"

	### OPTIONAL: You can align with BWA instead of bowtie2
	# printf "Building bwa index: $INDX\n"
	# printf "Buidling...\n"
	# sleep 1
	# bwa index -p $INDX $ITERGENOME
	# printf "Completed buidling bwa index.\n\n"
	#
	# printf "Concatenating Illumina reads for alignment...\n"
	# CATREADS="$RACONBASEDIR"/temp.fastq.gz
	# zcat $ILLREAD1 $ILLREAD2 | gzip -c > $CATREADS
	# printf "Completed concatenation.\n\n"
	#
	# ALIGNDREADS="$RACONBASEDIR"/"$GENOMEBASE"."$I".sam
	#
	# printf "Aligning $CATREADS to $INDX with bwa...\n"
	# bwa mem -t 14 -o $ALIGNDREADS $INDX $CATREADS
	# printf	"Completed bwa alignment.\n"
	# printf "Aligned reads ouput file: $ALIGNDREADS\n\n"

	#	Pass genome.iteration, illumina.reads, and iteration.alignment to racon
	printf "Begining racon polish iteration: $I...\n"
	POLISHEDGENOME="$RACONBASEDIR"/"$GENOMEBASE".racon.polish."$I".fasta
	racon -t 16 $CATREADS $ALIGNDREADS $ITERGENOME > $POLISHEDGENOME

	printf "Completed $GENOMEBASE racon polish iteration $I of 3.\n"
	printf	"Polished genome saved to: $POLISHEDGENOME\n\n"

	#	Update variables
	ITERGENOME=$(readlink -f $POLISHEDGENOME)
	rm "$RACONBASEDIR"/*.sam
	rm "$RACONBASEDIR"/*.bt2
done
sleep 2

### Pilon Genome Polish ###
printf "\n\n\n"
printf "#####################################################################\n"
printf "###                   5. Pilon Genome Polish                      ###\n"
printf "#####################################################################\n"
printf "\n"

for I in {1..3};
do
	printf "Pilon Polish $I of 3\n"
	printf "#####################################################################\n"
	# pilon --genome genome.fasta [--frags frags.bam] [--jumps jumps.bam] [--unpaired unpaired.bam] [...other options...]
	#	Build bowtie2 index for current iteration
	PILONBASEDIR=pilon"$I"
	mkdir $PILONBASEDIR
	INDX="$PILONBASEDIR"/"$GENOMEBASE""$I"
	printf "Building bowtie2 index: $INDX\n"
	printf "Buidling...\n"
	sleep 1
	bowtie2-build $ITERGENOME $INDX
	printf "Completed buidling bowtie2 index.\n\n"
	printf "Aligning paired reads "$(echo $READ1 | cut -d "." -f1)" to $INDX with bowtie2...\n"

	ALIGNDREADS="$PILONBASEDIR"/"$GENOMEBASE"."$I".bam
	bowtie2 -p 16 --non-deterministic --very-sensitive -x $INDX -1 $ILLREAD1 -2 $ILLREAD2 | samtools view -bS | samtools sort > $ALIGNDREADS
	samtools index	$ALIGNDREADS
	printf	"Completed bowtie2 alignment.\n\n"

	### OPTIONAL: You can align with BWA instead of bowtie2
	# printf "Building bwa index: $INDX\n"
	# printf "Buidling...\n"
	# sleep 1
	# bwa index -p $INDX $ITERGENOME
	# printf "Completed buidling bwa index.\n\n"
	#
	# printf "Concatenating Illumina reads for alignment...\n"
	# CATREADS="$RACONBASEDIR"/temp.fastq.gz
	# zcat $ILLREAD1 $ILLREAD2 | gzip -c > $CATREADS
	# printf "Completed concatenation.\n\n"
	#
	# ALIGNDREADS="$RACONBASEDIR"/"$GENOMEBASE"."$I".sam
	#
	# printf "Aligning $CATREADS to $INDX with bwa...\n"
	# bwa mem -t 14 -o $ALIGNDREADS $INDX $CATREADS
	# printf	"Completed bwa alignment.\n"
	# printf "Aligned reads ouput file: $ALIGNDREADS\n\n"

	#	Pass genome.iteration, illumina.reads, and iteration.alignment to racon
	printf "Begining Pilon polish iteration: $I...\n"
	POLISHEDGENOME="$GENOMEBASE".pilon.polish."$I"
	pilon --genome $ITERGENOME --frags $ALIGNDREADS --output $POLISHEDGENOME --outdir $PILONBASEDIR --threads 16 --verbose --changes --tracks

	printf "Correct assembly saved at: "$POLISHEDGENOME"/"$PILONBASEDIR"\n"
	printf "Completed $GENOMEBASE pilon polish iteration $I of 3.\n\n"

	#	Update variables
	ITERGENOME=$(readlink -f "$PILONBASEDIR"/"$POLISHEDGENOME".fasta)
	# printf "Refernce assembly for polish updated to: $ITERGENOME\n\n\n"
done
sleep 2

printf "#####################################################################\n"
printf "###              3. Correct Assembly Linearization                ###\n"
printf "#####################################################################\n"
printf "\n"
printf "Experimental: Adjusting assembly start position with circlator.\n\n"

CIRCLATORDIR="$GENOMEBASE"_Circlator
mkdir $CIRCLATORDIR
cd $CIRCLATORDIR
circlator fixstart --genes_fa /SCRATCH/benScratch/Meso_dnaA_repA.fasta --verbose "$ITERGENOME" $GENOMEBASE
cd ..

ITERGENOME="$CIRCLATORDIR"/"$GENOMEBASE".fasta
ITERGENOME=$(readlink -f $ITERGENOME)

printf "#####################################################################\n"
printf "###               4. Racon Genome Polish Round 2                  ###\n"
printf "#####################################################################\n"
printf "\n"

# racon [options ...] <sequences> <overlaps> <target sequences>
for I in {1..2};
do
	printf "Racon Polish 2 $I of 3\n"
	printf "#####################################################################\n"

	#	Build bowtie2 index for current iteration
	RACONBASEDIR=racoon20"$I"
	mkdir $RACONBASEDIR
	INDX="$RACONBASEDIR"/"$GENOMEBASE""$I"
	printf "Building bowtie2 index: $INDX\n"
	printf "Buidling...\n"
	sleep 1
	bowtie2-build $ITERGENOME $INDX
	printf "Completed buidling bowtie2 index.\n\n"

	ALIGNDREADS="$RACONBASEDIR"/"$GENOMEBASE"."$I".sam
	printf "Aligning $CATREADS to $INDX with bowtie2...\n\n"

	bowtie2 -p 16 --non-deterministic --very-sensitive -x $INDX -U $CATREADS -S $ALIGNDREADS
	printf	"\nCompleted bowtie2 alignment.\n\n"

	### OPTIONAL: You can align with BWA instead of bowtie2
	# printf "Building bwa index: $INDX\n"
	# printf "Buidling...\n"
	# sleep 1
	# bwa index -p $INDX $ITERGENOME
	# printf "Completed buidling bwa index.\n\n"
	#
	# printf "Concatenating Illumina reads for alignment...\n"
	# CATREADS="$RACONBASEDIR"/temp.fastq.gz
	# zcat $ILLREAD1 $ILLREAD2 | gzip -c > $CATREADS
	# printf "Completed concatenation.\n\n"
	#
	# ALIGNDREADS="$RACONBASEDIR"/"$GENOMEBASE"."$I".sam
	#
	# printf "Aligning $CATREADS to $INDX with bwa...\n"
	# bwa mem -t 14 -o $ALIGNDREADS $INDX $CATREADS
	# printf	"Completed bwa alignment.\n"
	# printf "Aligned reads ouput file: $ALIGNDREADS\n\n"

	#	Pass genome.iteration, illumina.reads, and iteration.alignment to racon
	printf "Begining racon polish iteration: $I...\n"
	POLISHEDGENOME="$RACONBASEDIR"/"$GENOMEBASE".racon.polish."$I".fasta
	racon -t 16 $CATREADS $ALIGNDREADS $ITERGENOME > $POLISHEDGENOME

	printf "Completed $GENOMEBASE racon polish iteration $I of 3.\n"
	printf	"Polished genome saved to: $POLISHEDGENOME\n\n"

	#	Update variables
	ITERGENOME=$(readlink -f $POLISHEDGENOME)

	rm "$RACONBASEDIR"/*.sam
        rm "$RACONBASEDIR"/*.bt2

done
sleep 2

printf "\n\n\n"
printf "#####################################################################\n"
printf "###                   5. Pilon Genome Polish                      ###\n"
printf "#####################################################################\n"
printf "\n"

for I in {1..2};
do
	printf "Pilon Polish $I of 3\n"
	printf "#####################################################################\n"
	# pilon --genome genome.fasta [--frags frags.bam] [--jumps jumps.bam] [--unpaired unpaired.bam] [...other options...]
	#	Build bowtie2 index for current iteration
	PILONBASEDIR=pilon20"$I"
	mkdir $PILONBASEDIR
	INDX="$PILONBASEDIR"/"$GENOMEBASE""$I"
	printf "Building bowtie2 index: $INDX\n"
	printf "Buidling...\n"
	sleep 1
	bowtie2-build $ITERGENOME $INDX
	printf "Completed buidling bowtie2 index.\n\n"
	printf "Aligning paired reads "$(echo $READ1 | cut -d "." -f1)" to $INDX with bowtie2...\n"

	ALIGNDREADS="$PILONBASEDIR"/"$GENOMEBASE"."$I".bam
	bowtie2 -p 16 --non-deterministic --very-sensitive -x $INDX -1 $ILLREAD1 -2 $ILLREAD2 | samtools view -bS | samtools sort > $ALIGNDREADS
	samtools index	$ALIGNDREADS
	printf	"Completed bowtie2 alignment.\n\n"

	### OPTIONAL: You can align with BWA instead of bowtie2
	# printf "Building bwa index: $INDX\n"
	# printf "Buidling...\n"
	# sleep 1
	# bwa index -p $INDX $ITERGENOME
	# printf "Completed buidling bwa index.\n\n"
	#
	# printf "Concatenating Illumina reads for alignment...\n"
	# CATREADS="$RACONBASEDIR"/temp.fastq.gz
	# zcat $ILLREAD1 $ILLREAD2 | gzip -c > $CATREADS
	# printf "Completed concatenation.\n\n"
	#
	# ALIGNDREADS="$RACONBASEDIR"/"$GENOMEBASE"."$I".sam
	#
	# printf "Aligning $CATREADS to $INDX with bwa...\n"
	# bwa mem -t 14 -o $ALIGNDREADS $INDX $CATREADS
	# printf	"Completed bwa alignment.\n"
	# printf "Aligned reads ouput file: $ALIGNDREADS\n\n"

	#	Pass genome.iteration, illumina.reads, and iteration.alignment to racon
	printf "Begining Pilon polish iteration: $I...\n"
	POLISHEDGENOME="$GENOMEBASE".pilon.polish."$I"
	pilon --genome $ITERGENOME --frags $ALIGNDREADS --output $POLISHEDGENOME --outdir $PILONBASEDIR --threads 16 --verbose --changes --tracks

	printf "Correct assembly saved at: "$POLISHEDGENOME"/"$PILONBASEDIR"\n"
	printf "Completed $GENOMEBASE pilon polish iteration $I of 3.\n\n"

	#	Update variables
	ITERGENOME=$(readlink -f "$PILONBASEDIR"/"$POLISHEDGENOME".fasta)
	# printf "Refernce assembly for polish updated to: $ITERGENOME\n\n\n"
done
sleep 2

printf "#######################################################################\n"
printf "#######################################################################\n"
printf "\n\n\nCompleted ONP genome assembly and polish at: $(date)\n"

GENOMEFINAL="$GENOMEBASE".polished.final.fasta
cp $ITERGENOME $GENOMEFINAL
printf "Polished genome saved at: $(readlink -f $GENOMEFINAL)\n\n"

END=`date +%s`
RUNTIME=$((END-START))
RUNTIME=$((RUNTIME/60))
printf "Total runtime: $RUNTIME m\n\n"
printf "K Bye.\n"
