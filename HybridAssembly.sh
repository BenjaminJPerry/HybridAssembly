# 2020 Benjamin J Perry - Attribution-NonCommercial-ShareAlike 4.0 International
# (CC BY-NC-SA 4.0)
# Version: 1.3.0
# Maintainer: Benjamin J Perry
# Email: benjamin.perry@postgrad.otago.ac.nz
# Status: Functional

clear
printf "\n"
printf "               ONT Illumina Hybrid Assembly Pipeline\n"
printf "#####################################################################\n"
printf "\n\n"
printf "2020 Benjamin J Perry - (CC BY-NC-SA 4.0)\n"
printf "Author: Benjamin .J Perry\n"
printf "Email: benjamin.perry@postgrad.otago.ac.nz\n"
printf "Version: v1.3.0\n\n"
printf "Revised: 7/03/2020\n\n"
sleep 5

printf "Begin Execution at: $(date)\n\n"
START=`date +%s`
sleep 1

###  Check for read files ###
if [ "$1" == "" ]; then
	printf "No ONT Reads Indicated...\n"
	printf "Usage: ./HybridAssembly.sh \$SAMPLE.chop.filt.*.fastq Illumina.*.R1.fastq.gz Illumina.*.R2.fastq.gz\n\n"
	exit 1
fi
if [ "$2" == "" ]; then
	printf "No Illumina R1 File Indicated...\n"
	printf "Usage: ./HybridAssembly.sh \$SAMPLE.chop.filt.*.fastq Illumina.*.R1.fastq.gz Illumina.*.R2.fastq.gz\n\n"
	exit 1
fi
if [ "$3" == "" ]; then
	printf "No Illumina R2 File Indicated...\n"
	printf "Usage: ./HybridAssembly.sh \$SAMPLE.chop.filt.*.fastq Illumina.*.R1.fastq.gz Illumina.*.R2.fastq.gz\n\n"
	exit 1
fi
if [ -f $1 ]; then
	ONTFILT=$1
	printf "ONT Read File: $ONTFILT\n\n"
    STRAIN=$(echo $ONTFILT | cut -d "." -f 1)
else
	printf "ERROR: No file designated for ONT reads.\n"
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

# Set Environment to python v3.6
source activate HybridAsBro

printf "\n"
printf "               Module 1: SPAdes k-mer Correction\n"
printf "#####################################################################\n"
printf "\n\n"

# Module 1
# SPAdes K-mer Correction of Illumina Data
spades.py -1 "$ILLUMINAR1" -2 "$ILLUMINAR2" -o "$SPADESCOROUT" --only-error-correction -t 14 -m 24
SPADESCOROUTR1="$SPADESCOROUT"corrected/`ls "$SPADESCOROUT"corrected | grep -e "R1"`
SPADESCOROUTR2="$SPADESCOROUT"corrected/`ls "$SPADESCOROUT"corrected | grep -e "R2"`

printf "\n"
printf "               Module 2: LoRDEC Long-Read k-mer Correction\n"
printf "#####################################################################\n"
printf "\n\n"

# Module 2
# LoRDEC K-mer Correction of the Nanopore Long Reads
mkdir "$LORDECCOROUT"

LORDECCOROUTFILEK19="$LORDECCOROUT""$STRAIN".lordec.k19.fasta
LORDECCOROUTFILEK31="$LORDECCOROUT""$STRAIN".lordec.k31.fasta
LORDECCOROUTFILEK41="$LORDECCOROUT""$STRAIN".lordec.k41.fasta
ILLUMINASPADESCOR="$SPADESCOROUT"corrected/"$STRAIN".cat.cor.fastq.gz
cat "$SPADESCOROUTR1" "$SPADESCOROUTR2" > "$ILLUMINASPADESCOR"

lordec-correct -i "$ONTFILT" -2 "$ILLUMINASPADESCOR" -k 19 -s 4 -T 14 -p -o "$LORDECCOROUTFILEK19"
lordec-correct -i "$LORDECCOROUTFILEK19" -2 "$ILLUMINASPADESCOR" -k 31 -s 3 -T 14 -p -o "$LORDECCOROUTFILEK31"
lordec-correct -i "$LORDECCOROUTFILEK31" -2 "$ILLUMINASPADESCOR" -k 41 -s 3 -T 14 -p -o "$LORDECCOROUTFILEK41"

printf "\n"
printf "               Module 3: Flye de novo Assembly with Long-Reads\n"
printf "#####################################################################\n"
printf "\n\n"

#0 Module 3
# Flye Genome Assembly of the LoRDEC Corrected Long Reads
flye --iterations 3 --nano-corr "$LORDECCOROUTFILEK41" -g 6.5m --out-dir "$FLYEOUT" -t 14

printf "\n"
printf "               Module 4: Final Unicycler Hybrid Assembly\n"
printf "#####################################################################\n"
printf "\n\n"

# Module 4
# Unicycler Hybrid Assembly with Corrected Illumina and Nanopore Reads py36
FLYEGFA="$FLYEOUT"assembly_graph.gfa
mkdir "$UNICYCLEROUT"
#FLYEASSEMBLY="$FLYEOUT"
unicycler -1 "$ILLUMINAR1" -2 "$ILLUMINAR2" --existing_long_read_assembly "$FLYEGFA" -l "$LORDECCOROUTFILEK41" --verbosity 2 -t 14 --keep 2 -o "$UNICYCLEROUT"

cp $UNICYCLEROUT"assembly.fasta" "$STRAIN".hybrid.complete.fasta

# Exit HybridAsBro Environment
conda deactivate

printf "Thank you for using the Hybrid Assembly Pipeline :D\n"
printf "If  you found this pipeline helpful please consider citing:\n"
printf "				TBA							\n"

printf "K bye.\n\n"

# Exit Pipeline
exit 0
