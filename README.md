# Hybrid Assembly Pipeline
This pipeline is for hybrid assembly of microbial genomes using ONT long read and Illumina PE read data.  

# Assumptions
This pipeline assumes you have adapter trimmed and filtered the ONT long reads (min length > 1000bp). Additionaly, this pipeline assumes you have pre-treated your Illumina reads but note they will be further kmer corrected.

# Overview
The pipeline uses k-mer correction, De bruijn grah assembly with long reads, and then hybrid assembly and polishing using Unicycler:
1. Module: SPAdes k-mer Correction
    + Conduct SPAdes mediate error correction of Illumina PE data.
2. Module: LoRDEC k-mer Correction
    + Conduct LoRDEC mediated error correction of ONT long read data.
3. Module: Flye De Novo Genome Assembly
    + Complete De Bruijn graph-mediated assembly using k-mer corrected data.  
4. Module:Unicycler Assembly and Polish
    + Complete hybrid assembly, polish, and ori-identification using Flye .gfa and k-mer corrected Illumina data.

# Dependencies
For managing the dependies depends on conda, the pipeline requires the following environments.    

conda create -n HybridAsBro python=3.6
source activate HybridAsBro
conda install spades lordec unicycler
conda deactivate

conda create -n Flye python=2.7
source activate Flye
conda install
conda deactivate

# Usage

./HybridAssembly.sh [SAMPLEID].chop.filt.fastq Illumina.R1.fastq.gz Illumina.R2.fastq.gz

# Citations

TBD
