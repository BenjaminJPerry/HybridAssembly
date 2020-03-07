# Hybrid Assembly Pipeline
This pipeline is for hybrid assembly of microbial genomes using ONT long read and Illumina PE data.  

# Assumptions
This pipeline assumes you have adapter trimmed and filtered the ONT long reads (min length > 1000bp).  
Additionaly, this pipeline assumes you have pre-treated your Illumina reads, but note they will be further kmer corrected using BayesHammer.

# Overview
1. Module: SPAdes k-mer Correction
    + SPAdes BayesHammer error correction of Illumina PE data.
    + SPAdes k-mer corrected reads can be found in SPAdes_Cor/corrected/ for alignment.

2. Module: Iterative LoRDEC Long-Read k-mer Correction with Increasing Kmer Size
    + LoRDEC error correction of ONT long read data with kmer corrected Illumina PE reads.
    + LoRDEC corrected long reads remain available for alignment in the LoRDEC subdirectory.

3. Module: Flye de novo Assembly with Long-Reads
    + de novo assembly with k-mer corrected ONT reads to generate assembly graph for Unicycler.
    + Flye iteratively polishes final assembly which remains available in the Flye subdirectory (expect SNVs and miss-assemblies).

4. Module: Final Unicycler Hybrid Assembly
    + Unicycler SPAdes de novo assembly of Illumina reads is closed/corrected with Flye de novo assembly graph.
    + Unicycler built in Pilon polish with chromosomal rotation and iteration corrects consensus assembly.

# Dependencies
Pipeline requires the following environments be created.    
```bash  
conda create -n HybridAsBro python=3.6  
source activate HybridAsBro  
conda install spades lordec flye unicycler  
conda deactivate  
```  
Additionally, you will need to increase the memmory usage allowed by Pilon.  
To do so, modify the following line in '/home/\<user\>/miniconda3/env/HybridAsBro/bin/pilon'  
```bash
default_jvm_mem_opts = ['-Xms512m', '-Xmx1g']
```
To  
```bash
default_jvm_mem_opts = ['-Xms1g', '-Xmx12g']
```
# Usage
```bash  
./HybridAssembly.sh [SAMPLEID].chop.filt.fastq Illumina.R1.fastq.gz Illumina.R2.fastq.gz
```  

# Citations

TBD
