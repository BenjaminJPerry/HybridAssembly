# Hybrid Assembly Pipeline
This pipeline is for hybrid assembly of microbial genomes using ONT long read and Illumina PE data.  

# Assumptions
This pipeline assumes you have adapter trimmed and filtered the ONT long reads (min length > 1000bp).  
Additionaly, this pipeline assumes you have pre-treated your Illumina reads, but note they will be further kmer corrected using BayesHammer.

# Overview
1. Module: SPAdes k-mer Correction
    + SPAdes BayesHammer error correction of Illumina PE data.

2. Module: LoRDEC Long-Read k-mer Correction
    + LoRDEC error correction of ONT long read data with kmer corrected Illumina PE reads.

3. Module: Flye de novo Assembly with Long-Reads
    + de novo assembly with k-mer corrected ONT reads to generate assembly graph for Unicycler.

4. Module: Final Unicycler Hybrid Assembly
    + Unicycler SPAdes de novo assembly of Illumina reads is closed/corrected with Flye de novo assembly graph.
    + Unicycler built in Pilon polish with chromosomal rotation and iteration corrects consensus assembly.

# Dependencies
Managing the dependencies depends on conda as the pipeline requires the following environments be created.    
```bash  
conda create -n HybridAsBro python=3.6  
source activate HybridAsBro  
conda install spades lordec unicycler  
conda deactivate  

conda create -n Flye python=2.7  
source activate Flye  
conda install flye  
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
