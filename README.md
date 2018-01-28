# Introduction

On Monday the 22th of January, my colleague [Sander Wuyts](https://github.com/swuyts) won the Davos Bitcoin Challenge issued by EMBL scientist Nick Goldman (see [here](https://swuyts.wordpress.com/2018/01/16/from-dna-to-bitcoin-how-i-won-the-davos-bitcoin-challenge/) for his story). The scripts he used to get from raw sequencing reads to the files encoded in the DNA are available [here](https://github.com/swuyts/Davos); this repository was cloned from there. Sander did an amazing job decoding the files hidden in the DNA, but some of his scripts are rather slow because they were written in R. In this forked repository, I added a small python tool that can perform assembly of the derandomized reads using a number of different strategies.

# Assembly strategies

In the simplest case, `assembler.py` can assemble the reads using a **consensus per position** strategy. This means that for each position in the assembly, all reads that overlap with this position are considered, the nucleotides corresponding to this position are extracted from the reads, their frequencies are multiplied by the frequencies of the reads they come from and a consensus nucleotide is determined by majority vote. 

```bash
./assembler.py \
  --fin_sequences seqs.tsv \
  --dout_assemblies assemblies_new
  ```
  
To make it more interesting, a **dynamic programming strategy** is also possible. In this case, the optimal read for each index will first be deteremined, after which assembly is performed by consensus per position. "Optimal" in this case means leading to the set of reads with the least number of total nucleotide conflicts between overlapping reads. 

```bash
./assembler.py \
  --fin_sequences seqs.tsv \
  --dout_assemblies assemblies_new \
  --dynamic_programming
  ```
  
If there are some indices with a very large number of unique reads, the dynamic programming approach might get slow. To speed up the process, it is also possible to **first select the n most frequent reads per position**:

```bash
./assembler.py \
  --fin_sequences seqs.tsv \
  --dout_assemblies assemblies_new \
  --max_sequences_per_index 10 \
  --dynamic_programming
  ```
