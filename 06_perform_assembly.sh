#!/bin/bash

./assembler.py \
  --fin_sequences seqs.tsv \
  --dout_assemblies assemblies_new \
  --max_sequences_per_index 10 \
  --dynamic_programming

