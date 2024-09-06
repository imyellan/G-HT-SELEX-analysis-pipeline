#!/bin/bash

# emulating https://link.springer.com/article/10.1186/s13059-020-01996-3#Sec16/,
# the MEX paper's motif derivation method, etc etc

# 1. de-duplicate reads

# 2. split into train-test sets randomly, 70-30

# 3. with the training set, use the MEX guys kmer enrichment script to extract 
# the top 10000 reads with the highest kmer enrichment
~/software/HT-SELEX-kmer-filtering/extract_topk.py --fastq --out...

# 4. use the top 10000 reads to derive motifs with the various tools

# 5. TBD testing on the test set