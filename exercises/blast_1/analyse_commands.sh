#!/bin/bash
# use docker://quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4

gunzip zebrafish.1.protein.faa.gz


makeblastdb -in zebrafish.1.protein.faa -dbtype prot

blastp -query P04156.fasta -db zebrafish.1.protein.faa -out results.txt
