#!/bin/bash


echo "Pipeline started..."


# step 1
cd ../reference
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
salmon index -t ggal_1_48850000_49020000.Ggal71.500bpflank.fa -i out_index &>log_index
cd ../data
echo " indexing over"


# step 2
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
salmon quant --libType=U -i ../reference/out_index -1 ggal_gut_1.fq -2 ggal_gut_2.fq -o ggal_gut &>log_quant
echo " quantification over"


# step 3
mkdir -p out_fastqc
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
fastqc -o out_fastqc -f fastq -q ggal_gut_1.fq ggal_gut_2.fq &>log_fq
echo " quality control over"


# step 4
mkdir -p out_multiqc
cd out_multiqc
ln -s ../ggal_gut .
ln -s ../out_fastqc .
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
multiqc -v . &>../log_mq
cd ..
echo " multiple quality control over"


echo "Pipeline finished.  Outcomes not checked, up to you."
