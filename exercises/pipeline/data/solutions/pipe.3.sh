#!/bin/bash

PATH="$WORK/exercises/pipeline/data/wrappers:$PATH"
export SINGULARITY_BINDPATH="$WORK/exercises/pipeline"


echo "Pipeline started..."


# step 1
cd ../reference
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
salmon index -t ggal_1_48850000_49020000.Ggal71.500bpflank.fa -i out_index &>log_index
cd ../data
echo " indexing completed"


# step 2
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
salmon quant --libType=U -i ../reference/out_index -1 ggal_gut_1.fq -2 ggal_gut_2.fq -o ggal_gut &>log_quant
echo " quantification completed"


# step 3
mkdir out_fastqc
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
fastqc -o out_fastqc -f fastq -q ggal_gut_1.fq ggal_gut_2.fq &>log_fq
echo " quality control completed"


# step 4
mkdir out_multiqc
cd out_multiqc
ln -s ../ggal_gut .
ln -s ../out_fastqc .
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
multiqc -v . &>../log_mq
cd ..
echo " multiple quality control completed"


echo "Pipeline finished!"
