#!/bin/bash

salmon_image="/data/bio-intro-containers/exercises/pipeline/data/salmon_1.2.1--hf69c8f4_0.sif"
fastqc_image="/data/bio-intro-containers/exercises/pipeline/data/fastqc_0.11.9--0.sif"
multiqc_image="/data/bio-intro-containers/exercises/pipeline/data/multiqc_1.9--pyh9f0ad1d_0.sif"

export SINGULARITY_BINDPATH="/data"


echo "Pipeline started..."


# step 1
cd ../reference
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
singularity exec \
  $salmon_image \
  salmon index -t ggal_1_48850000_49020000.Ggal71.500bpflank.fa -i out_index &>log_index
cd ../data
echo " indexing over"


# step 2
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
singularity exec \
  $salmon_image \
  salmon quant --libType=U -i ../reference/out_index -1 ggal_gut_1.fq -2 ggal_gut_2.fq -o ggal_gut &>log_quant
echo " quantification over"


# step 3
mkdir -p out_fastqc
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
singularity exec \
  $fastqc_image \
  fastqc -o out_fastqc -f fastq -q ggal_gut_1.fq ggal_gut_2.fq &>log_fq
echo " quality control over"


# step 4
mkdir -p out_multiqc
cd out_multiqc
ln -s ../ggal_gut .
ln -s ../out_fastqc .
# ONLY CHANGE THE NEXT LINE - EXECUTION LINE
singularity exec \
  $multiqc_image \
  multiqc -v . &>../log_mq
cd ..
echo " multiple quality control over"


echo "Pipeline finished.  Outcomes not checked, up to you."
