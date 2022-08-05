#!/bin/bash

set -exo pipefail #if any part goes wrong, job will fail

dx-download-all-inputs # download inputs from json

tar xvzf /home/dnanexus/in/sentieon_tar/sentieon-genomics-*.tar.gz -C /usr/local # unpack tar

source /home/dnanexus/license_setup.sh # run license setup script

export SENTIEON_INSTALL_DIR=/usr/local/sentieon-genomics-*

SENTIEON_BIN_DIR=$(echo $SENTIEON_INSTALL_DIR/bin)

export PATH="$SENTIEON_BIN_DIR:$PATH"

# resources folder will not exist on worker, so removed here.

# Quick start package has

#   sentieon_quickstart.sh: the sample shell script that drives the entire pipeline.
#   reference: a directory that contains human genome reference files and database files of known SNP sites.
#   FASTQ files: sample sequence files.

#tar xzvf quick_start.tar.gz 
#sh sentieon_quickstart.sh &

sentieon STAR --runThreadN NUMBER_THREADS --genomeDir STAR_REFERENCE \
    --readFilesIn SAMPLE SAMPLE2 --readFilesCommand "zcat" \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outBAMcompression 0 \
    --outSAMattrRGline ID:GROUP_NAME SM:SAMPLE_NAME PL:PLATFORM \
    --twopassMode Basic --twopass1readsN -1 --sjdbOverhang READ_LENGTH_MINUS_1 \
    | sentieon util sort -r REFERENCE -o SORTED_BAM -t NUMBER_THREADS -i -