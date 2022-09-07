#!/bin/bash

set -exo pipefail #if any part goes wrong, job will fail

dx-download-all-inputs # download inputs from json
#dx download "$left_fastq" -o left_fastq
#dx download "$right_fastq" -o right_fastq

# resources folder will not exist on worker, so removed here.
mkdir /home/dnanexus/genomeDir
mkdir /home/dnanexus/reference_genome
mkdir -p /home/dnanexus/out/output_bam
mkdir /home/dnanexus/out/output_bam_bai

tar xvzf /home/dnanexus/in/sentieon_tar/sentieon-genomics-*.tar.gz -C /usr/local # unpack tar
tar xvzf /home/dnanexus/in/genome_indexes/*.tar.gz -C /home/dnanexus/genomeDir #transcirpt data from that release of gencode
tar xvzf /home/dnanexus/in/reference_genome/*tar.gz -C /home/dnanexus/reference_genome

source /home/dnanexus/license_setup.sh # run license setup script

export SENTIEON_INSTALL_DIR=/usr/local/sentieon-genomics-*

SENTIEON_BIN_DIR=$(echo $SENTIEON_INSTALL_DIR/bin)

export PATH="$SENTIEON_BIN_DIR:$PATH"



NUMBER_THREADS=32
export STAR_REFERENCE=/home/dnanexus/genomeDir/output # Reference transcripts
echo $STAR_REFERENCE
export REFERENCE=/home/dnanexus/reference_genome/*.fa # Reference genome, standard GRCh38
echo $REFERENCE
SAMPLE=/home/dnanexus/in/left_fq/*.fastq.gz
SAMPLE2=/home/dnanexus/in/right_fq/*.fastq.gz
GROUP_NAME="test_group"
SAMPLE_NAME="test"
PLATFORM=ILLUMINA
READ_LENGTH_MINUS_1=100
SORTED_BAM='/home/dnanexus/out/bam_file'

sentieon STAR --runThreadN ${NUMBER_THREADS} --genomeDir ${STAR_REFERENCE} \
    --readFilesIn ${SAMPLE} ${SAMPLE2} --readFilesCommand "zcat" \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outBAMcompression 0 \
    --outSAMattrRGline ID:${GROUP_NAME} SM:${SAMPLE_NAME} PL:${PLATFORM} \
    --twopassMode Basic --twopass1readsN -1 --sjdbOverhang ${READ_LENGTH_MINUS_1} \
    | sentieon util sort -r ${REFERENCE} -o ${SORTED_BAM} -t ${NUMBER_THREADS} -i -

mv /home/dnanexus/out/bam_file /home/dnanexus/out/output_bam
mv /home/dnanexus/out/bam_file.bai /home/dnanexus/out/output_bam_bai

dx-upload-all-outputs