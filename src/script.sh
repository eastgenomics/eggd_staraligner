#!/bin/bash

set -exo pipefail #if any part goes wrong, job will fail

dx-download-all-inputs # download inputs from json
#dx download "$left_fastq" -o left_fastq
#dx download "$right_fastq" -o right_fastq

# path2fastqs=project-GFfV0qQ42QQVf6VX12zxxvx3:/eggd_staraligner/test_fastqs

# for sample in $(dx ls ${path2fastqs}/*_L001_R1_001.fastq.gz);do echo $sample |  cut -d '_' -f 1 | grep -v Undetermined;done > rna_samples

# for sample in $(less rna_samples);
# do dx run applet-G91YBp84fVkjv3GyK2219y22 \
# -ifiles=${path2fastqs}/${sample}*_L001_R1_001.fastq.gz \
# -ifiles=${path2fastqs}/${sample}*_L002_R1_001.fastq.gz \
# -ifiles2=${path2fastqs}/${sample}*_L001_R2_001.fastq.gz \
# -ifiles2=${path2fastqs}/${sample}*_L002_R2_001.fastq.gz \
# -ioutput_filename=${sample}_R1_concat.fastq.gz \
# -ioutput_filename2=${sample}_R2_concat.fastq.gz \
# --destination=/output/concatenated_fastqs -y;
# done

# resources folder will not exist on worker, so removed here.
mkdir /home/dnanexus/genomeDir
mkdir /home/dnanexus/reference_genome
mkdir -p /home/dnanexus/out/output_bam
mkdir /home/dnanexus/out/output_bam_bai
mkdir /home/dnanexus/out/chimeric_junctions

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
SORTED_BAM='/home/dnanexus/out/bam_file.bam'

sentieon STAR --runThreadN ${NUMBER_THREADS} --genomeDir ${STAR_REFERENCE} \
    --readFilesIn ${SAMPLE} ${SAMPLE2} --readFilesCommand "zcat" \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outBAMcompression 0 \
    --outSAMattrRGline ID:${GROUP_NAME} SM:${SAMPLE_NAME} PL:${PLATFORM} \
    --twopassMode Basic --twopass1readsN -1 --sjdbOverhang ${READ_LENGTH_MINUS_1} \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 8 \
    --chimOutJunctionFormat 1 \
    | sentieon util sort -r ${REFERENCE} -o ${SORTED_BAM} -t ${NUMBER_THREADS} -i -

mv /home/dnanexus/out/bam_file.bam /home/dnanexus/out/output_bam
mv /home/dnanexus/out/bam_file.bam.bai /home/dnanexus/out/output_bam_bai
mv /home/dnanexus/Chimeric.out.junction /home/dnanexus/out/chimeric_junctions

dx-upload-all-outputs