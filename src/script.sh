#!/bin/bash

set -exo pipefail #if any part goes wrong, job will fail

dx-download-all-inputs # download inputs from json

# resources folder will not exist on worker, so removed here.
mkdir /home/dnanexus/fastqs
mkdir /home/dnanexus/genomeDir
mkdir /home/dnanexus/reference_genome
mkdir -p /home/dnanexus/out/output_bam
mkdir /home/dnanexus/out/output_bam_bai
mkdir /home/dnanexus/out/chimeric_junctions
mkdir /home/dnanexus/R1/
mkdir /home/dnanexus/R2/

# Unpack tarred files 
tar xvzf /home/dnanexus/in/sentieon_tar/sentieon-genomics-*.tar.gz -C /usr/local
tar xvzf /home/dnanexus/in/genome_indices/*.tar.gz -C /home/dnanexus/genomeDir
tar xvzf /home/dnanexus/in/reference_genome/*tar.gz -C /home/dnanexus/reference_genome

# Move all the fastqs from subdirectories into one directory
find ~/in/fastqs -type f -name "*" -print0 | xargs -0 -I {} mv {} ~/fastqs

source /home/dnanexus/license_setup.sh # run license setup script

export SENTIEON_INSTALL_DIR=/usr/local/sentieon-genomics-*

SENTIEON_BIN_DIR=$(echo $SENTIEON_INSTALL_DIR/bin)

export PATH="$SENTIEON_BIN_DIR:$PATH"

# Concatenate the R1 and R2 files over multiple lanes into one file
cd /home/dnanexus/fastqs  # Move into fastqs directory to list fastqs
R1=($(ls *R1*))
R2=($(ls *R2*))

### Tests
# Check that there are the same number of files in each list
# There should be an equal number of R1 and R2 files
if [[ ${#R1[@]} -ne ${#R2[@]} ]]
  then echo "The number of R1 and R2 files for this sample are not equal"
  exit 1
fi

# Check R1 and R2 are paired correctly; for each R1 there should be a matching R2
# Create test arrays that are equal to the arrays for R1 and R2
R1_test=${R1[@]}
R2_test=${R2[@]}

# Define strings to remove from file name in test arrays
to_cut_R1="R1"
to_cut_R2="R2"
cut_fastq=".fastq.gz"

# Remove "R1" and "R2" and the file suffix ".fastq.gz" from all file names
for i in "${!R1_test[@]}"; do
  R1_test[$i]=${R1_test[$i]//$to_cut_R1/};
  R1_test[$i]=${R1_test[$i]//$cut_fastq/}
done

for i in "${!R2_test[@]}"; do
  R2_test[$i]=${R2_test[$i]//$to_cut_R2/};
  R2_test[$i]=${R2_test[$i]//$cut_fastq/}
done

# Test that when "R1" and "R2" are removed the two arrays have identical file names
for i in "${!R1_test[@]}"; do
  if [[ ! "${R2_test}" =~ "${R1_test[$i]}" ]];
  then echo "Each R1 FASTQ does not appear to have a matching R2 FASTQ"
  exit 1
  fi
done

# Extract sample name from input FASTQ file names
sample_name=$(echo $R1[0] | cut -d '_' -f 1)

# Concatenate all R1 files into single concatenated R1 file
for i in "${!R1[@]}"; do
  cat "${R1[$i]}" >> /home/dnanexus/R1/"${sample_name}_R1_concat.fastq.gz"
done

# Concatenate all R2 files into single concatenated R2 file
for i in "${!R2[@]}"; do
  cat "${R2[$i]}" >> /home/dnanexus/R2/"${sample_name}_R2_concat.fastq.gz"
done

cd /home/dnanexus

# Run STAR-aligner
NUMBER_THREADS=32
export STAR_REFERENCE=/home/dnanexus/genomeDir/output # Reference transcripts
export REFERENCE=/home/dnanexus/reference_genome/*.fa # Reference genome, standard GRCh38
SAMPLE=/home/dnanexus/R1/${sample_name}_R1_concat.fastq.gz
SAMPLE2=/home/dnanexus/R2/${sample_name}_R2_concat.fastq.gz
GROUP_NAME=${read_group}
SAMPLE_NAME=${sample_name}
PLATFORM=ILLUMINA
READ_LENGTH_MINUS_1=100
SORTED_BAM="/home/dnanexus/out/${sample_name}.bam"

sentieon STAR --runThreadN ${NUMBER_THREADS} --genomeDir ${STAR_REFERENCE} \
    --readFilesIn ${SAMPLE} ${SAMPLE2} --readFilesCommand "zcat" \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outBAMcompression 0 \
    --outSAMattrRGline ID:${GROUP_NAME} SM:${SAMPLE_NAME} PL:${PLATFORM} \
    --twopassMode Basic --twopass1readsN -1 --sjdbOverhang ${READ_LENGTH_MINUS_1} \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 8 \
    --chimOutJunctionFormat 1 \
    | sentieon util sort -r ${REFERENCE} -o ${SORTED_BAM} -t ${NUMBER_THREADS} -i -

# Move output files to /out directory so they will be uploaded
mv /home/dnanexus/out/${sample_name}.bam /home/dnanexus/out/output_bam
mv /home/dnanexus/out/${sample_name}.bam.bai /home/dnanexus/out/output_bam_bai
mv /home/dnanexus/Chimeric.out.junction /home/dnanexus/out/chimeric_junctions

dx-upload-all-outputs