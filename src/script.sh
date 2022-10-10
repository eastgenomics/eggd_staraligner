#!/bin/bash

set -exo pipefail #if any part goes wrong, job will fail

dx-download-all-inputs # download inputs from json

mkdir /home/dnanexus/fastqs
mkdir /home/dnanexus/genomeDir
mkdir /home/dnanexus/reference_genome
mkdir -p /home/dnanexus/out/output_bam
mkdir /home/dnanexus/out/output_bam_bai
mkdir /home/dnanexus/out/chimeric_junctions

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

# The LD_PRELOAD environment variable can be used to load the jemalloc library in Sentieon at run time.
# This escapes a repeated error where this fails to load during running of eggd_staraligner
export LD_PRELOAD=/usr/lib64/libjemalloc.so.2

cd /home/dnanexus/fastqs  # Move into fastqs directory to list fastqs
R1=($(ls *R1*))
R2=($(ls *R2*))

### Tests
## Check that there are the same number of files in each list
# There should be an equal number of R1 and R2 files
if [[ ${#R1[@]} -ne ${#R2[@]} ]]
  then echo "The number of R1 and R2 files for this sample are not equal"
  exit 1
fi

## Check that each R1 has a matching R2
# Remove "R1" and "R2" and the file suffix ".fastq.gz" from all file names
_trim_fastq_endings () {
  # Takes array of fastq files with their read number ("R1" or "R2"), 
  # trims the endings off every file, and returns as an array
  local fastq_array=("$@")
  # Define strings to remove from file name in test arrays
  local read_to_cut=$1
  cut_fastq=".fastq.gz"
  for i in "${!fastq_array[@]}"; do
    fastq_array[$i]=${fastq_array[$i]//$read_to_cut/};
    fastq_array[$i]=${fastq_array[$i]//$cut_fastq/};
  done
  echo ${fastq_array[@]}
}

R1_test=$(_trim_fastq_endings "R1" ${R1[@]})
R2_test=$(_trim_fastq_endings "R2" ${R2[@]})

# Test that when "R1" and "R2" are removed the two arrays have identical file names
for i in "${!R1_test[@]}"; do
  if [[ ! "${R2_test}" =~ "${R1_test[$i]}" ]];
  then echo "Each R1 FASTQ does not appear to have a matching R2 FASTQ"
  exit 1
  fi
done

### Running
# Extract sample name from input FASTQ file names
sample_name=$(echo $R1[0] | cut -d '_' -f 1)

# Convert array of R1 files into comma-separated list of R1 files
printf -v R1_list ',/home/dnanexus/fastqs/%s' "${R1[@]}"
R1_list=${R1_list:1}  # Remove leading comma

# Convert array of R2 files into comma-separated list of R1 files
printf -v R2_list ',/home/dnanexus/fastqs/%s' "${R2[@]}"
R2_list=${R2_list:1}  # Remove leading comma

# NUMBER_THREADS input to STAR-aligner needs the number of cores on the server node
# This can be extracted from the DNAnexus instance type
INSTANCE=$(dx describe --json $DX_JOB_ID | jq -r '.instanceType')  # Extract instance type

# --readFilesManifest input to STAR-aligner needs the read group information from the fastq
fq_arr=($(ls *fastq.gz)) # ls command is alphabetical so R1 should be before R2

# Create array of values for lane e.g. L001, L002, L003 etc.
for i in ${!fq_arr[@]};
    do fq_arr[$i]=$(cut -d'_' -f3 <<<${fq_arr[${i}]});
done

# Filter this array so it contains unique values only
IFS=" " read -r -a fq_arr <<< "$(tr ' ' '\n' <<< "${fq_arr[@]}" | sort -u | tr '\n' ' ')"

# Check all the fastqs and sort by lane
_check_for_string () {
    local lane_to_check=$1
    local arr=()
    arr=($(ls *fastq.gz))
    for fq in ${arr[@]}; do
        if [[ ${fq} == *$lane_to_check* ]];
            then echo ${fq}
        fi
    done
}

for L in ${fq_arr[@]}; do echo $(_check_for_string ${L}) >> file.txt; done

# Extract read group info for each lane and add to manifest
# Read group includes:
# ID: read group ID
# PL: platform, which is set to Illumina as this will use data from Illumina instruments
# SM: sample name
while IFS=' ' read -r R1 R2; do
    chopped=$(zgrep -m 1 '^@' $R1 | cut -d':' -f -4 || true)
    echo '/home/dnanexus/fastqs/'$R1 '/home/dnanexus/fastqs/'$R2 'ID:'$chopped 'PL:ILLUMINA' 'SM:'$sample_name>> newfile.tsv
done < file.txt

# Replace spaces with tabs
tr " " "\t" < newfile.tsv > manifest.tsv

cd /home/dnanexus

# Run STAR-aligner
NUMBER_THREADS=${INSTANCE##*_x}
export STAR_REFERENCE=/home/dnanexus/genomeDir/home/dnanexus/* # Reference transcripts have been unpacked here
export REFERENCE=/home/dnanexus/reference_genome/*.fa # Reference genome, standard GRCh38
SAMPLE=${R1_list}
SAMPLE2=${R2_list}
READ_LENGTH_MINUS_1=99   # 99 is recommended as the default for Illumina instruments
SORTED_BAM="/home/dnanexus/out/${sample_name}.star.bam"

sentieon STAR --runThreadN ${NUMBER_THREADS} \
    --genomeDir ${STAR_REFERENCE} \
    --readFilesIn ${SAMPLE} ${SAMPLE2} \
    --readFilesCommand "zcat" \
    --outStd BAM_Unsorted \
    --outSAMtype BAM Unsorted \
    --outBAMcompression 0 \
    --readFilesManifest "/home/dnanexus/fastqs/manifest.tsv" \
    --twopassMode Basic \
    --twopass1readsN -1 \
    --sjdbOverhang ${READ_LENGTH_MINUS_1} \
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 8 \
    --chimOutJunctionFormat 1 \
    | sentieon util sort -r ${REFERENCE} -o ${SORTED_BAM} -t ${NUMBER_THREADS} -i -

# Move output files to /out directory so they will be uploaded
mv /home/dnanexus/out/${sample_name}.star.bam /home/dnanexus/out/output_bam
mv /home/dnanexus/out/${sample_name}.star.bam.bai /home/dnanexus/out/output_bam_bai
mv /home/dnanexus/Chimeric.out.junction /home/dnanexus/out/chimeric_junctions/${sample_name}.chimeric.out.junction

dx-upload-all-outputs