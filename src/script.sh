#!/bin/bash

set -exo pipefail #if any part goes wrong, job will fail

dx-download-all-inputs # download inputs from json

mkdir /home/dnanexus/fastqs
mkdir /home/dnanexus/genomeDir
mkdir /home/dnanexus/genome_lib
mkdir /home/dnanexus/reference_genome
mkdir -p /home/dnanexus/out/output_bam
mkdir /home/dnanexus/out/output_bam_bai
mkdir /home/dnanexus/out/chimeric_junctions
mkdir /home/dnanexus/out/output_mark_duplicates_bam
mkdir /home/dnanexus/out/logs

# Unpack tarred files 
tar xvzf /home/dnanexus/in/genome_lib/*.tar.gz -C /home/dnanexus/genome_lib
tar xvzf /home/dnanexus/in/sentieon_tar/sentieon-genomics-*.tar.gz -C /usr/local

# Extract CTAT library filename
lib_dir=$(find /home/dnanexus/genome_lib -type d -name "*" -mindepth 1 -maxdepth 1 | rev | cut -d'/' -f-1 | rev)

# Move genome indices and reference genome to specific folders
mv /home/dnanexus/genome_lib/${lib_dir}/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/* /home/dnanexus/genomeDir/
mv /home/dnanexus/genome_lib/${lib_dir}/ctat_genome_lib_build_dir/ref_genome.fa /home/dnanexus/reference_genome
mv /home/dnanexus/genome_lib/${lib_dir}/ctat_genome_lib_build_dir/ref_genome.fa.fai /home/dnanexus/reference_genome

# Move all the fastqs from subdirectories into one directory
find ~/in/fastqs -type f -name "*" -print0 | xargs -0 -I {} mv {} ~/fastqs

source /home/dnanexus/license_setup.sh # run license setup script

export SENTIEON_INSTALL_DIR=/usr/local/sentieon-genomics-*

SENTIEON_BIN_DIR=$(echo $SENTIEON_INSTALL_DIR/bin)

export PATH="$SENTIEON_BIN_DIR:$PATH"

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
# Remove "R1" and "R2" and the file suffix from all file names
_trim_fastq_endings () {
  # Takes array of fastq files with their read number ("R1" or "R2"), 
  # trims the endings off every file, and returns as an array
  local fastq_array=("$@")
  # Define strings to remove from file name in test arrays
  # This app can take fastqs with .fq.gz and .fastq.gz suffixes so need to
  # identify which suffix the input files have
  # The variable fastq_suffix is exported so it can be used later in the script
  local read_to_cut=$1
  if [[ "${fastq_array[1]}" == *".fastq.gz" ]]; then
    fastq_suffix=".fastq.gz"
    export fastq_suffix
  elif [[ "${fastq_array[1]}" == *".fq.gz" ]]; then
    fastq_suffix=".fq.gz"
    export fastq_suffix
  else
    echo "Suffixes of fastq files not recognised as .fq.gz or .fastq.gz"
    exit 1
  fi
  
  for i in "${!fastq_array[@]}"; do
    fastq_array[$i]=${fastq_array[$i]//$read_to_cut/};
    fastq_array[$i]=${fastq_array[$i]//$fastq_suffix/};
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

# Test that there are no banned parameters in --parameters input string
banned_parameters=(--runThreadN --genomeDir --readFilesIn --readFilesCommand --readFilesManifest --outSAMattrRGline --sjdbOverhang --outStd --outSAMtype)
for parameter in ${banned_parameters[@]}; do
  if [[ "$parameters" == *"$parameter"* ]]; then
    echo "Ihe parameter ${parameter} was set as an input. This parameter is set within the app and cannot be set as an input. Please repeat without this parameter"
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
fq_arr=($(ls *$fastq_suffix)) # ls command is alphabetical so R1 should be before R2

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
    arr=($(ls *$fastq_suffix))
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

# STAR can run the aligner with the manifest and the input --readFilesManifest. Future development could implement this,
# however this version of eggd_staraligner will extract the read groups from the manifest and pass them to the aligner manually.

# Extract newline separated list of read groups from manifest.tsv
read_group_list=$(while IFS=' ' read -r R1 R2; do
    lane=$( echo $R1 | cut -d'_' -f3)  # The lanes are not necessarily L001 and L002, so this grabs the lane information from the fastq file name
    rg=$(grep _${lane}_ manifest.tsv | cut -f3-)  # Extract the read group
    echo "$rg"
done < file.txt)

# Convert the read groups into a space-comma-space delimited string so they can be passed to the --outSAMattrRGline input to the aligner
readarray -t read_groups <<< "$read_group_list"              # Convert to array
printf -v read_groups_delimited ' , %s' "${read_groups[@]}"  # Convert to space-comma-space delimited string
read_groups_delimited=${read_groups_delimited:3}             # Remove leading space-comma-space
echo -e  $read_groups_delimited

cd /home/dnanexus

# Run STAR-aligner
NUMBER_THREADS=${INSTANCE##*_x}
export STAR_REFERENCE=/home/dnanexus/genomeDir # Reference transcripts have been moved from the CTAT lib to here
export REFERENCE=/home/dnanexus/reference_genome/*.fa # Reference genome has been moved from the CTAT lib to here
SORTED_BAM="/home/dnanexus/out/${sample_name}.star.bam"
CTAT_GENOME_INDICES_READ_LENGTH_MINUS_1=$((${ctat_genome_indices_read_length}-1)) # Use read_length value from input JSON. The default is 151 bp, because that is the read length used in generation of the CTAT genome library

sentieon STAR --runThreadN ${NUMBER_THREADS} \
    --genomeDir ${STAR_REFERENCE} \
    --readFilesIn ${R1_list} ${R2_list} \
    --outSAMattrRGline ${read_groups_delimited} \
    --readFilesCommand "zcat" \
    --outStd BAM_Unsorted \
    --outSAMtype BAM Unsorted \
    --sjdbOverhang ${CTAT_GENOME_INDICES_READ_LENGTH_MINUS_1} \
    ${parameters} \
    | sentieon util sort -r ${REFERENCE} -o ${SORTED_BAM} -t ${NUMBER_THREADS} -i -

# Take the output bam file and run the STAR command to mark the duplicates. This generates a .mark_duplicates.star.Processed.out.bam file with the duplicates marked
sentieon STAR --runMode inputAlignmentsFromBAM --limitBAMsortRAM 32000000000 --inputBAMfile ${SORTED_BAM} --bamRemoveDuplicatesType UniqueIdentical --outFileNamePrefix /home/dnanexus/out/${sample_name}.mark_duplicates.star.

# Move output files to /out directory so they will be uploaded
mv /home/dnanexus/out/${sample_name}.star.bam /home/dnanexus/out/output_bam
mv /home/dnanexus/out/${sample_name}.star.bam.bai /home/dnanexus/out/output_bam_bai
mv /home/dnanexus/Chimeric.out.junction /home/dnanexus/out/chimeric_junctions/${sample_name}.chimeric.out.junction
mv /home/dnanexus/out/${sample_name}.mark_duplicates.star.Processed.out.bam /home/dnanexus/out/output_mark_duplicates_bam
mv /home/dnanexus/Log* /home/dnanexus/out/logs

dx-upload-all-outputs
