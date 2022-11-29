# eggd_staraligner

## What does this app do?
Runs Sentieon STAR-Aligner to align RNA sequence data to a reference 

## What inputs are required for this app to run?
* `--fastqs`: (array of files) Array of gzipped FASTQs for one sample
* `--sentieon_tar`: (file) Tarballed Sentieon package. Currently defaults to use Sentieon 202112.05
* `--ctat_genome_indices_read_length`: (int) The length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads used in construction of the indices. The reference genome and genome indices from the CTAT library used reads of length 151 bp, so this is the default
* `--genome_lib`:(file) A CTAT genome library, which is a reference file bundle required by Trinity CTAT tools. This contains the genome indices and reference genome that are needed to run Sentieon STAR aligner
* `--parameters`: (string) The parameters with which to run Sentieon STAR aligner. Should be space-delimited, and in the format `--parameter-name value` e.g. `--alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5`. Has default parameters set if not configured by the user. This input cannot be used to set the following parameters, as they are already set within the app:
    * `--runThreadN`: the app calculates this from the spec of the machine it is running on
    * `--genomeDir`: this is extracted from the `--genome_lib` input
    * `--readFilesIn`: the app sets this using the files in the `--fastqs` input
    * `--readFilesCommand`: hardcoded in the app as `zcat` as the fastqs have to be .gz compressed
    * `--readFilesManifest`: the manifest file is generated in the app
    * `--sjdbOverhang`: this is set by the input `--ctat_genome_indices_read_length`
    * `--outStd`: this is set in the app to ensure the bam file is the main output
    * `--outSAMtype`: this is set in the app to ensure it outputs a bam file

## How does this app work?
eggd_staraligner takes an input array of fastq.gz files and runs sentieon STAR aligner on the fastqs to align them to the reference genome. 

## What does this app output?
eggd_staraligner outputs a .bam and .bam.bai file for the alignmenent of the sample and a Chimeric.out.junctions file with junction information. It also outputs several Sentieon STAR aligner log files

## This app was made by East GLH