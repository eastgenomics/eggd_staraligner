# eggd_staraligner

## What does this app do?
Runs STAR-Aligner to align RNA sequence data to a reference 

## What inputs are required for this app to run?
* `--fastqs`: (array of files) Array of gzipped FASTQs for one sample
* `--sentieon_tar`: (file) Tarballed Sentieon package. Currently defaults to use Sentieon 202112.05
* `--read_length`: (int) The length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads. The reference genome and genome indices from the CTAT library needs this to be set at 150, so that is the default.
* `--genome_lib`:(file) A CTAT genome library, which is a reference file bundle required by Trinity CTAT tools. This contains the genome indices and reference genome that are needed to run Sentieon STAR aligner
* `--parameters`: (string) The parameters with which to run Sentieon's STAR aligner. Should be space-delimated, and in the format `--parameter-name value` e.g. `--alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5`. Has default parameters set if not configured by the user. This input cannot be used to set the following parameters, as they are already set within the app:
    * `--runThreadN`: the app calculates this from the spec of the machine it is running on
    * `--genomeDir`: this is extracted from the `--genome_lib` input
    * `--readFilesCommand`: hardcoded in the app as `zcat` as the fastqs have to be .gz compressed
    * `--readFilesManifest`: this file is generated in the app
    * `--sjdbOverhang`: this is set by the input `--read_length`

## How does this app work?
eggd_staraligner takes an input array of fastq.gz files and runs sentieon STAR aligner on the fastqs to align them to the reference genome. 

## What does this app output?
eggd_staraligner outputs a .bam and .bam.bai file for the alignmenent of the sample and a Chimeric.out.junctions file with junction information. It also outputs several Sentieon STAR aligner log files

## This app was made by East GLH