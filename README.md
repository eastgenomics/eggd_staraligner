# eggd_staraligner

## What does this app do?
Runs STAR-Aligner to align RNA sequence data to a reference 

## What inputs are required for this app to run?
* `--fastqs`: Array of gzipped FASTQs for one sample
* `--sentieon_tar`: Tarballed Sentieon package. Currently defaults to use Sentieon 202112.05
* `--genome_indices`: Tarballed genome indices for reference genome. Currently defaults to use indices generated using GENCODE v41 .gtf
* `--reference_genome`: Tarballed GRCh38 reference genome FASTA + index. Current defaults to use GRCh38.no_alt_analysis_set_chr_mask21.fasta-index.tar.gz in 001_Reference

## How does this app work?
eggd_staraligner takes an input array of fastq.gz files, with multiple for each R1 and R2. These are concatenated into one R1 fastq.gz and one R2 fastq.gz. STAR aligner is run on the concatenated fastqs to align them to the reference genome. 

## What does this app output?
eggd_staraligner outputs a .bam and .bam.bai file for the alignmenent of the sample and a Chimeric.out.junctions file with junction information

## Notes
* This app is not ready for production use

## This app was made by East GLH