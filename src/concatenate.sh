path2fastqs=project-GFfV0qQ42QQVf6VX12zxxvx3:/eggd_staraligner/test_fastqs

for sample in $(dx ls ${path2fastqs}/*_L001_R1_001.fastq.gz);do echo $sample |  cut -d '_' -f 1 | grep -v Undetermined;done > rna_samples

for sample in $(less rna_samples);
do dx run applet-G91YBp84fVkjv3GyK2219y22 \
-ifiles=${path2fastqs}/${sample}*_L001_R1_001.fastq.gz \
-ifiles=${path2fastqs}/${sample}*_L002_R1_001.fastq.gz \
-ifiles2=${path2fastqs}/${sample}*_L001_R2_001.fastq.gz \
-ifiles2=${path2fastqs}/${sample}*_L002_R2_001.fastq.gz \
-ioutput_filename=${sample}_R1_concat.fastq.gz \
-ioutput_filename2=${sample}_R2_concat.fastq.gz \
--destination=/output/concatenated_fastqs -y;
done