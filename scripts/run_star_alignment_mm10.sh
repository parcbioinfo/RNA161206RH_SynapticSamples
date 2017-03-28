#!/bin/bash
#Priscila Darakjian
#Script for Alignment with STAR 
# STAR version - 2.5.2b
# samtools version - 1.4

for file in *.fastq.gz
do
  file_path="./"$file
  echo $file_path
  STAR --genomeDir /lawrencedata/sharedRNASeq/Mouse_reference_genome/unmasked/genome_for_STAR_mm10/ --readFilesIn $file_path --readFilesCommand zcat --outFileNamePrefix ${file}_ --runThreadN 12 --outFilterMismatchNmax 3 --outFilterMultimapNmax 1             
  samtools view -bS ${file}_Aligned.out.sam > ${file}_Aligned.out.bam
  samtools sort ${file}_Aligned.out.bam -o ${file}_Aligned_sorted.bam
  samtools index ${file}_Aligned_sorted.bam
  file_path=""
 done
