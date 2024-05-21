# Intro 

## Generated ATAC-seq in DJ-L1 KD or Control naive hESC for T2T alignment and analysis of DJ chromatin landscape
  # Use Bowtie2 to align to genome.  
  # T2T indexes are retrieved from UCSC   
  # Deduplicate and Index with SAMtools 
  # Summarize reports with MultiQC.  


##### alignment.sh #####


#!/bin/bash
#SBATCH --job-name=alignment
#SBATCH -N 1
#SBATCH -n 6

# Define path to genome indices
genome_index=/home/ltri/santoslab/share/genome_files/bowtie2_indexes/t2t_v1.1/chm13_v1.1

# to simplify names and remove the extensions
for file in *_L002_R1_001.fastq.gz; do
sample=${file%_L002_R1_001.fastq.gz}

echo $sample

# Align with Bowtie2
  # Comma separate input files, no white space between
  # Local alignment to softclip the ends
  # Increase pair distance to 1000 to avoid the 500bp blip on fragment length plot
  # Write out alignment stats (stderr) to file for multiqc summary

bowtie2 -p 6 --local -X 1000 -x $genome_index -1 $file -2 ${sample}_L002_R2_001.fastq.gz -S ${sample}.sam 2>${sample}.align.stats

echo "wc -l aligned"
wc -l ${sample}.sam

# Count reads assigned to chrM, Un, random, then remove
echo "wc -l chrM"
grep chrM ${sample}.sam | wc -l
echo "wc -l chrUn"
grep chrUn ${sample}.sam | wc -l
echo "wc -l random"
grep random ${sample}.sam | wc -l
sed '/chrM/d;/random/d;/chrUn/d' < ${sample}.sam > ${sample}removedchrs.sam

# Convert to BAM, sort with samtools
samtools view -bu ${sample}removedchrs.sam | samtools sort -o sorted_${sample}.bam

# Remove duplicates - don't need any of these files later, so name intermediate files redundantly
samtools collate -o namecollate.bam sorted_${sample}.bam
samtools fixmate -m namecollate.bam fixmate.bam
samtools sort -o positionsort.bam fixmate.bam
samtools markdup -rs positionsort.bam filtered_${sample}.bam

# Index bam
samtools index filtered_${sample}.bam

# For completeness sake, get a line count of the number of reads passing all filters, then remove the file
samtools view filtered_${sample}.bam -o filtered_${sample}.sam
echo "wc -l filtered"
wc -l filtered_${sample}.sam
rm filtered_${sample}.sam

# Remove intermediate sam files	
rm ${sample}.sam
rm ${sample}removedchrs.sam 	
rm namecollate.bam
rm fixmate.bam
rm positionsort.bam

done

# Compile alignment statistics with multiqc
multiqc . --ignore "FastQC/"
exit 0


##### make bigwig from bam normalized by RPKM ######

#!/bin/bash
#SBATCH --job-name=bigwig

sample=${file:9}
sample2=${sample%.bam}

for file in filtered*.bam; do

bamCoverage -bam $file -o ${sample2}.bw -v -normalizeUsingRPKM

done
exit 0 


# Merge bam and bigwig

#!/bin/bash
#SBATCH --job-name=merge
#SBATCH -N 1
#SBATCH -n 6

samtools merge Control_ATAC_merge.bam filtered_Control_1_S5.bam filtered_Control_2_S6.bam filtered_Control_3_S7.bam filtered_Control_4_S8.bam

samtools merge DJL1KD_ATAC_merge.bam filtered_DJL1_1_S1.bam filtered_DJL1_2_S2.bam filtered_DJL1_3_S3.bam filtered_DJL1_4_S4.bam

for file in *merge.bam; do

samtools index $file

bamCoverage -b Control_ATAC_merge.bam -o RControl_ATAC_merge.bw -v --normalizeUsing RPKM
bamCoverage -b DJL1KD_ATAC_merge.bam -o DJL1KD_ATAC_merge.bw -v --normalizeUsing RPKM

done

exit 0
