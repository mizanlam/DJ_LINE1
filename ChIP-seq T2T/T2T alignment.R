### Raw Fastq files were retrieved from Zhang et al. (Submitted) 2024 ### 
# Aligning Histone ChIP-seq from naive and primed hESC to T2T for DJ analyses #


# alignment.sh

#!/bin/bash
#SBATCH --job-name=alignment
#SBATCH -N 1
#SBATCH -n 6



# Define path to genome indices
genome_index=/home/ltri/santoslab/share/genome_files/bowtie2_indexes/t2t_v1.1/chm13_v1.1


for file in *_L003_R1_001.fastq.gz; do
sample=${file%_L003_R1_001.fastq.gz}

echo $sample

bowtie2 -p 6 --local -X 1000 -x $genome_index -1 $file -2 ${sample}_L003_R2_001.fastq.gz -S ${sample}.sam 2>${sample}.align.stats

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

# Convert to BAM, sort
samtools view -bu ${sample}removedchrs.sam | samtools sort -o sorted_${sample}.bam

# Remove duplicates - don't need any of these files later, so name intermediate files redundantly
samtools collate -o namecollate.bam sorted_${sample}.bam
samtools fixmate -m namecollate.bam fixmate.bam
samtools sort -o positionsort.bam fixmate.bam
samtools markdup -rs positionsort.bam filtered_${sample}.bam
#markdup flag the duplicates and -r argument remove them

# Index bam
samtools index filtered_${sample}.bam

# For completeness sake, get a line count of the number of reads passing all filters, then remove the file
samtools view filtered_${sample}.bam -o filtered_${sample}.sam
echo "wc -l filtered"
wc -l filtered_${sample}.sam
rm filtered_${sample}.sam

# Remove intermediate sam files	
# Only uncomment this if you are sure this script is working for your samples; this is just here to save space
# rm ${sample}.sam
# rm ${sample}removedchrs.sam 	
rm namecollate.bam
rm fixmate.bam
rm positionsort.bam

done

# Compile alignment statistics with multiqc
multiqc . --ignore "FastQC/"
exit 0



#!/bin/bash
#SBATCH --job-name=bigwig

bamCoverage -b filtered_RSeT_1_H3K4me3_S8.bam -o RSeT_1_H3K4me3_T2T.bw -v --normalizeUsing RPKM
bamCoverage -b filtered_RSeT_2_H3K4me3_S8.bam -o RSeT_2_H3K4me3_T2T.bw -v --normalizeUsing RPKM

# continue for each sample ##
done
exit 0
