# Making heatmap from matrix of Control ATAC signal across ATAC-Seq and ChIP-seq peaks

#!/bin/bash
#SBATCH --job-name=deeptools
#SBATCH -N 1
#SBATCH -n 8

computeMatrix reference-point -p 8 -S Control_ATAC_merge.bw 
                                      DJL1KD_ATAC_merge.bw 
                                      RSeT_DT_H3K4me1_merge.bw 
                                      RSeT_DT_H3K4me3_merge.bw 
                                      RSeT_DT_H3K9me3_merge.bw 
                                      RSeT_DT_H3K27ac_merge.bw 
                                      RSeT_DT_H3K27me3_merge.bw 
                                      -R downDJ_log15.txt 
                                      --referencePoint center 
                                      -a 5000 
                                      -b 5000 
                                      --skipZeros  
                                      -o Control_peaks_centre.matrix.gz

                                      #!/bin/bash
                                      #SBATCH --job-name=deeptools
                                      #SBATCH -N 1
                                      #SBATCH -n 8
                                      
plotHeatmap -m Control_peaks_centre.matrix.gz -out ATAC-histone-heatmap.png --colorMap Blues

done
exit 0