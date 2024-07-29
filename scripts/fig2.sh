# to make the pielup vizualisations we need the bw files for chip and atac data plus the peak files.
K4=Comprehensive_analysis_wt_vs_4dpi_schMedS3h1/swap/005_MACS2/wt_H3K4me3_SPMR/wt_H3K4me3_SPMR_treat_pileup.bw
K27=Comprehensive_analysis_wt_vs_4dpi_schMedS3h1/swap/005_MACS2/wt_H3K27ac_SPMR/wt_H3K27ac_SPMR_treat_pileup.bw 
ATAC=Comprehensive_analysis_wt_vs_4dpi_schMedS3h1/swap/002_merged_bam/merged_SM/SM_SPMR/SM_SPMR_treat_pileup.bw

s1190-schmidtea:/projects/mivanko/Bioinfo/Comparative_Genomics_Schmidtea/Comprehensive_analysis_wt_vs_4dpi_schMedS3h1/swap/005_MACS2/wt_H3K27ac_SPMR/wt_H3K27ac_SPMR_treat_pileup.bw
s1190-schmidtea:/projects/mivanko/Bioinfo/Comparative_Genomics_Schmidtea/Comprehensive_analysis_wt_vs_4dpi_schMedS3h1/swap/002_merged_bam/merged_SM/SM_SPMR/SM_SPMR_treat_pileup.bw  
# supporting_information/BEDlike_files_comparative_genomics/Smed/SMED_accessible_chromatin_wt_H3K27ac_withH3K4me3_summits.bed
# supporting_information/BEDlike_files_comparative_genomics/Smed/SMED_accessible_chromatin_wt_H3K27ac_noH3K4me3_summits.bed

K4=wt_H3K4me3_SPMR_treat_pileup.bw
K27=wt_H3K27ac_SPMR_treat_pileup.bw 
ATAC=SM_SPMR_treat_pileup.bw


declare -A associative
associative=([supporting_information/BEDlike_files_comparative_genomics/Smed/SMED_accessible_chromatin_wt_H3K27ac_noH3K4me3_summits.bed]=ATAC_H3K27ac_noH3K4me3 \
[supporting_information/BEDlike_files_comparative_genomics/Smed/SMED_accessible_chromatin_wt_H3K27ac_withH3K4me3_summits.bed]=ATAC_H3K27ac_withH3K4me3 \
[supporting_information/BEDlike_files_comparative_genomics/Smed/SMED_accessible_chromatin.narrowPeak]=ATAC_all)

for i in ${!associative[@]}; do \
computeMatrix reference-point \
--referencePoint TSS \
-b 1000 -a 1000 \
--scoreFileName ${K4} ${K27} ${ATAC} \
--regionsFileName $i \
--outFileName ${associative[$i]}_summit_wt_SPMR_computeMatrix.mat.gz \
--outFileNameMatrix ${associative[$i]}_summit_wt_SPMR_computeMatrix_Matrix \
--outFileSortedRegions ${associative[$i]}_summit_wt_SPMR_computeMatrix_SortedRegions \
--samplesLabel wt_H3K4me3 wt_H3K27ac wt_ATACseq \
--numberOfProcessors max \
&& \
plotHeatmap \
--matrixFile ${associative[$i]}_summit_wt_SPMR_computeMatrix.mat.gz \
--outFileName ${associative[$i]}_summit_wt_SPMR_plotHeatmap.svg \
--averageType mean \
--plotType lines \
--legendLocation upper-right \
--whatToShow 'plot, heatmap and colorbar' \
--zMin 0 0 \
--zMax 5 5 \
--yMin 0 0 \
--yMax 5 5 \
--legendLocation upper-right \
--colorList '#ffffff,#cc79a7' '#ffffff,#0072b2' '#ffffff,#e69f00' \
;done