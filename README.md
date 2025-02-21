# CRISPRi-seq for Gallium Synergy Target Identification in Pseudomonas aeruginosa
Author: Tingting Zhang [zhangtingting9218@163.com]

This repository contains an R-based software package developed for the analysis of CRISPR interference (CRISPRi) and genome-wide screening sequencing data in Pseudomonas aeruginosa. The method enables the quantitative measurement of gene essentiality and the identification of genes that interact with non-antibiotic therapies, such as gallium.

sgRNA Design:
The sgRNA design for CRISPR interference (CRISPRi) in this repository follows the pipeline outlined in the CRISPRi-seq repository (https://github.com/vincentdebakker/CRISPRi-seq). This pipeline automates the design of highly specific and efficient single-guide RNAs (sgRNAs) for gene knockdown in bacterial genomes.

Scripts:
The scripts used to generate individual figure panels, as well as for data processing and analysis.

--Deseq2_analysis_from_2fastq2_result_01.R   This R script is used for CRISPRi-seq data analysis, primarily leveraging the DESeq2 package for differential gene expression analysis between different experimental conditions. The script also generates volcano plots to visualize the results.

--Gene_clustering_analysis_02.R This R script is designed for identifying and classifying essential genes based on CRISPRi-seq data.

--Gene_Essentiality_Visualization_lineplot_03.R This R script conducts comprehensive data processing and analysis for clustering results based on log2 Fold Change values across different gene clusters. It integrates multiple data sources, including gene expression data and COG (Cluster of Orthologous Groups) annotations, through a series of preprocessing steps. it visualizes the results by generating detailed line plots, which illustrate the trends of log2 Fold Change values over generations for each cluster, and bar plots that highlight the distribution of COG categories within each gene cluster. 

--Essential_gene_VennPlot.04.R The script then generates Venn diagrams to compare gene sets between the current study and previous studies.

--gene_classification_circos_05.R The script uses the circlize package to plot genomic regions, with custom color coding for genes from various datasets.
