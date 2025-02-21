# CRISPRi-seq for Gallium Synergy Target Identification in Pseudomonas aeruginosa
Author: Tingting Zhang [zhangtingting9218@163.com]

This repository contains an R-based software package developed for the analysis of CRISPR interference (CRISPRi) and genome-wide screening sequencing data in Pseudomonas aeruginosa. The method enables the quantitative measurement of gene essentiality and the identification of genes that interact with non-antibiotic therapies, such as gallium.

__sgRNA Design__:
The sgRNA design for CRISPR interference (CRISPRi) in this repository follows the pipeline outlined in the CRISPRi-seq repository (https://github.com/vincentdebakker/CRISPRi-seq). This pipeline automates the design of highly specific and efficient single-guide RNAs (sgRNAs) for gene knockdown in bacterial genomes.

__Scripts__:
The scripts used to generate individual figure panels, as well as for data processing and analysis.

__--Deseq2_analysis_from_2fastq2_result_01.R__ Used the DESeq2 package for differential gene expression analysis between different experimental conditions.

__--Gene_clustering_analysis_02.R__  Designed for identifying and classifying essential genes based on CRISPRi-seq data.

__--Gene_Essentiality_Visualization_lineplot_03.R__ Comprehensive data processing and analysis for clustering results based on log2 Fold Change values across different gene clusters. It integrates multiple data sources, including gene expression data and COG (Cluster of Orthologous Groups) annotations, through a series of preprocessing steps. it visualizes the results by generating detailed line plots, which illustrate the trends of log2 Fold Change values over generations for each cluster, and bar plots that highlight the distribution of COG categories within each gene cluster. 

__--Essential_gene_VennPlot.04.R__ Generates Venn diagrams to compare gene sets between the current study and previous studies.

__--gene_classification_circos_05.R__ Uses the circlize package to plot genomic regions, with custom color coding for genes from various datasets.


