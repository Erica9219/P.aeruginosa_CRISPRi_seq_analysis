# CRISPRi-seq data analysis with DESeq2 
# zhangttingting9218@163.com 2023.11.09

# PRELIMINARIES
setwd("\\\\Kjzb\\张婷婷\\08.liu_lab\\10.zhangyu\\data for R\\FC_1")

Alpha <- 0.05 
FCthreshold <- 1

# Packages
library(DESeq2)
library(ggplot2)
library(writexl)
library(dplyr)
library(ggrepel)

Group <- "mymeta.csv"
head(Group)

set <- c("set1", "set2")
file <- c("Compiled-7g-high.csv", "Compiled-14g-high.csv", "Compiled-21g-high.csv", 
          "Compiled-7g-low.csv", "Compiled-14g-low.csv", "Compiled-21g-low.csv")

# Loop through files and sets
for(i in 1:length(file)){
  for(s in 1:length(set)){
    Set <- set[s]
    Data <- file[i]
    deseqout <- paste0(sub("\\.csv$", "", file[i]), ".", set[s],".",FCthreshold, ".res_raw.xlsx")
    volcano <- paste0(sub("\\.csv$", "", file[i]), ".", set[s],".",FCthreshold, ".volcano.pdf")
    
    # DATA PREPARATION
    raw <- read.csv(file[i], header = TRUE)
    
    names(raw) <- c("sgRNA",rep("control",(ncol(raw)-1)/2),rep("treatment",(ncol(raw)-1)/2) )
    
    raww <- raw[grepl(Set , raw$sgRNA), ]
    Non <- raw[grepl("NT" , raw$sgRNA), ]
    
    # Check if Non (NT) data exists; if not, skip to the next iteration
    if(nrow(Non) == 0) {
      next  # Skip this iteration if NT data is not available
    }
    
    raw <- rbind(raww, Non)
    cts <- as.matrix(raw[, -1])  # Exclude first column (Feature)
    rownames(cts) <- raw$sgRNA
    
    # Read coldata and create DESeq2 object
    coldata <- read.csv(Group)
    rownames(coldata) <- colnames(cts)
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~des)

    # Run DESeq2
    dds <- DESeq(dds)
    
    # RESULTS
    TMP <- resultsNames(dds)[2]
    res <- results(dds, name = TMP, alpha = Alpha, lfcThreshold = FCthreshold)
    raw_res <- data.frame(sgRNA = rownames(res), res)
    
    # Classify significant results
    raw_res %>%
      mutate(significant = case_when(
        log2FoldChange >= FCthreshold & padj <= Alpha ~ "UP",
        log2FoldChange <= -FCthreshold & padj <= Alpha ~ "DOWN",
        TRUE ~ "NOT_CHANGE"
      )) -> raw_res_1
    
    # Output results to Excel
    write_xlsx(raw_res_1, path = paste0(getwd(), "/", deseqout))
    
    # Volcano plot
    raw_res_11 <- na.omit(raw_res_1)
    p <- ggplot(raw_res_11, aes(log2FoldChange, -log10(padj), col = significant)) + 
      geom_point(size = 2) + 
      geom_vline(xintercept = c(-FCthreshold, FCthreshold), lty = 2, col = "blue") + 
      geom_hline(yintercept = -log10(Alpha), lty = 2, col = "red") + 
      scale_color_manual(values = c("DOWN" = "#56B4E9", "NOT_CHANGE" = "grey", "UP" = "#E69F00")) + 
      labs(x = "log2 Fold Change", y = expression(-log[10](padj))) + 
      theme_bw(base_size = 20)
    
    # Save volcano plot
    ggsave(volcano, p, path = paste0(getwd(), "/"), width = 10, height = 8)
  }
}
