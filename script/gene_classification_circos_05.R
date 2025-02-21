library(dplyr)       # Data manipulation
library(ggvenn)     # Venn diagram visualization
library(readxl)     # Reading Excel files
library(circlize)   # Circular visualization
library(grid)       # Combining circular plots with legends
library(purrr)
# Function to process multiple 'res_raw.xlsx' files
process_res_raw_files <- function(xlsx_files, n_files, FoldChange, Padj) {
  result <- NULL  # Initialize an empty data frame
  
  for (i in 1:n_files) {
    # Read Excel file and classify genes based on fold change and p-adjusted values
    set1 <- read_excel(xlsx_files[i]) %>%
      mutate(significant = case_when(
        log2FoldChange >= FoldChange & padj <= Padj ~ "UP",    # Upregulated genes
        log2FoldChange <= -FoldChange & padj <= Padj ~ "DOWN",  # Downregulated genes
        TRUE ~ "NOT_CHANGE"                                      # No significant change
      ))
    
    result <- rbind(result, set1)  # Append data
  }
  
  # Separate 'sgRNA' column into 'Gene' and 'Set' if applicable
  if (any(grepl("_set", result$sgRNA, fixed = TRUE))) {
    result <- result %>% separate(sgRNA, into = c("Gene", "Set"), sep = "_set", remove = FALSE)
  } else {
    colnames(result)[colnames(result) == "sgRNA"] <- "Gene"
  }
  
  # Identify genes with log2FoldChange > 1 in "UP" group
  genes_to_remove <- result %>%
    filter(log2FoldChange > 1) %>%
    pull(Gene) %>%
    unique()
  
  # Keep only "DOWN" genes that are not in genes_to_remove
  filtered_data <- result %>%
    filter(significant == "DOWN" & !Gene %in% genes_to_remove)
  
  return(filtered_data)
}

# Set working directory
setwd("\\\\kjzb\\张婷婷\\08.liu_lab\\10.zhangyu\\必需基因分类\\tmp")

# Get list of Excel files
xlsx <- list.files(pattern = "res_raw.xlsx$")
set1 <- xlsx[grep("set1", xlsx)]
set2 <- xlsx[grep("set2", xlsx)]

# Process set1 and set2 files
set1_data <- process_res_raw_files(set1, 3, 1, 0.05)
set2_data <- process_res_raw_files(set2, 3, 1, 0.05)

# Combine unique genes from set1 and set2
set1_set2_data <- unique(c(set2_data$Gene, set1_data$Gene))
set <- data.frame(Gene = set1_set2_data)
head(set)
# Load gene start-end positions
gene <- read.csv("gene.start.end.txt", sep = "\t", col.names = c("Gene", "start", "end"))

# Load genome file
genome <- read.csv("tmp1.txt", sep = "\t", header = TRUE)

# Load Venn diagram data
venn_dat <- read.csv("veen_and_set1_set2.txt", sep = "\t", header = TRUE)
names(venn_dat) <- c("2019", "2006", "2013")

# Merge data for circular plot
create_circos_data <- function(venn_dat,venn_col, gene_data, type_value) {
  merged_data <- merge(venn_dat, gene_data, by.x = venn_col, by.y = "Gene")
  merged_data <- merged_data[, c(ncol(merged_data) - 1, ncol(merged_data))]
  merged_data$seq_id <- "genome"
  merged_data$type <- type_value
  merged_data <- merged_data[,c(3,1,2,4)]
  colnames(merged_data) <- c("seq_id", "start", "end", "type")
  return(merged_data)
}

M1 <- create_circos_data(set, "Gene", gene, 1)

M2 <- create_circos_data(venn_dat,"2006", gene, 2)
M3 <- create_circos_data(venn_dat,"2013", gene, 3)
M4 <- create_circos_data(venn_dat,"2019", gene, 4)

# Initialize circular plot
circos.genomicInitialize(genome, plotType = c('axis', 'labels'), major.by = 1000000, track.height = 0.04)

# Define color mapping
color_assign <- colorRamp2(breaks = c(1, 2, 3, 4), color = c('#BC80BD', '#80B1D3', '#FB8072', '#8DD3C7'))

# Function to add track to circular plot
add_circos_track <- function(data) {
  circos.genomicTrackPlotRegion(
    data, track.height = 0.1, stack = TRUE, bg.border = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
    }
  )
}

# Add tracks
# add_circos_track(M4)
# add_circos_track(M3)
# add_circos_track(M2)
# add_circos_track(M1)


walk(list(M4, M3, M2, M1), add_circos_track)
