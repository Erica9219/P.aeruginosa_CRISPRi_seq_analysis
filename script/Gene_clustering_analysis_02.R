
# Set the working directory
set_working_directory <- function(directory_path) {
  setwd(directory_path)  # Change the working directory to the specified path
}

# Load necessary libraries
load_libraries <- function() {
  library(readxl)  # For reading Excel files
  library(writexl)  # For writing Excel files
  library(dplyr)  # For data manipulation
  library(tidyr)  # For data reshaping
  library(rio)  # For importing and exporting data
  library(ggplot2)  # For visualization
  library(reshape2)  # For reshaping data
}

# Process multiple 'res_raw.xlsx' files and classify the data based on log2FoldChange and padj
process_res_raw_files <- function(file_pattern, n_files,FoldChange,Padj ) {
  # Get all files matching the pattern
  xlsx_files <- list.files(pattern = file_pattern)
  xlsx_files 
    # Initialize an empty data frame to store results
  result <- NULL
  
  for (i in 1:n_files) {
    # Read the current Excel file and classify the data based on fold change and adjusted p-value
    set1 <- read_excel(xlsx_files[i]) %>%
      mutate(group = case_when(
        log2FoldChange >= FoldChange & padj <= Padj ~ "UP",  # Mark up-regulated genes
        log2FoldChange <= FoldChange & padj <=  Padj ~ "DOWN",  # Mark down-regulated genes
        TRUE ~ "NOT_CHANGE"  # Mark genes with no significant change
      ))
    
    result <- rbind(result, set1)  # Append the data to the result
  }
  
  # Separate sgRNA column into Gene and Set columns
  # Check if 'sgRNA' column contains '_set' before applying separate
  if (any(grepl("_set", result$sgRNA, fixed = TRUE))) {
    result <- result %>%
      separate(sgRNA, into = c("Gene", "Set"), sep = "_set", remove = FALSE)
  } else {
    print("No '_set' found in sgRNA column, renaming column to 'Gene'.")
    colnames(result)[colnames(result) == "sgRNA"] <- "Gene"  # Rename 'sgRNA' to 'Gene'
  }
  
  return(result)  # Return the final result data frame
 }



# Merge data from multiple files
merge_files_data <- function(file_pattern, n_files,FoldChange,Padj) {
  xlsx_files <- list.files(pattern = file_pattern)
  xlsx_files 
  # Create a data frame with unique sgRNAs
  gene_id <- as.data.frame(unique(process_res_raw_files(file_pattern, n_files,FoldChange,Padj)$Gene))
  head(gene_id)
  names(gene_id) <- "Gene"
  
#  tmp <- gene_id  # Initialize with unique sgRNAs
  
  for (i in 1:n_files) {
    set <- read_excel(xlsx_files[i])
    set <- if (any(grepl("_set", set$sgRNA))) {
      separate(set, sgRNA, into = c("Gene", "Set"), sep = "_set", remove = FALSE)
    } else {
      rename(set, Gene = sgRNA)
    }
    set <- set %>%
      select(Gene, log2FoldChange) %>%
      rename_with(~ paste0(strsplit(xlsx_files[i], ".res_raw")[[1]][1], ".log2FoldChange"), .cols = 2)
    gene_id <- merge(gene_id, set, by = "Gene", all.x = TRUE)
  }
  gene_id # Return the merged data
}


# Merge essential gene data and handle missing values
merge_essential_data <- function(tmp_data, essential_file) {
  essential <- read_excel(essential_file)  # Read essential gene file
  
  # Merge the essential data and add a type column for essential or non-essential genes
  result <- merge(tmp_data, essential, by.x = "Gene", by.y = "Gene", all.x = TRUE) %>%
    mutate(type = ifelse(is.na(type), "non-essential", type)) 
  
  return(result)  # Return the result
}

# Merge essential gene data and handle missing values
merge_essential_data <- function(tmp_data, essential_file) {
  essential <- read_excel(essential_file)  # Read essential gene file
  
  # Merge the essential data and add a type column for essential or non-essential genes
  result <- merge(tmp_data, essential, by.x = "Gene", by.y = "Gene", all.x = TRUE) %>%
    mutate(type = ifelse(is.na(type), "non-essential", type))
 #   select(Gene, type, starts_with("log2FoldChange"))  # Select relevant columns
  result <- result[,c(1,6,2,4,7,3,5,8)]
  
  return(result)  # Return the result
}

# # # Example usage
#  essential_file = "This_study_essential.xlsx"
#  tmp_data <- gene_id  # Assuming `gene_id` is already defined
# # 
# # # Call the function and store the result in 'merged'
# merged <- merge_essential_data(tmp_data, essential_file)
# # 
# # # View the result
# head(merged)  # Use the 'merged' variable to inspect the result
# # 
# # 


head(M)
# Perform clustering based on gene expression patterns
cluster_data <- function(M) {
  for(i in 1:2){
    MT <- M[M$type == "essential", ]
    a=(i-1)*3+2
    b=i*3+1
    M1 <- MT[,c(a:b)]
    head(MT)
    dim(M1)
    names(M1) <- c("g7","g14","g21")
    rownames(M1)
    M1$name <- MT$Gene
    head(M1)
    M2 <- M1
    head(M2)
    #######Very essential 
    C1  <- subset(M2,M2$g7 < -4  |M2$g14< -4| M2$g21< -4)
    dim(C1)
    head(C1)
    ###########Very essential and quick response
    C11 <- subset(C1,C1$g7< -2 )
    head(C11)
    
    C11$cluster <- "C1"
    head(C11)
    dim(C11) 
    ###########Very essential and slow response
    C12 <- subset(C1,C1$g7 >= -2 )
    dim(C12)
    head(C12)
    C12$cluster <- "C2"
    
    ####essential
    M2 <- M2[!rownames(M2) %in% rownames(C1) , ]
    dim(M2)
    head(M2,20)
    ######
    C2 <- subset(M2,M2$g7 < 1 & M2$g14 < 1 & M2$g21 < 1)
    head(C2)
    C2 <-  subset(C2,C2$g7 < -1 | C2$g14 < -1 |C2$g21 < -1)
    head(C2)
    dim(C2)
    ########Essential and quick response
    C21 <- subset(C2,C2$g7< -1 )
    dim(C21)
    C21$cluster <- "C3"
    #####Essential and slow response
    C22 <- subset(C2,C2$g7 >= -1 )
    dim(C22)
    C22$cluster <- "C4"
    head(C22)
    ###
    
    
    M1_1  <-  M[!rownames(M) %in% C1$name , ]
    head(M1_1)
    M1_2  <-  M1_1[!rownames(M1_1 ) %in% C2$name , ]
    head(M1_2)
    C5 <- M1_2[,c(a:b)]
    head(C5)
    names(C5) <- c("g7","g14","g21")
    C5$name <- M1_2$Gene
    tail(C5)
    
    
    C5$cluster <- "C5"
    dim(C5)
    tail(C5)
    
    CC <- rbind(C11,C12,C21,C22,C5)
    head(CC)
    CC1 <- as.data.frame(CC)
    head(CC1)
    write_xlsx(CC1,paste0(i,".newCluster.xlsx"))
    #write.table(CC1,paste0(i,".newCluster.txt"),sep="\t", row.names = F, col.names = T,quote = FALSE)
  }  
}



# Select the best cluster from the results
select_best_cluster <- function(input_file = "1.newCluster.xlsx",input_file1 = "2.newCluster.xlsx", output_file = "Select_the_best_cluster_in_different_set.xlsx") {
  # Read the data from two sheets in the Excel file
  set1 <- read_excel(input_file, sheet = 1)
  set2 <- read_excel(input_file1, sheet = 1)
  
  set1$set <- "set1"  # Label the first set
  set2$set <- "set2"  # Label the second set
  
  tmp <- rbind(set1, set2)  # Combine both sets
  head(tmp)
  Result <- NULL  # Initialize an empty data frame for the result
  unique_names <- unique(tmp$name)  # Get unique gene names
  
  for (b in 1:length(unique_names)) {
    tt <- tmp[tmp$name == unique_names[b], ]  # Filter by gene name
    if (nrow(tt) == 2) {
      # Select the best cluster based on the cluster number
      if (as.numeric(str_split(tt[1, 5], "C")[[1]][2]) >= as.numeric(str_split(tt[2, 5], "C")[[1]][2])) {
        Result <- rbind(Result, tt[2, ])
      } else {
        Result <- rbind(Result, tt[1, ])
      }
    } else {
      Result <- rbind(Result, tt)  # If only one result, include it directly
    }
  }
  
  # Write the selected results to a text file
  write_xlsx(Result, output_file)
}
#######################################################################################
# Main program entry
main <- function() {
  Padj <- 0.05
  FoldChange <- 1
  # Set working directory and load libraries
  set_working_directory("\\\\Kjzb\\张婷婷\\08.liu_lab\\10.zhangyu\\必需基因分类\\tmp")
  load_libraries()
  
  # Process data from 'res_raw.xlsx' files
  result <- process_res_raw_files("res_raw.xlsx$", 6,FoldChange,Padj )
  essential <- data.frame(Gene= result[result$significant=="DOWN",]$Gene, type="essential")
  head( essential)
  write_xlsx(essential,"PA14_all_essential_Genes.xlsx")
  
  # Merge data from multiple files
  tmp_data <- merge_files_data("res_raw.xlsx$", 6,FoldChange,Padj )
  
  # Merge essential gene data
  merged_data <- merge_essential_data(tmp_data,"PA14_all_essential_Genes.xlsx")
  merged_data <- unique(merged_data)
  # Save merged data to an Excel file
  write_xlsx(merged_data, "PA14_log2FC_data.xlsx")
  
  # Perform clustering on the data
  M <- read_xlsx("PA14_log2FC_data.xlsx", sheet = 1) %>% as.data.frame()
  head(M)
  rownames(M) <- M$Gene
  table(M$Gene)
  cluster_data(M)
  # Select the best cluster and output the results
  select_best_cluster()
}

# Run the main program
main()

