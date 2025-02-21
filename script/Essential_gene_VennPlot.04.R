
# Process multiple 'res_raw.xlsx' files and classify the data based on log2FoldChange and padj
process_res_raw_files <- function(xlsx_files, n_files,FoldChange,Padj ) {
  # Get all files matching the pattern
  # Initialize an empty data frame to store results
  result <- NULL

  for (i in 1:n_files) {
    # Read the current Excel file and classify the data based on fold change and adjusted p-value
    set1 <- read_excel(xlsx_files[i]) %>%
      mutate(significant = case_when(
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
  
  # 找到 group == "UP" 且 log2FoldChange > 1 的基因
  genes_to_remove <- result %>%
    filter(log2FoldChange > 1 ) %>%
    pull(Gene) %>%
    unique()
  
  # 筛选出 group == "DOWN" 且 Gene 不在 genes_to_remove 中的行
  filtered_data <- result %>%
    filter(significant  == "DOWN" & !Gene %in% genes_to_remove)
 
  return(filtered_data)  # Return the final result data frame
}




library(dplyr) ##过滤包含特定字符串的行
library("ggvenn")
setwd("\\\\Kjzb\\张婷婷\\08.liu_lab\\10.zhangyu\\必需基因分类\\tmp")
# Get list of Excel files matching pattern
xlsx <- list.files(pattern = "res_raw.xlsx$")
xlsx

set1 <- xlsx[grep("set1", xlsx)]
set1

set2 <- xlsx[grep("set2", xlsx)]
set2

###############################set1 and set2 
set1_data <- process_res_raw_files(set1, 3 , 1 , 0.05)
head(set1_data)

set2_data <- process_res_raw_files(set2, 3 , 1 , 0.05)
head(set2_data)

length(unique(set1_data$Gene))


library(ggvenn) #加载ggvenn包
b <- list(`set1` = set1_data$Gene,
          `set2` = set2_data$Gene
          )

p <- ggvenn(b, show_elements = FALSE,fill_color = c('#BC80BD', '#FB8072','#80B1D3'),
            label_sep = "\n",stroke_size = 0.1,set_name_size = 7,
            text_size =5)

p
ggsave("set1_set2.pdf",p,width=6,height=6)



########################
###############################set1 and set2 

set1_set2_data <- unique(c(set2_data$Gene,set1_data$Gene))
length(set1_set2_data )
venn_dat  <- read.csv("veen_and_set1_set2.txt", check.names = FALSE,sep = "\t",header=T)


###set1
b <- list(`This study 618 genes` = set1_set2_data ,
          `Liberati et al. 2006 335 genes` = na.omit(venn_dat$tn2),
          `Skurnik et al. 2013 634 genes` = na.omit(venn_dat$Tn1),
          `Poulsen et al. 2019 437 genes` = na.omit(venn_dat$`2019`))

names(b) <- c(
  "This study\n618 genes", 
  "Liberati et al. 2006\n335 genes",
  "Skurnik et al. 2013\n634 genes",
  "Poulsen et al. 2019\n437 genes"
)


p <- ggvenn(b, show_elements = FALSE,fill_color = c('#BC80BD','#80B1D3', '#FB8072','#8DD3C7'),
            label_sep = "\n",stroke_size = 0.1,set_name_size = 4,
            text_size =4)

p
ggsave("set12_Tn12.成比例_2019.pdf",p,width=7,height=7)
