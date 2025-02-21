setwd("\\\\Kjzb\\张婷婷\\08.liu_lab\\10.zhangyu\\必需基因分类\\tmp")
library(xlsx)
library(readxl)
library(writexl)
library(plyr)
library(reshape2)
library(ggplot2)
library(ggpmisc)
library(dplyr)

Result <- read_excel("Select_the_best_cluster_in_different_set.xlsx")
Result[Result $name=="PA14_03410",]
# 
COG <- read.table("Pseudomonas_aeruginosa_UCBPP-PA14_109.cog.txt", sep = "\t", header = T,quote="\"",check.names=FALSE)
cc <- NULL
for(b in 1:nrow(COG)){
  if ( length(str_split(COG[b,3],"")[[1]])==1){
    cc <- rbind(cc,COG[b,])
    
  }else{
    for (a in 1:length(str_split(COG[b,3],"")[[1]])){
      cc1 <-  c(COG[b,1],COG[b,2],str_split(COG[b,3],"")[[1]][a])
      cc <- rbind(cc,cc1)
    }
  }
}
head(cc,10)  


m1 <- merge(Result,cc,by.x="name",by.y="query",all=T)

# 
write_xlsx(m1, path = paste0("tmp",1,"_",2,"cog.xlsx"))
# 

m1 <- melt( m1)
m1  <-m1  %>%
   mutate(variable = gsub("g", "", variable))  # 去除 'g'
head(m1)
write_xlsx(m1, path = paste0("tmp",1,"_",2,"cog.xlsx"))
#

Result <- m1

Result <- Result[!grepl("NT", Result$name), ]

plot_line_and_bar <- function(cluster_num, color) {
  # cluster_num = 1
  # color = "#F09148"

  name1 <- paste0("C", cluster_num)
  name1
  C <- unique(Result[Result$cluster == name1, c(1, 2, 6,7)])
  name2 <- paste0("Cluster", cluster_num, "_", length(unique(na.omit(C$name))))
 

  name2 
  C_TMP <- data.frame(name = C$name, cluster = rep("C1", nrow(C)), value = rep(0, nrow(C)), variable  = rep(0, nrow(C)))
  
  C_TMP <- unique(C_TMP)
  c12 <- rbind(C, C_TMP)
  head(c12)
  p <- ggplot(c12, aes(x = as.numeric(variable), y = value)) +
    geom_line(aes(group = name), color = color, alpha = 0.5) +
    geom_point() +
    geom_smooth(method = lm, se = TRUE, color = "black", fill = "#8B8989") +
    stat_poly_eq(formula = y ~ x, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE) +
    theme_bw() +
    labs(title = name2, x = "Generation", y = "log2 Fold Change") +
    scale_y_continuous(limits = c(-13, 4)) +
    scale_x_continuous(limits = c(0, 21), breaks = seq(0, 21, 7)) +
    theme(axis.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          panel.grid.major.x = element_blank())
  p

  C <- unique(Result[Result$cluster == name1, ][, c(1, 2, 5)])
  F11 <- as.data.frame(table(C$COG_category))
  F1 <- head(F11[order(F11$Freq, decreasing = TRUE), ], 5)
  F1 <- F1[order(F1$Freq, decreasing = FALSE), ]
  
  p1 <- ggplot(F1, aes(x = factor(Var1, levels =Var1 ), y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", width = 0.5) +
    coord_flip() +
    theme_bw() +
    theme(axis.title = element_text(size = 13), 
          axis.text = element_text(size = 10, color = "black"),
          plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
          legend.position = "none") +
    labs(y = "count", x = "")
  p1
  g <- ggplotGrob(p1)
  p2 <- p + annotation_custom(g, xmin = -0.7, xmax = 6.5, ymin = -14, ymax = -5)
  
  p2
 # ggsave(paste0(name1, "tmp.jpg"), plot = p2, width = 6, height = 6)
  ggsave(paste0(name1, "tmp.pdf"), plot = p2, width = 6, height = 6)
}

# 主循环
for (c in 1:4) {
  plot_line_and_bar(c, "#F09148")
}
plot_line_and_bar(5, "#427AB2")


