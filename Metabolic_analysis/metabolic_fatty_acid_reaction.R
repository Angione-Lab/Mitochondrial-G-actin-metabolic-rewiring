library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggplotify) 
library(heatmaply)

df = read_excel('Results/high_glucose_results.xls',sheet = 'reactions')
df <- df[df$Pathway %in% c('Fatty acid oxidation', 'Fatty acid biosynthesis','Fatty acid biosynthesis (even-chain)'),]

df1 <- df[df$`log2(Foldchange m4_vs_wildtype)` < -0.001, ]
df2 <- df[df$`log2(Foldchange m4_vs_wildtype)` > 0.001, ]
df <- data.frame(rbind(df1, df2))

df <- df[!duplicated(df$Reaction.name),]

df <- data.frame(df, check.names = FALSE)

df <- df[order(df$Pathway), ]


row.names(df) = df$Reaction.name
pathways <- data.frame(pathways = df$Pathway)
rownames(pathways) <- rownames(df)

df <- df %>% select(-c('Reaction.ID','Reaction.name','Pathway'))
df <- df %>% select(c("APC_WT1" ,"APC_WT2","APC_WT3","APC_m4_1","APC_m4_2", "APC_m4_3" ))

df <- data.frame(t(df), check.names = FALSE)
df$Condition <- rownames(df)
# Convert id and time into factor variables

df$Condition<-gsub("_1","",as.character(df$Condition))
df$Condition<-gsub("_2","",as.character(df$Condition))
df$Condition<-gsub("_3","",as.character(df$Condition))

df$Condition<-gsub("1","",as.character(df$Condition))
df$Condition<-gsub("2","",as.character(df$Condition))
df$Condition<-gsub("3","",as.character(df$Condition))

annoCol <- colorRampPalette(rev(c("#4BDC98","#FDAE61" , "#667761","#DB6A4B", "#E0F3F8" ,"#DBBC4B","#4575B4")))(length(unique(pathways$pathways)))
names(annoCol) <- unique(pathways$pathways)
annoCol <- list(pathways = annoCol)

ann_colors = list(
  Condition = c(APC_WT='#5DADE2', APC_m4='#D6EAF8'))

group <- df %>% dplyr::select('Condition')

df <- df %>% select(-c(Condition))

df <- t(normalize(df))
p = pheatmap(df,
             annotation_col = group,
             cluster_row = FALSE, 
             #scale = "row",
             cluster_cols = FALSE,
             #annotation_colors = annoCol, 
             #angle_col = 315,
             annotation_row = pathways, 
             annotation_colors = ann_colors,
             color=colorRampPalette(c("navy", "white", "red"))(10000)
             )
print(p)

ggsave(
  paste("fatty_acid_reactions", ".pdf", sep = ""),
  plot = p,
  device = NULL,
  path = 'Figures/',
  scale = 1,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  width = 15, height = 15,
  bg = NULL
)
