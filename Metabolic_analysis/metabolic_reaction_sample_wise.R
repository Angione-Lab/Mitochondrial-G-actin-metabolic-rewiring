library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggplotify) 
library(heatmaply)

setwd("C:/Users/u0034546/OneDrive - Teesside University/1.PhDworks/Maria_project/MetabolicModelling/Simulate_FVA")

df = read_excel('C:/Users/u0034546/OneDrive - Teesside University/1.PhDworks/Maria_project/MetabolicModelling/Simulate_FVA/Results/high_glucose_results.xls',sheet = 'reactions')

fasn_rxns = read.csv("Results/FASN_reactions.csv")

df <- df[df$`log2(Foldchange m4_vs_wildtype)` !=0, ]

df <- as.data.frame(df)
df <- df[df$`Reaction ID` %in% fasn_rxns$FASN_rxns1, ]
df <- df[!duplicated(df$`Reaction name`),]

#df$`Reaction ID` <- ordered(df$`Reaction ID`,levels = glycolysis_reactions)
df <- df[order(df$`Reaction name`),]
df <- data.frame(df, check.names = FALSE)

row.names(df) = df$`Reaction name`
df <- df %>% select(-c('Reaction ID','Reaction name','Pathway'))


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
             #angle_col = 315,
             annotation_colors = ann_colors,
             color=colorRampPalette(c("navy", "white", "red"))(10000))
print(p)

ggsave(
  paste("FASN_reactions", ".pdf", sep = ""),
  plot = p,
  device = NULL,
  path = 'Figures/',
  scale = 1,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  width = 15, height = 15,
  bg = NULL
)
