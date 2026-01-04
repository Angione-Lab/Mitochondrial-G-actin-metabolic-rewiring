library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggplotify) 
library(heatmaply)

setwd("C:/Users/u0034546/OneDrive - Teesside University/1.PhDworks/Maria_project/MetabolicModelling/Simulate_FVA")

df = read_excel('C:/Users/u0034546/OneDrive - Teesside University/1.PhDworks/Maria_project/MetabolicModelling/Simulate_FVA/Results/high_glucose_results.xls',sheet = 'reactions')
df <- df[df$Pathway == "Glycolysis / Gluconeogenesis",]
df <- df[df$`log2(Foldchange m4_vs_wildtype)` !=0, ]
#df <- df[!duplicated(df$`Reaction name`),]


glycolysis_reactions <- c(
  'MAR04394', 'MAR04381', 'MAR04379', 'MAR04375', 'MAR04391', 'MAR04373', 'MAR04368', 
  'MAR04371', 'MAR04363', 'MAR04360', 'MAR04358', 'MAR04355', 'MAR04365', 'MAR04372', 
  'MAR04370', 'MAR04301', 'MAR07747', 'MAR04388', 'MAR04283', 'MAR08357', 'MAR04097', 
  'MAR04099', 'MAR04108', 'MAR04133', 'MAR04137', 'MAR06410', 'MAR06412'
)

df <- df %>%
  filter(`Reaction ID` %in% glycolysis_reactions)

df$`Reaction ID` <- ordered(df$`Reaction ID`,levels = glycolysis_reactions)
df <- df[order(df$`Reaction ID`),]
df <- data.frame(df, check.names = FALSE)

row.names(df) = df$`Reaction name`
df <- df %>% select(-c('Reaction.ID','Reaction name','Pathway'))


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
  paste("glycolysis_reactions", ".pdf", sep = ""),
  plot = p,
  device = NULL,
  path = 'Figures/',
  scale = 1,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  width = 15, height = 7,
  bg = NULL
)
