library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggplotify) 
library(heatmaply)


setwd("C:/Users/u0034546/OneDrive - Teesside University/1.PhDworks/Maria_project/MetabolicModelling/Simulate_FVA")

df = read_excel('C:/Users/u0034546/OneDrive - Teesside University/1.PhDworks/Maria_project/MetabolicModelling/Simulate_FVA/Results/high_glucose_results.xls',sheet = 'reactions')
df <- df %>% select(c('Pathway',"Foldchange m4_vs_wildtype"))

pathways_df <- df %>%
  group_by(Pathway) %>%
  summarize(Foldchange = mean(`Foldchange m4_vs_wildtype`))
pathways_df <- data.frame(pathways_df)


row.names(pathways_df) <- pathways_df$Pathway
pathways_df <- pathways_df %>% select(-c('Pathway'))

pathways_df <- pathways_df %>%
  mutate(
    LogFoldChange = log2(`Foldchange`)
  )

pathways_df <- pathways_df %>% select(-c('Foldchange'))


upregulated <- pathways_df %>%
  filter(
    LogFoldChange > quantile(LogFoldChange, 0.80) & LogFoldChange > 0.001)

downregulated <- pathways_df %>%
  filter(
    LogFoldChange < quantile(LogFoldChange, 0.15) & LogFoldChange < 0)

logfc_df <- rbind(upregulated, downregulated)

breaks = unique(c(seq(min(logfc_df), 0, length.out = 500), seq(0, max(logfc_df), length.out = 500)))

p = pheatmap(logfc_df,
             #annotation_col = group,
             cluster_row = FALSE, 
             #scale = "column",
             cluster_cols = FALSE,
             angle_col = 0,
             annotation_colors = ann_colors,
             border_color = NA,
             breaks = breaks,
             color=colorRampPalette(c("navy", "white", "red"))(1000))
print(p)
ggsave(
  paste("pathways", ".pdf", sep = ""),
  plot = p,
  device = NULL,
  path = 'Figures/',
  scale = 1,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  width = 5, height = 5,
  bg = NULL
)


