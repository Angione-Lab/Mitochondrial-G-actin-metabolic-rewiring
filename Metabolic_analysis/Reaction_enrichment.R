library(readxl)
library(dplyr)  
library(pheatmap)
library(RColorBrewer)

foldchangeRxn = read_excel('../GSMM/Simulate_FVA/Results_fpkm/reaction_foldchange.xls',sheet = 'reactions')
foldchangeRxn <-foldchangeRxn[order(foldchangeRxn$pathways, decreasing=FALSE), ]
foldchangeRxn <- foldchangeRxn[!duplicated(foldchangeRxn$`Reaction Names`), ]
foldchangeRxn <- foldchangeRxn[foldchangeRxn$pathways != "Transport reactions",]

df1 <- foldchangeRxn[foldchangeRxn$`FC M4 vs wildtype` < -0.01, ]
df2 <- foldchangeRxn[foldchangeRxn$`FC M4 vs wildtype` > 0.01, ]

foldchangeRxn <- data.frame(rbind(df1, df2))


foldchange_df <- foldchangeRxn %>% select(c("FC.M4.vs.wildtype"))
rownames(foldchange_df) <- foldchangeRxn$`Reaction Names`



pathways <- data.frame(pathways = foldchangeRxn$pathways)
rownames(pathways) <- foldchangeRxn$`Reaction Names` 

# creat colours for each group
#newCols <- colorRampPalette(grDevices::rainbow(length(unique(pathways$pathways))))

annoCol <- colorRampPalette(rev(c("#4BDC98","#FDAE61" , "#667761","#DB6A4B", "#E0F3F8" ,"#DBBC4B","#4575B4")))(length(unique(pathways$pathways)))
#annoCol <- newCols(length(unique(pathways$pathways)))
names(annoCol) <- unique(pathways$pathways)
annoCol <- list(pathways = annoCol)

paletteLength = length(foldchangeRxn$Reactions)

myColor1 <- colorRampPalette(c('blue', 'white' , 'red'))(paletteLength)

myBreaks <- c(seq(min(foldchange_df), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(foldchange_df)/paletteLength, max(foldchange_df), length.out=floor(paletteLength/2)))



# add row annotations
#pheatmap(foldchange_df, cluster_cols = F, cluster_rows = F, annotation_row = pathways)

xx<- pheatmap(foldchange_df, 
         cluster_cols = F, 
         cluster_rows = F, 
         annotation_row = pathways, 
         color = myColor1,
         breaks = myBreaks,
         annotation_colors = annoCol, 
         cellheight = 10, 
         cellwidth = 30,
         
         ) 
save_pheatmap_pdf <- function(x, filename, width=20, height=40) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(xx, "Reaction_Enrichment_exp1.pdf")
