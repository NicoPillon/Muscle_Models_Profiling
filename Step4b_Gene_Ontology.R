#=============================================================================================================================
# Gene Ontology using ClusterProfiler https://github.com/GuangchuangYu/clusterProfiler/issues/32
#=============================================================================================================================
setwd("C:/ownCloud/Projects/MuscleModels/Data/Transcriptomics")
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
keytypes(org.Hs.eg.db) # use keytypes to list all supporting types

#genes have to be ENTREZID!
Annotation <- read.delim("C:/ownCloud/R/Annotation_All.txt", sep="") # This is my homemade annotation file
MouseC2C12 <- read.table("../../Stats/MouseC2C12_Specific.txt")
MouseC2C12 <- merge(Annotation, MouseC2C12, by.x=3, by.y=0, all=F)
MouseC2C12 <- data.frame(aggregate(MouseC2C12[,6:ncol(MouseC2C12)],by = list(MouseC2C12$ENTREZID), FUN=mean, na.rm=TRUE), row.names = 1)
RatL6 <- read.table("../../Stats/RatL6_Specific.txt")
RatL6 <- merge(Annotation, RatL6, by.x=3, by.y=0, all=F)
RatL6 <- data.frame(aggregate(RatL6[,6:ncol(RatL6)],by = list(RatL6$ENTREZID), FUN=mean, na.rm=TRUE), row.names = 1)
HumanCell <- read.table("../../Stats/HumanCell_Specific.txt")
HumanCell <- merge(Annotation, HumanCell, by.x=3, by.y=0, all=F)
HumanCell <- data.frame(aggregate(HumanCell[,6:ncol(HumanCell)],by = list(HumanCell$ENTREZID), FUN=mean, na.rm=TRUE), row.names = 1)

#Make matrix
HumanCell$entrez <- rownames(HumanCell)
HumanCell$group <- "Human Primary"
HumanCell$direction <- "up"
HumanCell$direction[HumanCell$LogFC_mean <0] <- "down"
MouseC2C12$entrez <- rownames(MouseC2C12)
MouseC2C12$group <- "Mouse C2C12"
MouseC2C12$direction <- "up"
MouseC2C12$direction[MouseC2C12$LogFC_mean <0] <- "down"
RatL6$entrez <- rownames(RatL6)
RatL6$group <- "Rat L6"
RatL6$direction <- "up"
RatL6$direction[RatL6$LogFC_mean <0] <- "down"
mydf <- rbind(HumanCell[,6:9], MouseC2C12[,6:9], RatL6[,6:9])
rm(Annotation, HumanCell, MouseC2C12, RatL6)

#Gene Ontology
bp <- compareCluster(entrez~group+direction, data=mydf,
                              fun="enrichGO", ont = "BP",
                              OrgDb = org.Hs.eg.db, pvalueCutoff=0.05)
saveRDS(bp, "../Stats/enrichGO.Rds")


#=============================================================================================================================
#Gene Ontology plot and simplify
#=============================================================================================================================
setwd("C:/ownCloud/Projects/MuscleModels/Data/Transcriptomics")
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
keytypes(org.Hs.eg.db) # use keytypes to list all supporting types
bp <- readRDS("../../Stats/enrichGO.Rds")
bp.result <- bp@compareClusterResult
dotplot(bp, x=~direction) + ggplot2::facet_grid(~group)

#save original GO
png(filename="../../Figures/enrichGO.png", #print graph
    units="cm", width=24, height=18, 
    pointsize=12, res=600)
dotplot(bp, x=~direction, showCategory = 5, font.size=9) + ggplot2::facet_grid(~group)
dev.off()

#manually drop redundant GO to simplify
bp2 <- dropGO(bp, level=NULL, term=c('GO:0009167', 'GO:0009126', 'GO:0009205', 'GO:0009199', #purine metabolism
                                     'GO:0000398', 'GO:0000377', 'GO:0000375', #RNA splicing
                                     'GO:0033141', #STAT
                                     'GO:0003206', 'GO:0003007',#cardiac
                                     'GO:0000070', 'GO:0140014', #mitotic
                                     'GO:0060541') #respiratory
              )
dotplot(bp2, x=~direction, show=3) + ggplot2::facet_grid(~group)

#save cleaned GO
png(filename="../../Figures/enrichGO_simplified.png", #print graph
    units="cm", width=24, height=12, 
    pointsize=12, res=600)
dotplot(bp2, x=~direction, showCategory = 3, font.size=9) + ggplot2::facet_grid(~group)
dev.off()

